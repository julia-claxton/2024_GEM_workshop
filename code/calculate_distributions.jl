println("\n========= calculate_distributions.jl START =========\n")
println("Compiling libraries... ")
cd(@__DIR__)
include("./Julia_ELFIN_Tools/Events.jl") 
include("./Julia_ELFIN_Tools/Visualization.jl") 
include("./Julia_ELFIN_Tools/Event_Detection.jl")
using Statistics
using LinearAlgebra
using HDF5

#=
calculate_distributions.jl

This script calculates flux recorded by ELFIN and surveillance time in an L/MLT grid. The results
are saved in <project top level>/data/<distribution name>_distributions.hdf5. For the GEM workshop,
this is used to get the coverage of ELFIN in L/MLT, which is then converted to lat/long for normalization
of the atmospheric energy input distribution.

Written by Julia Claxton (julia.claxton@colorado.edu)
Released under MIT License (see <project top level>/LICENSE.txt for full license)
=#

###################################
# Functions                       #
###################################
function elfin_all_flux(; MLT_nbins = 24 * 6,
                          L_nbins = 25 * 6,
                          remove_noise = true,
                          pitch_angle_range = "full",
                          dst_range = (-Inf, Inf),
                          result_file_prefix = ""
                       )
    # Status message
    println("Creating Dataset: $result_file_prefix")
    println("\tdst_range = $(dst_range)")
    println("\tpitch_angle_range = \"$(pitch_angle_range)\"")

    # Put args in dict for clarity
    args = Dict("remove_noise" => remove_noise, 
                "pitch_angle_range" => pitch_angle_range,
                "dst_range" => dst_range)

    # Write string with detection parameters
    detection_parameters = """
                            MLT_nbins = $(MLT_nbins),
                            L_nbins = $(L_nbins),
                            remove_noise = $(remove_noise),
                            pitch_angle_range = "$(pitch_angle_range)",
                            dst_range = $(dst_range)
                            """
    detection_parameters_path = dirname(@__DIR__) * "/results/" * result_file_prefix * "_detection_parameters.txt"
    rm(detection_parameters_path, force = true)
    touch(detection_parameters_path)
    write(open(detection_parameters_path, "w"), detection_parameters)

    # Get all ELFIN data files in the data directory
    data_path = dirname(@__DIR__) * "/data/processed_scientific_data"
    days = []
    sats = []
    for filename in readdir(data_path)
        if contains(filename, ".npz")
            append!(days, [Date(filename[1:8], dateformat"yyyymmdd")])
            append!(sats, [string(filename[end-4])])
        end
    end
    
    # Set up result arrays
    MLT_bin_edges = LinRange(0, 24, MLT_nbins + 1)
    L_bin_edges = LinRange(0, 25, L_nbins + 1)
    
    coverage = zeros(L_nbins, MLT_nbins)
    fluence = zeros(16, L_nbins, MLT_nbins)
    

    #############################################################
    #                  CALCULATE DISTRIBUTIONS                  #
    #############################################################
    N = length(days)
    for i = 1:N
        print("\r\t$(i)/$(N)")
        fullday_event = create_event(days[i], sats[i])
        if fullday_event == nothing; continue; end

        for obs = 1:fullday_event.n_observations
            start_idx = fullday_event.observation_edge_idxs[obs]
            stop_idx = fullday_event.observation_edge_idxs[obs+1]-1
            obs_slice = start_idx:stop_idx

            bin_observation_period!(fullday_event, obs_slice, coverage, fluence, MLT_bin_edges, L_bin_edges, args)
        end

        if i % 30 == 0; GC.gc(); end # Get rid of ghosts in the system (undeleted Event objects) taking up loads of memeory without bound
    end

    # Get rid of MeV^-1 dimension in fluence data
    e = create_event(Date("2020-01-18"), "a") # To get energy bin values
    energy_bin_widths_MeV = (e.energy_bins_max .- e.energy_bins_min) ./ 1000
    for i = 1:16
        fluence[i,:,:] .*= energy_bin_widths_MeV[i]
    end
    
    ##############################################################
    #                     SAVE DISTRIBUTIONS                     #
    ##############################################################
    save_results(coverage, fluence, MLT_bin_edges, L_bin_edges, result_file_prefix)

    return coverage, fluence
end

function bin_observation_period!(fullday_event, obs_slice, coverage, fluence, MLT_bin_edges, L_bin_edges, args)
    # Make sure data isn't broken
    if check_for_bad_events(fullday_event, obs_slice) == true
        return
    end
    
    # Slice data from full day into just the event we're looking at
    time = fullday_event.time[obs_slice]
    MLT = fullday_event.MLT[obs_slice]
    L = fullday_event.L[obs_slice]
    dst = fullday_event.dst[obs_slice]
    fullday_e_flux, fullday_n_flux = integrate_flux(fullday_event, pitch_angle = true, pitch_angle_range = args["pitch_angle_range"])
    e_flux = fullday_e_flux[obs_slice, :]
    n_flux = fullday_n_flux[obs_slice, :] 

    if args["remove_noise"] == true
        # Thresholds determined from thresholding.jl
        e_thresholds = 10 .^[8.2, 7.56, 7.56, 7.78, 7.78, 7.78, 7.78, 7.56, 7.56, 7.56, 7.56, 7.78, 7.99, 7.99, 8.2, 8.41]
        e_mask_to_kill = [e_flux[t, channel] .< e_thresholds[channel] for t = 1:length(obs_slice), channel = 1:16]
        e_flux[e_mask_to_kill] .= 0

        n_thresholds = 10 .^[5.61, 5.61, 5.38, 5.38, 5.38, 5.15, 5.15, 4.7, 4.7, 4.47, 4.47, 4.47, 4.47, 4.47, 4.47, 4.7]
        n_mask_to_kill = [n_flux[t, channel] .< n_thresholds[channel] for t = 1:length(obs_slice), channel = 1:16]
        n_flux[e_mask_to_kill] .= 0
    end

    # Get fluence of each datapoint for the time period of interest
    datapoint_lengths = diff(time)
    append!(datapoint_lengths, 2)
    obs_period_fluence = [datapoint_lengths[t] .* n_flux[t, e] for t = 1:length(time), e = 1:16] # Rectangular integration

    # Add each datapoint's time and fluence to the result arrays
    for t_obs = 1:length(obs_slice)
        # Check if this datapoint is within the user-specified dst range
        if (args["dst_range"][1] <= dst[t_obs] < args["dst_range"][2]) == false; continue; end

        # Check if time is sorted
        if datapoint_lengths[t_obs] < 0; continue; end

        # Get fluence at this datapoint
        datapoint_fluence = obs_period_fluence[t_obs,:]
        if isnan(datapoint_fluence[1]); continue; end # Skip over if we are in loss cone or anti loss cone mode and there is no coverage of the cone

        # Get MLT and L bin idxs of this datapoint
        datapoint_MLT = MLT[t_obs]
        datapoint_L = L[t_obs]
        MLT_idx = findlast(MLT_bin_edges[1:end-1] .< datapoint_MLT)
        L_idx = findlast(L_bin_edges[1:end-1] .< datapoint_L)

        # Guard against NaN L or MLT values that result in no bin placement
        if (L_idx == nothing) || (MLT_idx == nothing); continue; end

        # Store data to results array
        fluence[:, L_idx, MLT_idx] .+= datapoint_fluence
        coverage[L_idx, MLT_idx] += datapoint_lengths[t_obs]
    end
end

function save_results(coverage_distribution, fluence_distribution, MLT_bin_edges, L_bin_edges, result_file_prefix)
    # Create event to store the energy bin values. Which particular event it is is not important.
    e = create_event(Date("2020-01-18"), "a")

    # Normalize fluence distribution to coverage
    normalized_flux_distribution = [fluence_distribution[i,:,:] ./ coverage_distribution for i = 1:16]
    normalized_flux_distribution = permutedims(cat(normalized_flux_distribution[1:end]..., dims = 3), [3, 1, 2]) # Flatten back to 3D matrix

    # Create normalized MLT and L distributions
    MLT_distro = [dropdims(sum(fluence_distribution[i,:,:], dims = 1), dims = 1) ./ dropdims(sum(coverage_distribution, dims = 1), dims = 1) for i = 1:16]
    MLT_distro = permutedims(cat(MLT_distro[1:end]..., dims = 2), [2, 1])

    L_distro = [dropdims(sum(fluence_distribution[i,:,:], dims = 2), dims = 2) ./ dropdims(sum(coverage_distribution, dims = 2), dims = 2) for i = 1:16]
    L_distro = permutedims(cat(L_distro[1:end]..., dims = 2), [2, 1])

    # Save results of analysis to file
    results_path = dirname(@__DIR__) * "/results/" * result_file_prefix * "_distributions.hdf5"
    rm(results_path, force = true)
    results_file = h5open(results_path, "w")
        results_file["/data/coverage_distribution"] = coverage_distribution
        write_attribute(results_file["/data/coverage_distribution"], "dimensions", "(1)L-shell, (2)MLT")
        write_attribute(results_file["/data/coverage_distribution"], "units", "seconds")

        results_file["/data/fluence_distribution"] = fluence_distribution
        write_attribute(results_file["/data/fluence_distribution"], "dimensions", "(1)Energy channel, (2)L-shell, (3)MLT")
        write_attribute(results_file["/data/fluence_distribution"], "units", "electrons/(cm^2 MeV s)")
        
        results_file["/data/normalized_flux_distribution"] = normalized_flux_distribution
        write_attribute(results_file["/data/normalized_flux_distribution"], "dimensions", "(1)Energy channel, (2)L-shell, (3)MLT")
        write_attribute(results_file["/data/normalized_flux_distribution"], "units", "electrons/(cm^2 MeV s)")
        
        results_file["/data/normalized_mlt_distribution"] = MLT_distro
        write_attribute(results_file["/data/normalized_mlt_distribution"], "dimensions", "(1)Energy channel, (2)MLT")
        write_attribute(results_file["/data/normalized_mlt_distribution"], "units", "electrons/(cm^2 s)")

        results_file["/data/normalized_l_shell_distribution"] = L_distro
        write_attribute(results_file["/data/normalized_l_shell_distribution"], "dimensions", "(1)Energy channel, (2)L-shell")
        write_attribute(results_file["/data/normalized_l_shell_distribution"], "units", "electrons/(cm^2 s)")

        results_file["/labels/energy_bins_min_keV"] = e.energy_bins_min
        write_attribute(results_file["/labels/energy_bins_min_keV"], "units", "keV")

        results_file["/labels/energy_bins_mean_keV"] = e.energy_bins_mean
        write_attribute(results_file["/labels/energy_bins_mean_keV"], "units", "keV")

        results_file["/labels/energy_bins_max_keV"] = e.energy_bins_max
        write_attribute(results_file["/labels/energy_bins_max_keV"], "units", "keV")

        results_file["/labels/mlt_bin_edges"] = collect(MLT_bin_edges)
        write_attribute(results_file["/labels/mlt_bin_edges"], "units", "hour")

        results_file["/labels/mlt_bin_means"] = [mean([MLT_bin_edges[i], MLT_bin_edges[i+1]]) for i = 1:(length(MLT_bin_edges)-1)]
        write_attribute(results_file["/labels/mlt_bin_means"], "units", "hour")

        results_file["/labels/l_shell_bin_edges"] = collect(L_bin_edges)
        write_attribute(results_file["/labels/l_shell_bin_edges"], "units", "Re")

        results_file["/labels/l_shell_bin_means"] = [mean([L_bin_edges[i], L_bin_edges[i+1]]) for i = 1:(length(L_bin_edges)-1)]
        write_attribute(results_file["/labels/l_shell_bin_means"], "units", "Re")

    close(results_file)
end

###################################
# Script                          #
###################################
println("\nCalculating Distributions...")
elfin_all_flux(MLT_nbins = 24,
                L_nbins = 25,
                result_file_prefix = "lowres",
                pitch_angle_range = "full",
                dst_range = (-Inf, Inf)
              )
println("\nCalculation Complete")
println("\n========== calculate_distributions.jl END ==========")