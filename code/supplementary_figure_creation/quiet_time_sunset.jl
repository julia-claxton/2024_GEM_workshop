println("\n========= analyze.jl START =========\n")
println("Compiling libraries... ")
cd(@__DIR__)
include("../Julia_ELFIN_Tools/Events.jl") 
include("../Julia_ELFIN_Tools/Visualization.jl") 
include("../Julia_ELFIN_Tools/Event_Detection.jl")
using Statistics
using LinearAlgebra
using Dates
using Plots
using Glob
using DelimitedFiles

#=
[name].jl


this is bad code




Written by Julia Claxton (julia.claxton@colorado.edu)
Released under MIT License (see <project top level>/LICENSE.txt for full license)
=#


function get_fluence_histogram(E_bin_edges, pa_bin_edges, event::Event)
    _, n_fluence = integrate_flux(event, time = true)

    # Reverse flux readings in SH to make loss cones match
    if event.avg_loss_cone_angle > 90
        reverse!(n_fluence, dims = 2)
    end

    # Put data into proper format for 2D histogram making and then make the histogram
    α = []
    E = []
    fluence = []
    for pa_idx = 1:16
        for E_idx = 1:16
            append!(α, event.avg_pitch_angles[pa_idx])
            append!(E, event.energy_bins_mean[E_idx])
            append!(fluence, n_fluence[E_idx, pa_idx])
        end
    end
    α = float.(α)
    E = float.(E)
    fluence = float.(fluence)

    fluence_histogram = my_2d_histogram(E, α, E_bin_edges, pa_bin_edges, weights = fluence)
    return fluence_histogram
end

function get_quiet_J_over_J90()
    # Find quiet (no precipitation) events
    data_path = "../../data/processed_scientific_data"
    files = glob("*.npz", data_path)
    n_days = length(files)

    pa_nbins = 16

    # Create an event to get ELFIN's energy bins
    energy_event = create_event(DateTime("2020-01-08T05:44:00"), DateTime("2020-01-08T06:00:00"), "a")
    E_bin_edges = cat(energy_event.energy_bins_min, energy_event.energy_bins_max[end], dims = 1)
    E_nbins = length(E_bin_edges) - 1

    # Create pitch angle bins
    pa_bin_edges = LinRange(0, 180, pa_nbins + 1)
    pa_bin_means = [mean([pa_bin_edges[i], pa_bin_edges[i+1]]) for i = 1:pa_nbins]

    # Create results arrays
    all_fluence = zeros(E_nbins, pa_nbins)
    all_time = 0
    loss_cone_angles = []
    anti_loss_cone_angles = []

    # We pick random days because holding all 1500 days of data in memory is too much
    sample_size = n_days
    days_to_get_data_from = sample(1:n_days, sample_size, replace = false)

    for i = 1:sample_size
        day_number = days_to_get_data_from[i]
        print("\r\t$(i)/$(sample_size)")
        date = Date(files[day_number][end-15:end-8], dateformat"yyyymmdd")
        sat = string(files[day_number][end-4])
        event = create_event(date, sat)
        if event == nothing; continue; end

        # For each observation period in the day
        for obs_number = 1:event.n_observations
            start_idx = event.observation_edge_idxs[obs_number]
            stop_idx = event.observation_edge_idxs[obs_number+1]-1
            slice = start_idx:stop_idx

            quiet_start_idxs, quiet_stop_idxs = quiet_detector(event, time_range_idx = slice, show_plot = false)
            n_events = length(quiet_start_idxs)
            
            for detection_number = 1:n_events
                detection_event = create_event(event.time_datetime[quiet_start_idxs[detection_number]], event.time_datetime[quiet_stop_idxs[detection_number]], event.satellite)
                
                # Flux information
                all_fluence .+= get_fluence_histogram(E_bin_edges, pa_bin_edges, detection_event)
                all_time += detection_event.duration
            
                # Loss cone information
                α_lc = detection_event.avg_loss_cone_angle
                α_alc = detection_event.avg_anti_loss_cone_angle    
                if detection_event.avg_loss_cone_angle > 90
                    α_lc = 180 - α_lc
                    α_alc = 180 - α_alc
                end
                append!(loss_cone_angles, α_lc)
                append!(anti_loss_cone_angles, α_alc)
            end
        end
        if day_number % 5 == 0; GC.gc(); end
    end
    println()


    # Get flux and normalize it to J90
    all_flux = all_fluence ./ all_time
    J_over_J90 = copy(all_flux) .* 0

    _, normalizing_idx = findmin(abs.(90 .- pa_bin_means))
    for e = 1:16 # Normalize each energy bin (row)
        J90 = all_flux[e, normalizing_idx]
        J_over_J90[e,:] = all_flux[e,:] ./ J90
    end

    # Get average LC and ALC angles
    α_lc_avg = mean(loss_cone_angles)
    α_alc_avg = mean(anti_loss_cone_angles)

    # Return
    return pa_bin_edges, E_bin_edges, α_lc_avg, α_alc_avg, J_over_J90
end


pa_bin_edges, E_bin_edges, α_lc_avg, α_alc_avg, J_over_J90 = get_quiet_J_over_J90()





















# Get nbins and bin means
pa_nbins = length(pa_bin_edges) - 1
E_nbins = length(E_bin_edges) - 1
pa_bin_means = [mean([pa_bin_edges[i], pa_bin_edges[i+1]]) for i = 1:pa_nbins]
E_bin_means = [mean([E_bin_edges[i], E_bin_edges[i+1]]) for i = 1:E_nbins]

# Trim off the top two energy bins because they're not very useful
J_over_J90_trimmed = J_over_J90[1:end-2,:]
E_bin_means_trimmed = E_bin_means[1:end-2]
E_bin_edges_trimmed = E_bin_edges[1:end-2]

# Plot
contour(pa_bin_means, E_bin_means_trimmed, log10.(J_over_J90_trimmed),
    title = "Quiet Pitch Angle Distributions",
    levels = 7,
    fill = true,
    contour_labels = true,

    xlabel = "Pitch Angle, deg",
    xlims = (0, 180),
    xticks = 0:30:180,
    xminorticks = 3,

    ylabel = "Energy, eV",
    ylims = (E_bin_means_trimmed[1], E_bin_means_trimmed[end]),
    yscale = :log10,
    yminorticks = true,

    colorbar_title = "\nLog10 J/J90",
    clims = (-2, 0),
    colormap = :ice,

    rightmargin = 5mm,

    aspect_ratio = (180/(E_bin_means_trimmed[end]-E_bin_means_trimmed[1])) * .7,
    tickdirection = :out,

    grid = false,
    background_color_inside = :black,
    size = (1.5, 1) .* 400,
    dpi = 350
)
plot!(twinx(), [1, 1],
    zcolor = [1, 1],
    label = false,

    xlims = (0, 180),
    xticks = 0:30:180,
    xminorticks = 3,

    ylims = (E_bin_means_trimmed[1], E_bin_means_trimmed[end]),
    yscale = :log10,
    yminorticks = true,
    ytickfontcolor = :transparent,

    clims = (-2, 0),
    colormap = :ice,

    tickdirection = :out,
    aspect_ratio = (180/(E_bin_means_trimmed[end]-E_bin_means_trimmed[1])) * .7,
)
plot!(twiny(), [1, 1],
    zcolor = [1, 1],
    label = false,

    xlims = (0, 180),
    xticks = 0:30:180,
    xminorticks = 3,
    xtickfontcolor = :transparent,

    ylims = (E_bin_means_trimmed[1], E_bin_means_trimmed[end]),
    yscale = :log10,
    yminorticks = true,
    ytickfontcolor = :transparent,

    clims = (-2, 0),
    colormap = :ice,

    tickdirection = :out,
    aspect_ratio = (180/(E_bin_means_trimmed[end]-E_bin_means_trimmed[1])) * .7,
)
vline!([α_lc_avg],
    linecolor = :red,
    linewidth = 2,
    linestyle = :solid,
    label = "Loss Cone"
)
vline!([α_alc_avg],
    linecolor = :red,
    linewidth = 2,
    linestyle = :dash,
    label = "Anti Loss Cone",
    legendposition = :top
)
png("../../results/figures/supplementary/quiet_sunset.png")
display(plot!())














#=
header = ["start_time" "stop_time" "satellite" "MLT" "L"]
data = hcat(start_times, stop_times, satellites, MLT, L)
to_write = vcat(header, data)
writedlm("../results/quiet_detections.csv", to_write, ',')

println("\nAnalysis Complete")
println("\n========== analyze.jl END ==========")
=#