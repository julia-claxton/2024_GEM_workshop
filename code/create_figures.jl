cd(@__DIR__)
println("Compiling libraries...")
include("./Julia_ELFIN_Tools/Events.jl")
include("./Julia_ELFIN_Tools/Visualization.jl")
include("./Julia_ELFIN_Tools/Event_Detection.jl")
using Statistics
using LinearAlgebra
using Dates
using DelimitedFiles
using NPZ
using GeoDatasets
using Glob

#=
create_figures.jl

This script creates the figures for the 2024 GEM workshop poster. Each figure can be 
toggled using the booleans in the section below. Note that this code is not optimized
to create all figures at once – some calculations may be repeated when generating all
figures, so be patient.

Written by Julia Claxton (julia.claxton@colorado.edu)
Released under MIT License (see <project top level>/LICENSE.txt for full license)
=#

##############################################
# What plots to create                       #
##############################################
ELFIN_DATA        = true
EMIC_DETECTION    = true
EIP_DISTROS       = true
SUNSET_CONTOURS   = true
J_OVER_J90_DISTRO = true

BG_COLOR = RGB(0x252137) # Poster background color

##############################################
# Get EMIC Detections                        #
##############################################
data = readdlm("../results/EMIC_detections.csv", ',')
start_times = DateTime.(data[2:end, 1])
stop_times = DateTime.(data[2:end, 2])
sat = data[2:end, 3]
MLT = data[2:end, 4]
L = data[2:end, 5]
n_emic = length(start_times)


##############################################
# ELFIN Data Visualization                   #
##############################################
if ELFIN_DATA == true
    example_event = create_event(DateTime("2021-03-06T07:02:00"), DateTime("2021-03-06T07:08:00"), "b")

    tick_idxs = Int.(round.(LinRange(1, example_event.n_datapoints, 6)))
    ticks = Dates.format.(example_event.time_datetime[tick_idxs], "HH:MM:SS")

    energy = example_event.energy_bins_mean ./ 1000 # MeV

    heatmap(example_event.time, log10.(energy), log10.(example_event.Jprec_over_Jperp'),
        title = "Jprec/Jtrap",

        xlabel = "Time (UTC), $(Date(example_event.time_datetime[1]))",
        xticks = (example_event.time[tick_idxs], ticks),

        ylabel = "Energy (MeV)",
        yticks = (log10.(energy[1:3:16]), round.(energy[1:3:16], sigdigits = 2)),

        colorbar_title = "\nLog10 Jprec/Jperp",
        colormap = :ice,
        clims = (-1.25, .25),

        background = :black,
        grid = false,
        margin = 5mm,
        size = (1.3, .33) .* 500,
        dpi = 700
    )
    p1 = plot!(background = BG_COLOR)

    e_flux, _ = integrate_flux(example_event, pitch_angle = true)
    heatmap(example_event.time, log10.(energy), log10.(e_flux'),
        title = "Pitch Angle Integrated Flux",

        xlabel = "Time (UTC), $(Date(example_event.time_datetime[1]))",
        xticks = (example_event.time[tick_idxs], ticks),

        ylabel = "Energy (MeV)",
        yticks = (log10.(energy[1:3:16]), round.(energy[1:3:16], sigdigits = 2)),

        colorbar_title = "Log10 Energy Flux\n(keV/(cm² MeV s))",
        colormap = :ice,
        clims = (5, 11.5),

        background = :black,
        grid = false,
        margin = 5mm,
        size = (1.3, .33) .* 500,
        dpi = 700
    )
    p2 = plot!(background = BG_COLOR)


    function pad_image(event::Event)
        time = event.time
        x_label = "Time (UTC), $(Date(event.time_datetime[1]))"

        e_flux, n_flux = integrate_flux(event, energy = true)
        flux = log10.(e_flux) # Flux to plot
    
        # Get pitch angle bin edges
        distance_between_bins = diff(event.pitch_angles, dims = 2)
        pitch_angle_bin_edges = event.pitch_angles[:, 1:end-1] .+ (distance_between_bins ./ 2)
        # Add first and last edge
        pitch_angle_bin_edges = cat(event.pitch_angles[:,1] .- (distance_between_bins[:,1] ./ 2), pitch_angle_bin_edges, dims = 2)
        pitch_angle_bin_edges = cat(pitch_angle_bin_edges, event.pitch_angles[:,end] .+ (distance_between_bins[:,end] ./ 2), dims = 2)


        # Create grid for the image
        deltas = diff(time)
        append!(time, time[end] + deltas[end]) # Add another index so that the time grid represents edges
        image_grid_time = LinRange(0, time[end], length(time) * 2)
        image_grid_pitch_angle = LinRange(0, 180, 180 * 2)
        image = fill(NaN, length(image_grid_time), length(image_grid_pitch_angle))
    
        # Fill the image with the flux data
        for t = 1:length(time)-1
            for α = 1:16
                image_t_slice = time[t] .<= image_grid_time .< time[t+1]
                image_pa_slice = pitch_angle_bin_edges[t, α] .<= image_grid_pitch_angle .< pitch_angle_bin_edges[t, α+1]
    
                image[image_t_slice, image_pa_slice] .= flux[t, α]
            end
        end
    
        # Plot image
        replace!(image, -Inf => -100) # Makes -Inf plot as black instead of not plotting at all. Lets us differentiate between no flux vs. no measurements
        plot(
            title = "Pitch Angle Distribution",
            xlabel = x_label,
            ylabel = "Pitch Angle (deg)",
            colorbar_title = "Log10 Energy Flux\nkev/(s cm² str)",
    
            xlims = (time[1], time[end-1]),
            ylims = (0, 180),
    
            framestyle = :box,
            background = :black,
            grid = false,
            markershape = :square,
            markersize = (distance_between_bins ./ 5) .+ 1,
            markerstrokewidth = 0,
    
            aspect_ratio = ((time[end-1] - time[1]) / 180) * .15,
            leftmargin = 8mm,
            size = (1, .3) .* 600,
            dpi = 700
        )
        heatmap!(image_grid_time, image_grid_pitch_angle, image',
            colormap = :ice,
            clims = (4.0, 8.5)
        )
        tick_idxs = Int.(round.(LinRange(1, event.n_datapoints, 6)))
        ticks = Dates.format.(event.time_datetime[tick_idxs], "HH:MM:SS")
        plot!(xticks = (event.time[tick_idxs], ticks))
        deleteat!(time, length(time)) # Remove extra time index we added before for the heatmap plot to work
    
        # Plot loss cones
        plot!(time, event.loss_cone_angles,
            linewidth = 2,
            linecolor = :white,
            linestyle = :solid,
            label = ""
        )
        plot!(time, event.anti_loss_cone_angles,
            linewidth = 2,
            linecolor = :white,
            linestyle = :dash,
            label = ""
        )
        return plot!(background = BG_COLOR)
    end
    p3 = pad_image(example_event)


    layout = @layout [ plot1
                    plot2
                    plot3
    ]

    plot(p1, p2, p3,
        layout = layout,
        background = BG_COLOR,
        size = (1, .7) .* 800
    )
    png("../results/figures/ELFIN_data.png")
    display(plot!())
end

##############################################
# EMIC Detection Method                      #
##############################################
if EMIC_DETECTION == true
    function emic_detection(event::Event; time_range_idx = 1:event.n_datapoints)    
        high_e_score = high_energy_precipitation_score(event, time_range_idx = time_range_idx)
        low_e_score  = low_energy_precipitation_score(event, time_range_idx = time_range_idx)
    
        # Look at the region around the low-E score and set score to the highest in that region
        low_e_smoothing_radius = 4
        smoothed_low_e = copy(low_e_score) .* 0
        for i = 1:length(smoothed_low_e)
            lower_idx = clamp(i - low_e_smoothing_radius, 1, length(time_range_idx))
            upper_idx = clamp(i + low_e_smoothing_radius, 1, length(time_range_idx))
            smoothed_low_e[i] = max(low_e_score[lower_idx:upper_idx]...)
        end
    
        # Apply conditions for EMIC detection
        emic_mask = (high_e_score .>= .75) .&& (smoothed_low_e .< 4) .&& (smoothed_low_e .<= high_e_score .+ 1.5)
        emic_mask[end] = 0 # If the last datapoint is high, call it the end of an event. If it's low, no impact.
        event_edges = diff(emic_mask)
        if event_edges[1] == -1
            event_edges[1] = 0 # If first datapoint is high and second low, discard detection
        else
            event_edges[1] = emic_mask[1] # If the first datapoint is high, call it the start of an event
        end
    
        # Get info on detected events
        emic_start_idxs = findall(event_edges .== 1)
        emic_stop_idxs = findall(event_edges .== -1)
    
        # Plotting
        yticks = Int.(round.(event.energy_bins_mean))
        yticks = yticks[1:3:16]
        p1 = heatmap(time_range_idx, 1:16, log10.(event.Jprec_over_Jperp[time_range_idx,:]'),
            background = :black,
            ylabel = "Energy Channel",
            colorbar = false,
            clims = (-1.25, .25),
            xlims = (time_range_idx[1], time_range_idx[end]),
            ylims = (1,16),
            yticks = (1:3:16, yticks),
            aspect_ratio = ((time_range_idx[end] - time_range_idx[1])/16) * .15,
            c = :ice
        )
        # Precipitation Scores
        plot(
            title = "Precipitation Scores",
            xlabel = "Data Index",
            ylabel = "Score (arbitrary units)",
            xlims = (time_range_idx[1], time_range_idx[end]),
            ylims = (0, 4),
            aspect_ratio = ((time_range_idx[end] - time_range_idx[1])/4) * .15,
            background = :black
        )
        plot!(time_range_idx, high_e_score, label = "High-Energy Precip. Score")
        plot!(time_range_idx, smoothed_low_e, label = "Low-Energy Precip. Score")
        p2 = plot!()

        # Detection area
        plot(time_range_idx, 1 .* emic_mask,
            fill = true,
            background = :black,
            linetype = :steppre,
            fillcolor = :red,
            linecolor = :transparent,
            xlims = (time_range_idx[1], time_range_idx[end]),
            ylims = (0, 1),
            framestyle = :none,
            grid = false,
            label = false,
            topmargin = 5mm,
            bottommargin = -10mm
        )
        p3 = plot!()

        layout = @layout [a{.1h}; b; c]
        plot(p3, p1, p2,
            layout = layout,
            background = BG_COLOR,
            suptitle = "$(event.name)",
            dpi = 350
        )
    end


    example_event = create_event(DateTime("2021-03-06T07:02:00"), DateTime("2021-03-06T07:08:00"), "b")
    p1 = emic_detection(example_event)
    display(p1)
    png("../results/figures/EMIC_detection.png")
end

##############################################
# EIP Source and Atmo Distributions          #
##############################################
if EIP_DISTROS == true
    function heatmap_with_histograms(distribution, MLT_bin_edges, L_bin_edges, clims)
        ############## HEATMAP ##############
        heatmap(MLT_bin_edges, L_bin_edges, distribution,
            xlabel = "MLT (T87)",
            ylabel = "L (T87)",

            xlims = (0,24),
            xticks = 0:3:24,
            xminorticks = 3,

            ylims = (0,25),
            yticks = 0:2:25,
            yminorticks = 2,

            clims = clims,
            colorbar_ticks = :none, 
            colormap = :ice,
            colorbar = false,

            leftmargin = 2mm,
            rightmargin = -5mm,

            background = :black, # turns text white
            tickdirection = :out,
            framestyle = :box,
            grid = false
        )
        plot!(twiny(),
            xlims = (0,24),
            xticks = (0:3:24, repeat("", length(0:3:24))),
            xminorticks = 3,
            tickdirection = :out
        )
        plot!(twinx(),
            ylims = (0,25),
            yticks = (0:2:25, repeat("", length(0:2:25))),
            yminorticks = 2,
            tickdirection = :out
        )
        main = plot!()

        ############## COLORBAR ##############
        colorbar = scatter([0,0], [0,1],
            zcolor = [0,3],
            xlims = (1,1.1),
            label= "",
            framestyle = :none,
            colormap = :ice,
            colorbar_title = "\nEMIC-Driven Events/Second • 10³, #/s",
            background = :black,
            clims = clims,
            leftmargin = -15mm,
            rightmargin = 5mm,
            topmargin = 5mm
        )


        ############## MLT HISTOGRAM ##############
        distribution_finite_only = replace(distribution, NaN => 0, Inf => 0, -Inf => 0)
        MLT_frequencies = dropdims(sum(distribution_finite_only, dims = 1), dims = 1)
        MLT_frequencies = cat(MLT_frequencies, MLT_frequencies[end], dims = 1) # Make match number of bin edges
        integral = sum([(MLT_bin_edges[i+1] - MLT_bin_edges[i]) * MLT_frequencies[i] for i = 1:length(MLT_bin_edges)-1])
        MLT_frequencies ./= integral # Normalize to a pdf
        plot(MLT_bin_edges, MLT_frequencies,
            linetype = :steppost,
            linecolor = :transparent,
            fill = true,
            fillcolor = RGB(0x4073B4),

            xlims = (0, 24),
            xticks = false,

            ylims = (0, .15),
            yticks = false,

            label = false,
            grid = false,
            framestyle = :none,

            background = :black,
            bottommargin = -4mm
        )
        mlt = plot!()


        ############## L HISTOGRAM ##############
        L_frequencies = dropdims(sum(distribution_finite_only, dims = 2), dims = 2)
        L_frequencies = cat(L_frequencies, L_frequencies[end], dims = 1) # Make match number of bin edges
        integral = sum([(L_bin_edges[i+1] - L_bin_edges[i]) * L_frequencies[i] for i = 1:length(L_bin_edges)-1])
        L_density = L_frequencies ./ integral # Normalize to a pdf
        L_cdf = cumsum(L_density)
        
        # The `fill` plot option is behaving badly with the axis permutation, so I have to plot the histogram shape myself
        x = [L_bin_edges L_bin_edges]'[:]
        y = circshift([L_density L_density]'[:], -1)
        hist_shape = Shape(x, y)

        plot(hist_shape,
            permute = (:x, :y),
            linetype = :steppost,
            linecolor = :transparent,
            fill = true,
            fillcolor = RGB(0x4073B4),

            xlims = (0, 25),

            ylims = (0, .2),
            yticks = false,

            topmargin = 12mm,
            bottommargin = 9.25mm,
            leftmargin = -31mm,
            rightmargin = 0mm,

            background = :black,
            label = false,
            grid = false,
            framestyle = :none
        )
        l = plot!()


        ############## FULL PLOT ##############
        layout = @layout [ [p1{.8w}
                            p2{.9h}] p3{.1w} p4{.12w}
        ]
        plot(mlt, main, l, colorbar,
            layout = layout,
            background = BG_COLOR,
            suptitle = "Origin of EMIC-Driven Precipitation (N = $(n_emic))",
            dpi = 450
        )
        return plot!()
    end

    # Settings
    latitude_nbins = 45
    longitude_nbins = 45
    MLT_nbins = 24
    L_nbins = 25

    latitude_bin_edges = LinRange(-90, 90, latitude_nbins + 1)
    longitude_bin_edges = LinRange(-180, 180, longitude_nbins + 1)
    MLT_bin_edges = LinRange(0, 24, MLT_nbins + 1)
    L_bin_edges = LinRange(0, 25, L_nbins + 1)


    # MLT/lat coverage distribution
    # Get the time spent in every MLT/latitude bin
    println("Gathering coverage data...")
    data_path = "../data/processed_scientific_data"
    files = glob("*.npz", data_path)
    N = length(files)

    x = []
    y = []
    z = []
    coverage_MLT = []
    coverage_duration = []

    # For each day of data
    for i = 1:N
        print("\r\t$(i)/$(N)")
        date = Date(files[i][end-15:end-8], dateformat"yyyymmdd")
        satellite = string(files[i][end-4])
        event = create_event(date, satellite)
        if event == nothing; continue; end

        if isfile("../data/all_coverage_latlong.npz") == false # Only do this if needed
            append!(x, event.position["x"])
            append!(y, event.position["y"])
            append!(z, event.position["z"])
        end

        # But we always need MLT and duration
        append!(coverage_MLT, event.MLT)

        datapoint_lengths = diff(event.time)
        datapoint_lengths[datapoint_lengths .>= 30] .= 3 # Observation edges
        append!(datapoint_lengths, 3) # Last datapoint

        append!(coverage_duration, datapoint_lengths)

        if i % 30 == 0; GC.gc(); end
    end
    if isfile("../data/all_coverage_latlong.npz") == false
        # Convert GEI to lat/long
        println("Converting coverage to MLT/latitude...")
        npzwrite("../data/TEMP_locations_gei.npy", hcat(x,y,z))
        println()
        run(`python3.9 ./convert_gei_to_lla.py`, wait = true)
        mv("../data/locations_latlong.npz", "../data/all_coverage_latlong.npz") # Rename file so we can keep it
        x = y = z = []; GC.gc() # Clear up memory
        println("Coverage coordinate conversion complete\n")
    end
    coverage_duration = float.(coverage_duration)
    converted_coverage_coords = npzread("../data/all_coverage_latlong.npz")
    coverage_latitude = converted_coverage_coords["lat"]
    coverage_longitude = converted_coverage_coords["lon"]

    # Distribution of EMIC energy in MLT/latittude
    # Calculate energy input of each event, as well as getting ELFIN's position for field tracing
    energy_input = zeros(n_emic)
    durations = zeros(n_emic)
    positions = zeros(n_emic, 3)
    for i = 1:n_emic
        # Be careful here with loss cone correction. This will need to be dealt with more rigorously in the future.
        print("\r\t$(i)/$(n_emic)")
        event = create_event(start_times[i], stop_times[i], sat[i])
        energy_input[i], _ = integrate_flux(event, time = true, energy = true, pitch_angle = true, pitch_angle_range = "loss cone") # Units = keV/cm^2
        durations[i] = event.duration
        positions[i,1] = event.position["x"][1]
        positions[i,2] = event.position["y"][1]
        positions[i,3] = event.position["z"][1]
    end
    GC.gc()

    # Convert GEI to lat/long
    npzwrite("../data/TEMP_locations_gei.npy", positions)
    run(`python3.9 ./convert_gei_to_lla.py`, wait = true)
    converted_coords = npzread("../data/locations_latlong.npz")
    emic_latitude = converted_coords["lat"]
    emic_longitude = converted_coords["lon"]
    rm("../data/locations_latlong.npz")
    println("EMIC coordinate conversion complete\n")
    

    # Calculate distributions
    # Atmospheric input figure (MLT/lat)
    energy_input_distro = my_2d_histogram(emic_longitude, emic_latitude, longitude_bin_edges, latitude_bin_edges, weights = energy_input)
    coverage_distro = my_2d_histogram(coverage_longitude, coverage_latitude, longitude_bin_edges, latitude_bin_edges, weights = coverage_duration)
    emic_flux_input_distro = energy_input_distro ./ coverage_distro
    replace!(emic_flux_input_distro, 0 => NaN)

    # EMIC origin figure (MLT/L)
    # Get coverage for normalization
    coverage_path = dirname(@__DIR__) * "/results/lowres_distributions.hdf5"
    results = h5open(coverage_path, "r")
    coverage_distribution = read(results["data/coverage_distribution"])
    coverage_distribution[coverage_distribution .< 300] .= 0 # To avoid artificially blowing up the ratio

    # Get EMIC distro
    emics = my_2d_histogram(L, MLT, L_bin_edges, MLT_bin_edges)
    normalized_emic_distribution = emics ./ coverage_distribution
    
    # Get coastlines for plotting
    coastlines_lon, coastlines_lat, coastlines = GeoDatasets.landseamask(;resolution='c',grid=5)
    coastlines_MLT = ((collect(coastlines_lon) ./ 360) .+ .5) .* 24

    # Plot results
    # Plot the world overlaid with energy input
    heatmap(longitude_bin_edges, latitude_bin_edges, log10.(emic_flux_input_distro'),
        title = "EMIC-Driven Atmospheric Energy Input",

        xlabel = "Longitude, deg",
        xlims = (-180, 180),
        xticks = -180:60:180,
        xminorticks = 6,

        ylabel = "Latitude, deg",
        ylims = (-90, 90),
        yticks = -90:30:90,
        yminorticks = 3,

        clims = (2,7),

        background = :black, # turn text white
        tickdirection = :out,
        colormap = :ice,
        colorbar = false,

        grid = true,
        minorgrid = true,
        foreground_color_grid = :black,
        foreground_color_minor_grid = :black,
        aspect_ratio = 1.2,
    )
    plot!(
        background_color_inside = RGB(0xd4d2dc),
        background_color_outside = BG_COLOR
    )
    plot!(twinx(),
        xlims = (-180, 180),
        xticks = -180:60:180,
        xminorticks = 6,

        ylims = (-90, 90),
        yticks = -90:30:90,
        yminorticks = 3,

        aspect_ratio = 1.2,
        grid = true,
        foreground_color_grid = :black,
        tickdirection = :out
    )
    plot!(twiny(),
        xlims = (-180, 180),
        xticks = -180:60:180,
        xminorticks = 6,

        ylims = (-90, 90),
        yticks = -90:30:90,
        yminorticks = 3,

        aspect_ratio = 1.2,
        grid = true,
        foreground_color_grid = :black,
        tickdirection = :out
    )
    to_plot_coastlines = replace(coastlines', 2 => 0) .* 8
    contour!(coastlines_lon, coastlines_lat, to_plot_coastlines,
        levels = 5:11,
        clims = (2,7),
        linecolor = :black
    )
    main = plot!()
    # Plot the colorbar
    colorbar = scatter([0,0], [0,1],
        zcolor = [0,3],
        xlims = (1,1.1),
        label= "",
        framestyle = :none,
        colormap = :ice,
        colorbar_title = "Log10 Average EMIC-Driven Flux, keV/(cm² s)",
        clims = (2,7),
        leftmargin = -15mm,
        topmargin = 0mm
    )
    layout = @layout [a b{.15w}]
    plot(main, colorbar,
        layout = layout,
        background = BG_COLOR,
        dpi = 350
    )
    png("../results/figures/emic_atmosphere_input.png")
    display(plot!())


    # Plot where EMIC events originate
    heatmap_with_histograms(normalized_emic_distribution .* 10^3, MLT_bin_edges, L_bin_edges, (0, 1))
    png("../results/figures/emic_origin_cartesian.png")
    display(plot!())


    # Polar plot
    L_max = 15
    L_final_idx = findlast(L_bin_edges .< L_max) # Because polar doesn't trim to a circle, it crops a square
    heatmap(MLT_bin_edges .* (2π/24), L_bin_edges[1:L_final_idx+1], normalized_emic_distribution[1:L_final_idx, :] .* 10^3,
        proj = :polar,
        axis = false,
        grid = false,
        background = :black,

        ylims = (0, L_max),

        colorbar = false,
        colormap = :ice,
        clims = (0, 1),

        rightmargin = 10mm,
        framestyle = :none,
        dpi = 350
    )
    plot!(
        background = BG_COLOR
    )
    # L labels
    plot_ring!.(5:5:L_max)
    annotate!(0, (5/L_max) - .09, text("5", 12, RGBA(1,1,1,.9)))
    annotate!(0, (10/L_max) - .09, text("10", 12, RGBA(1,1,1,.9)))
    annotate!(0, (15/L_max) - .09, text("15", 12, RGBA(1,1,1,.9)))
    # MLT labels
    annotate!(1.2, 0, text("0000", 12, :white))
    annotate!(-1.2, 0, text("1200", 12, :white))
    annotate!(0, 1.1, text("0600", 12, :white))
    annotate!(0, -1.1, text("1800", 12, :white))
    # MLT lines
    mlt_line!(MLT) = plot!([MLT, MLT] .* (2π/24), [0, L_max],
        linecolor = RGBA(.5,.5,.5,.9),
        label = false
    )
    mlt_line!.(setdiff(0:3:24, [6]))
    # Earth
    plot!(Shape(vcat(Plots.partialcircle(0, 2π, 100, 0), reverse(Plots.partialcircle(0, 2π, 100, 1/L_max)))),
        fillcolor = :white,
        linecolor = RGBA(0,0,0,0),
        label = ""
    )
    plot!(Shape(vcat(Plots.partialcircle(-π/2, π/2, 100, 0), reverse(Plots.partialcircle(-π/2, π/2, 100, 1/(L_max + 2) )))), # 1/27 makes it smaller to leave a white border around the earth
        fillcolor = :black,
        linecolor = RGBA(0,0,0,0),
        label = ""
    )
    png("../results/figures/emic_origin_polar.png")
    display(plot!())
end

##############################################
# EMIC Sunset Contours                       #
##############################################
if SUNSET_CONTOURS == true
    function get_all_sunset()
        pa_nbins = 16

        # Create an event to get ELFIN's energy bins
        energy_event = create_event(start_times[1], stop_times[1], sat[1])
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

        for i = 1:n_emic
            print("\r$(i)/$(n_emic)")
            event = create_event(start_times[i], stop_times[i], sat[i])
            _, n_fluence = integrate_flux(event, time = true)

            # Reverse flux readings in SH to make loss cones match
            α_lc = event.avg_loss_cone_angle
            α_alc = event.avg_anti_loss_cone_angle
            if event.avg_loss_cone_angle > 90
                reverse!(n_fluence, dims = 2)
                α_lc = 180 - α_lc
                α_alc = 180 - α_alc
            end
            append!(loss_cone_angles, α_lc)
            append!(anti_loss_cone_angles, α_alc)

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
            all_fluence .+= my_2d_histogram(E, α, E_bin_edges, pa_bin_edges, weights = fluence)
            all_time += event.duration
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

        return pa_bin_edges, E_bin_edges, α_lc_avg, α_alc_avg, J_over_J90
    end

    pa_bin_edges, E_bin_edges, α_lc_avg, α_alc_avg, J_over_J90 = get_all_sunset()
    
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
        levels = 9,
        fill = true,
        contour_labels = true,

        title = "EIP Pitch Angle Distributions",
        titlefontcolor = :white,

        xlabel = "Pitch Angle, deg",
        xlims = (0, 180),
        xticks = 0:30:180,
        xminorticks = 3,
        xguidefont = :white,
        xtickfontcolor = :white,
        x_foreground_color_axis = :white,

        ylabel = "Energy, eV",
        ylims = (E_bin_means_trimmed[1], E_bin_means_trimmed[end]),
        yscale = :log10,
        yminorticks = true,
        yguidefont = :white,
        ytickfontcolor = :white,
        y_foreground_color_axis = :white,

        colorbar_title = "\nLog10 J/J90",
        clims = (-2, 0),
        colormap = :ice,
        colorbar_titlefontcolor = :white,

        rightmargin = 5mm,

        aspect_ratio = (180/(E_bin_means_trimmed[end]-E_bin_means_trimmed[1])) * .7,
        tickdirection = :out,

        background = :black, # sets all text to be white
        grid = false,
        size = (1.5, 1) .* 400,
        dpi = 350
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

        clims = (-2, 0),
        colormap = :ice,

        background = :black, # sets all text to be white
        tickdirection = :out,
        aspect_ratio = (180/(E_bin_means_trimmed[end]-E_bin_means_trimmed[1])) * .7,
    )
    plot!(twinx(), [1, 1],
        zcolor = [1, 1],
        label = false,

        xlims = (0, 180),
        xticks = 0:30:180,
        xminorticks = 3,
        xtickfontcolor = :transparent,

        ylims = (E_bin_means_trimmed[1], E_bin_means_trimmed[end]),
        yscale = :log10,
        yminorticks = true,

        clims = (-2, 0),
        colormap = :ice,

        background = :black, # sets all text to be white
        tickdirection = :out,
        aspect_ratio = (180/(E_bin_means_trimmed[end]-E_bin_means_trimmed[1])) * .7,
    )
    plot!(
        background_color_outside = BG_COLOR,
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


    #=
    x = [22, 32, 40, 48, 51, 51, 51]
    E_min = [140, 170, 170, 180, 310, 530, 1020] 
    plot!(x, E_min,
        linewidth = 1,
        marker = true,
        markercolor = RGBA(1,0,0),
        markeralpha = .5,
        linecolor = RGBA(1,0,0,.5),
    )
    =#

    png("../results/figures/EIP_sunset_contours.png")
    display(plot!())
end

##############################################
# J/J90 Distribution                         #
##############################################
if J_OVER_J90_DISTRO == true
    function get_all_E_and_J_over_J90()
        all_E = []
        all_Jprec_over_J90 = []
        for i = 1:n_emic
            print("\r$(i)/$(n_emic)")
            event = create_event(start_times[i], stop_times[i], sat[i])
            _, n_fluence = integrate_flux(event, time = true)
            n_flux = n_fluence ./ event.duration
            J_over_J90 = copy(n_flux) .* 0

            lc_idxs = event.lc_idxs[1]

            for e = 1:16
                _, normalizing_idx = findmin(abs.(90 .- event.avg_pitch_angles))
                J90 = n_flux[e, normalizing_idx]
                J_over_J90[e,:] = n_flux[e,:] ./ J90


                append!(all_E, [event.energy_bins_mean[e] for α = lc_idxs])
                append!(all_Jprec_over_J90, [J_over_J90[e, α] for α = lc_idxs])
            end
        end
        return all_E, all_Jprec_over_J90
    end

    # Create an event to get ELFIN's energy bins
    energy_event = create_event(start_times[1], stop_times[1], sat[1])
    E_bin_means = energy_event.energy_bins_mean

    all_E, all_Jprec_over_J90 = get_all_E_and_J_over_J90()

    plot(
        title = "Precipitating Flux Distribution by \nEnergy Channel\n",

        xlabel = "Energy, keV",
        xlims = (50, 8500),
        xscale = :log10,
        xminorticks = true,

        ylabel = "Jprec/Jperp",
        ylims = (0, 1.2),
        yticks = 0:.25:1,
        yminorticks = 5,
        
        grid = true,
        minorgrid = true,
        foreground_color_grid = :white,
        gridalpha = .3,
        minorgridalpha = .15,

        background = BG_COLOR,
        size = (1.5, 1.1) .* 400,
        thickness_scaling = 1.2,
        dpi = 350
    )
    plot!(twinx(),
        xlims = (50, 8500),
        xscale = :log10,
        xminorticks = true,
        xtickfontcolor = :transparent,

        ylims = (0, 1.2),
        yticks = 0:.25:1,
        yminorticks = 5,
        ytickfontcolor = :transparent,

        grid = true,
        foreground_color_grid = :white,
        gridalpha = .3
    )
    plot!(twiny(),
        xlims = (50, 8500),
        xscale = :log10,
        xminorticks = true,
        xtickfontcolor = :transparent,

        ylims = (0, 1.2),
        yticks = 0:.25:1,
        yminorticks = 5,
        ytickfontcolor = :transparent,

        grid = true,
        foreground_color_grid = :white,
        gridalpha = .3
    )
    for e = 1:16
        channel_idxs = all_E .== energy_event.energy_bins_mean[e]
        ratios = all_Jprec_over_J90[channel_idxs]
        idxs_to_remove = findall(isnan.(ratios) .|| isinf.(ratios) .|| (ratios .== 0))
        deleteat!(ratios, idxs_to_remove)

        quantiles = quantile(ratios, [.25, .50, .75])

        plot!(repeat([E_bin_means[e]], 3), quantiles,
            label = false,

            marker = true,
            markersize = 5,
            markercolor = [RGB(0x3E4B96), :white, RGB(0x3E4B96)],
            
            linecolor = :grey,
            linewidth = 2
        )
    end
    # Plot dummy data for legend
    scatter!(twinx(), repeat([[-1]], 3),
        markercolor = [RGB(0x3E4B96) :white RGB(0x3E4B96)],
        label = ["75th Percentile" "50th Percentile" "25th Percentile"],

        xlims = (50, 8500),
        xscale = :log10,
        xminorticks = true,
        xtickfontcolor = :transparent,

        ylims = (0, 1.2),
        yticks = 0:.25:1,
        yminorticks = 5,
        ytickfontcolor = :transparent,


        legend_position = :topleft,
        background_color_legend = :black,

        axis = false,
        grid = false
    )

    png("../results/figures/energy_profile.png")
    display(plot!())
end