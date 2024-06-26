include("./Events.jl")
using Statistics
using LinearAlgebra
using Dates
using Plots
using Plots.PlotMeasures
using NumericalIntegration
using ColorSchemes

#=
Visualization.jl

This file provides several methods for quickly visualizing the data stored in an Event object. The
most important method is quicklook(), which takes in an Event and displays a plot containing an
overview of the event.

Written by Julia Claxton (julia.claxton@colorado.edu)
Released under MIT License (see <project top level>/LICENSE.txt for full license)
=#

function quicklook(event::Event; by = "index")
# `by` kwarg determines what the x-axis of the plots should be and can take values of "index", "time", or "date".
    if event.n_datapoints == 1; @warn "n_datapoints = 1, no plot"; return nothing; end
    
    energy_heatmap = energy_time_series(event,         by = by, show_plot = false)
    lc_ratio       = J_over_J90_time_series(event, by = by, show_plot = false)
    pad_heatmap    = pad_time_series(event,            by = by, show_plot = false)
    L_shell        = L_MLT_time_series(event,              by = by, show_plot = false)
    abs_sunset     = absolute_sunset(event,                     show_plot = false)
    rel_sunset     = relative_sunset(event,                   show_plot = false)
    mlt_l_path     = MLT_L_track(event,                         show_plot = false)

    blank = plot(yticks = nothing, framestyle = :none)

    l = @layout [ plot1{.14h}
                  plot2{.14h}
                  plot3{.14h}
                  blank plot4{.93w} blank{.1125w}
                  plot5 plot6
    ]

    plot(lc_ratio, energy_heatmap, pad_heatmap, blank, L_shell, blank, abs_sunset, rel_sunset,
        layout = l,
        plot_title = "$(Date(event.time[1])) $(Time(event.time[1])) ELFIN-$(event.satellite)\nDuration = $(event.duration) s / $(event.n_datapoints) observations",
        background_color = :transparent,
        size = (1,1.4) .* 900,
        format = :png,
        dpi = 700
    )
    display("image/png", plot!())
end

##############################################################################
#                                                                            #
#                            SUPPORTING FUNCTIONS                            #
#                                                                            #
##############################################################################

function quicklook(event::Nothing; by = "index")
    @warn "No event"
    return nothing
end

function strip_heatmap(z; colormap = cgrad(:inferno))
    dims = size(z)
    x = 1:dims[2]
    y = 1:dims[1]

    return strip_heatmap(x, y, z, colormap = colormap)
end

function strip_heatmap(x, y, z; colormap = cgrad(:inferno), bg = :black)
    n_x = length(x)
    x_min = min(x...)
    x_max = max(x...)
    x_range = x_max - x_min

    n_y = length(y)
    y_min = min(y...)
    y_max = max(y...)
    y_range = y_max - y_min

    heatmap(x, y, z,
        background_color_inside = bg,
        colormap = colormap,
        framestyle = :grid,
        grid = true,
        xlims = (x_min, x_max),
        ylims = (y_min, y_max),
        aspect_ratio = (x_range / y_range) * .15,
        dpi = 700
    )
    return plot!()
end

function energy_time_series(event::Event; by = "index", save = false, show_plot = true)
    if by == "index"
        x = eachindex(event.time)
        x_label = "Data Index"
        obs_edges = event.observation_edge_idxs
    elseif by == "time"
        x = event.time
        x_label = "Time (s)"
        obs_edges = event.time[event.observation_edge_idxs[1:end-1]]
    elseif by == "date"
        x = event.time
        x_label = "Time (UTC), $(Date(event.time_datetime[1]))"
        obs_edges = event.time[event.observation_edge_idxs[1:end-1]]
    else
        error("Keyword argument by=\"$(by)\" not recognized.")
    end

    energy = event.energy_bins_mean ./ 1000 # MeV
    e_flux, _ = integrate_flux(event, pitch_angle = true)
    strip_heatmap(x, log10.(energy), log10.(e_flux'), colormap = :ice)

    plot!(
        title = "Pitch Angle Integrated Flux",
        xlabel = x_label,
        ylabel = "Energy (MeV)",
        colorbar_title = "Log10 Energy Flux\n(keV/(cm² MeV s))",
        clims = (5, 11.5),
        yticks = (log10.(energy[1:3:16]), round.(energy[1:3:16], sigdigits = 2)),
        margin = 5mm,
        size = (1.3, .33) .* 500
    )
    # Add date ticks if needed
    if by == "date"
        tick_idxs = Int.(round.(LinRange(1, event.n_datapoints, 6)))
        ticks = Dates.format.(event.time_datetime[tick_idxs], "HH:MM:SS")
        plot!(xticks = (event.time[tick_idxs], ticks))
    end
    vline!(obs_edges,
        linewidth = 1,
        linecolor = :white,
        linestyle = :dash,
        label = ""
    )
    if show_plot == true
        display(plot!())
    end
    if save == true
        png("../results/figures/time_series_$(event.name)")
    end
    return plot!()
end

function J_over_J90_time_series(event::Event; by = "index", save = false, show_plot = true)
    if by == "index"
        x = eachindex(event.time)
        x_label = "Data Index"
        obs_edges = event.observation_edge_idxs
    elseif by == "time"
        x = event.time
        x_label = "Time (s)"
        obs_edges = event.time[event.observation_edge_idxs[1:end-1]]
    elseif by == "date"
        x = event.time
        x_label = "Time (UTC), $(Date(event.time_datetime[1]))"
        obs_edges = event.time[event.observation_edge_idxs[1:end-1]]
    else
        error("Keyword argument by=\"$(by)\" not recognized.")
    end

    # Create diverging colorbar centered on 0
    colorbar_range = (-1.25, .25)
    color_anchors = [colorbar_range[1], 0, colorbar_range[2]]
    color_anchors = (color_anchors .- min(color_anchors...)) ./ (max(color_anchors...) - min(color_anchors...)) # Map to [0,1]
    colorbar = cgrad(:grays, color_anchors)
    c = cgrad(:cherry, rev = true)
    

    energy = event.energy_bins_mean ./ 1000 # MeV
    strip_heatmap(x, log10.(energy), log10.(event.Jprec_over_Jperp'), bg = RGB(.8,.8,.8), colormap = :ice)

    plot!(
        title = "Jprec/Jperp",
        xlabel = x_label,
        ylabel = "Energy (MeV)",
        colorbar_title = "Log10 Jprec/Jperp",
        background = :transparent,
        grid = false,
        clims = (-1.25, .25),
        yticks = (log10.(energy[1:3:16]), round.(energy[1:3:16], sigdigits = 2)),
        margin = 5mm,
        size = (1.3, .33) .* 500,
        dpi = 700
    )
    # Add date ticks if needed
    if by == "date"
        tick_idxs = Int.(round.(LinRange(1, event.n_datapoints, 6)))
        ticks = Dates.format.(event.time_datetime[tick_idxs], "HH:MM:SS")
        plot!(xticks = (event.time[tick_idxs], ticks))
    end
    vline!(obs_edges,
        linewidth = 1,
        linecolor = :white,
        linestyle = :dash,
        label = ""
    )
    if show_plot == true
        display("image/png", plot!())
    end
    if save == true
        png("../results/figures/time_series_prec_over_perp_$(event.name)")
    end
    return plot!()
end

function pad_time_series(event::Event; by = "index", show_plot = true)
    if by == "index"
        time = collect(eachindex(event.time))
        x_label = "Data Index"
        obs_edges = event.observation_edge_idxs
    elseif by == "time"
        time = event.time
        x_label = "Time (s)"
        obs_edges = event.time[event.observation_edge_idxs[1:end-1]]
    elseif by == "date"
        time = event.time
        x_label = "Time (UTC), $(Date(event.time_datetime[1]))"
        obs_edges = event.time[event.observation_edge_idxs[1:end-1]]
    else
        error("Keyword argument by=\"$(by)\" not recognized.")
    end
    
    # Get pitch angle bin edges
    distance_between_bins = diff(event.pitch_angles, dims = 2)
    pitch_angle_bin_edges = event.pitch_angles[:, 1:end-1] .+ (distance_between_bins ./ 2)
    # Add first and last edge
    pitch_angle_bin_edges = cat(event.pitch_angles[:,1] .- (distance_between_bins[:,1] ./ 2), pitch_angle_bin_edges, dims = 2)
    pitch_angle_bin_edges = cat(pitch_angle_bin_edges, event.pitch_angles[:,end] .+ (distance_between_bins[:,end] ./ 2), dims = 2)

    # Get fluxes
    e_flux, n_flux = integrate_flux(event, energy = true)
    flux = log10.(e_flux) # Flux to plot

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
        background_color_inside = RGB(.8,.8,.8),
        background_color_outside = :transparent,
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
    # Add date ticks if needed
    if by == "date"
        tick_idxs = Int.(round.(LinRange(1, event.n_datapoints, 6)))
        ticks = Dates.format.(event.time_datetime[tick_idxs], "HH:MM:SS")
        plot!(xticks = (event.time[tick_idxs], ticks))
    end
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
    vline!(obs_edges,
        linewidth = 1,
        linecolor = :white,
        linestyle = :dash,
        label = ""
    )

    if show_plot == true
        display("image/png", plot!())
    end
    return plot!()
end


function L_MLT_time_series(event::Event; by = "index", save = false, show_plot = true)
    if by == "index"
        x = eachindex(event.time)
        x_label = "Data Index"
        obs_edges = event.observation_edge_idxs
    elseif by == "time"
        x = event.time
        x_label = "Time (s)"
        obs_edges = event.time[event.observation_edge_idxs[1:end-1]]
    elseif by == "date"
        x = event.time
        x_label = "Time (UTC), $(Date(event.time_datetime[1]))"
        obs_edges = event.time[event.observation_edge_idxs[1:end-1]]
    else
        error("Keyword argument by=\"$(by)\" not recognized.")
    end

    n_x = length(x)
    x_min = min(x...)
    x_max = max(x...)
    x_range = x_max - x_min

    y = event.L
    y_range = 15 - 0

    plot(
        xlabel = x_label,
        xlims = (x_min, x_max),

        ylabel = "L (T87)",
        ylims = (0, 15),
        yticks = 0:5:15,
        minorgrid = true,
        yminorticks = 5, 

        aspect_ratio = (x_range / y_range) * .15,
        framestyle = :box,
        leftmargin = -5.2mm
    )
    # Add date ticks if needed
    if by == "date"
        tick_idxs = Int.(round.(LinRange(1, event.n_datapoints, 6)))
        ticks = Dates.format.(event.time_datetime[tick_idxs], "HH:MM:SS")
        plot!(xticks = (event.time[tick_idxs], ticks))
    end


    # Plot each observation period
    for i = 1:event.n_observations
        slice = event.observation_edge_idxs[i]:(event.observation_edge_idxs[i+1]-1)
        plot!(x[slice], y[slice],
            label = "",
            linewidth = 1.5,
            linecolor = :black
        )
    end



    for i = 1:event.n_observations
        slice = event.observation_edge_idxs[i]:(event.observation_edge_idxs[i+1]-1)
        plot!(twinx(), x[slice], event.MLT[slice],
            xlims = (x_min, x_max),
            ylims = (0, 24),
            ylabel = "MLT",
            yguidefontcolor = :grey,
            y_foreground_color_axis = :grey,
            y_foreground_color_text = :grey,
            aspect_ratio = (x_range / 24) * .15,
            label = "",
            linewidth = 1.5,
            linecolor = RGBA(0,0,0,.3)
        )
    end




    # Plot belt L ranges
    plot!(Shape([x_min, x_min, x_max, x_max], [1, 2.5, 2.5, 1]),
        fillcolor = RGBA(0,0,0,.1),
        linecolor = RGBA(0,0,0,0),
        label = ""
    )
    plot!(Shape([x_min, x_min, x_max, x_max], [3, 6, 6, 3]),
        fillcolor = RGBA(0,0,0,.1),
        linecolor = RGBA(0,0,0,0),
        label = ""
    )
    vline!(obs_edges,
        linewidth = 1,
        linecolor = :black,
        linestyle = :dash,
        label = ""
    )
    if show_plot == true
        display(plot!())
    end
    if save == true
        png("../results/figures/time_series_$(event.name)")
    end
    return plot!()
end

function absolute_sunset(event::Event; show_plot = true)
    x = event.avg_pitch_angles
    y = log10.(event.energy_bins_mean)
    e_fluence, n_fluence = integrate_flux(event, time = true)
    e_flux = e_fluence ./ event.duration

    n_x = length(x)
    x_min = min(x...)
    x_max = max(x...)
    x_range = 180

    n_y = length(y)
    y_min = min(y...)
    y_max = max(y...)
    y_range = y_max - y_min

    plot(
        title = "Absolute PAD",
        xlabel = "Pitch Angle (deg)",
        ylabel = "Energy (keV)",
        colorbar_title = "Avg. Flux, \nlog10(kev/(cm² str Mev))",
        xlims = (0, 180),
        ylims = (y_min, y_max),
        clims = (3.5, 8),
        yticks = (y[1:2:end], Int.(round.(10 .^ y[1:2:end], sigdigits = 2))),
        aspect_ratio = (x_range / y_range) * 1,
        background_color_inside = :black,
        background_color_outside = :transparent,
        framestyle = :grid,
        grid = false,
        dpi = 700
    )
    heatmap!(x, y, log10.(e_flux),
        colormap = :ice #:cherry
    )
    vline!([event.avg_loss_cone_angle],
        linewidth = 5,
        linecolor = RGB(.5,.5,.5),
        linestyle = :solid,
        label = ""
    )
    vline!([event.avg_anti_loss_cone_angle],
        linewidth = 5,
        linecolor = RGB(.5,.5,.5),
        linestyle = :dash,
        label = ""
    )
    if show_plot == true
        display(plot!())
    end
    return plot!()
end

function relative_sunset(event::Event; show_plot = true)
    # Normalize fluence to 90º fluence (J/J90º)
    e_fluence, _ = integrate_flux(event, time = true)
    e_flux = e_fluence ./ event.duration
    _, normalizing_idx = findmin(abs.(90 .- event.avg_pitch_angles))
    for e = 1:16 # Normalize each energy bin (row)
        J90 = e_flux[e, normalizing_idx]
        if J90 < 10^5.25
            e_flux[e,:] .= 0
        else
            e_flux[e,:] ./= J90
        end
    end
    e_flux = log10.(e_flux)

    #=
    # Create diverging colorbar centered on 0
    colorbar_range = (-2.5, .5)
    color_anchors = [colorbar_range[1], 0, colorbar_range[2]]
    color_anchors = (color_anchors .- min(color_anchors...)) ./ (max(color_anchors...) - min(color_anchors...)) # Map to [0,1]
    colorbar = cgrad(:vik, color_anchors)
    =#
    colorbar = :ice
    colorbar_range = (-2, 0)

    x = event.avg_pitch_angles
    n_x = length(x)
    x_min = min(x...)
    x_max = max(x...)
    x_range = 180

    y = log10.(event.energy_bins_mean)
    finite_idxs = isfinite.(y)
    y_min = min(y[finite_idxs]...)
    y_max = max(y[finite_idxs]...)
    y_range = y_max - y_min

    plot(
        title = "J/J₉₀ PAD",
        xlabel = "Pitch Angle (deg)",
        ylabel = "Energy (keV)",
        colorbar_title = "log10(J/J₉₀)",
        xlims = (0, 180),
        ylims = (y_min, y_max),
        clims = colorbar_range,
        yticks = (y[1:2:end], Int.(round.(10 .^ y[1:2:end], sigdigits = 2))),
        aspect_ratio = (x_range / y_range) * 1,
        background_color_inside = :black,
        background_color_outside = :transparent,
        framestyle = :grid,
        grid = false,
        dpi = 700
    )
    heatmap!(x, y, e_flux,
        colormap = colorbar
    )
    vline!([event.avg_loss_cone_angle],
        linewidth = 5,
        linecolor = RGB(.5,.5,.5),
        linestyle = :solid,
        label = ""
    )
    vline!([event.avg_anti_loss_cone_angle],
        linewidth = 5,
        linecolor = RGB(.5,.5,.5),
        linestyle = :dash,
        label = ""
    )
    if show_plot == true
        display(plot!())
    end
    return plot!()
end

function MLT_L_track(event::Event; index = false, show_plot = true)
    net_e_flux, _ = integrate_flux(event, energy = true, pitch_angle = true)
    net_e_flux = log10.(net_e_flux)
    net_e_flux[isinf.(net_e_flux)] .= 0

    plot(
        proj = :polar,
        xaxis = false,
        yaxis = false,
        ylims = (0, 8),
        background_color = :transparent,
        grid = false,
        dpi = 700
    )
    # Black background
    plot!(Shape(vcat(Plots.partialcircle(0, 2π, 100, 0), reverse(Plots.partialcircle(0, 2π, 100, 1)))), fillcolor = :black, label = "")
    # L-shells
    for L = 2:2:8
        plot!(LinRange(0,2π,100), zeros(100) .+ L,
            linecolor = RGBA(.5,.5,.5),
            label = ""
        )
    end
    # Radiation belts
    inner_belt = [1, 2.5] ./ 8    # Because plotting circles only works with relative units, we need to get the radiation belt ranges in relative units too
    outer_belt = [3, 6]   ./ 8
    plot!(Shape(vcat(Plots.partialcircle(0, 2π, 100, inner_belt[1]), reverse(Plots.partialcircle(0, 2π, 100, inner_belt[2])))),
        fillcolor = RGBA(1,1,1,.25),
        linecolor = RGBA(0,0,0,0),
        label = ""
    )
    plot!(Shape(vcat(Plots.partialcircle(0, 2π, 100, outer_belt[1]), reverse(Plots.partialcircle(0, 2π, 100, outer_belt[2])))),
        fillcolor = RGBA(1,1,1,.25),
        linecolor = RGBA(0,0,0,0),
        label = ""
    )
    # Earth
    plot!(Shape(vcat(Plots.partialcircle(0, 2π, 100, 0), reverse(Plots.partialcircle(0, 2π, 100, 1/15)))),
        fillcolor = :black,
        linecolor = RGBA(0,0,0,0),
        label = ""
    )
    plot!(Shape(vcat(Plots.partialcircle(3π/2, π/2, 100, 0), reverse(Plots.partialcircle(3π/2, π/2, 100, 1/15)))),
        fillcolor = :white,
        linecolor = RGBA(0,0,0,0),
        label = ""
    )
    # Spacecraft path & flux
    MLT_radians = (event.MLT/24) * 2π
    scatter!(MLT_radians, event.L,
        linewidth = 4,
        zcolor = net_e_flux,
        color = :inferno,
        markerstrokewidth = 0,
        label = "",
        colorbar = false,
        #clims = (4,8)
    )
    # MLT Labels
    annotate!(1.4, -.07, text("0000", 12))
    annotate!(-1.35, -.07, text("1200", 12))
    annotate!(0, 1.1, text("0600", 12))
    annotate!(0, -1.1, text("1800", 12))

    if show_plot == true
        display(plot!())
    end
    return plot!()
end


##############################################################################
#                                                                            #
#                              UNUSED FUNCTIONS                              #
#                                                                            #
##############################################################################

function energy_spectrum(event::Event; show_plot = true)
    energy = log10.(event.energy_bins_mean)

    fluence, _ = integrate_flux(event, time = true, pitch_angle = true)
    #=
    fluence_normalizer = 1 / integrate(energy, fluence, Trapezoidal())
    fluence .*= fluence_normalizer
    =#
    fluence = log10.(fluence)

    x = fluence
    n_x = length(x)
    x_min = min(x...)
    x_max = max(x...)
    x_range = x_max - x_min

    y = energy
    n_y = length(y)
    y_min = min(y...)
    y_max = max(y...)
    y_range = y_max - y_min

    plot(energy, fluence,
        marker = true,
        markercolor = :black,
        linecolor = :black,
        background_color = :transparent,
        legend = false,
        grid = false,
    #=
        title = "Energy Spectrum",
        xlabel = "log10 Electron Fluence (MeV⁻¹)\n",
        ylabel = "Energy (keV)",
        #xlims = (x_min - (x_pixel_size/2), x_max + (x_pixel_size/2)),
        #ylims = (y_min - (y_pixel_size/2), y_max + (y_pixel_size/2)),
        #aspect_ratio = (x_range / y_range),
        #size = (.15, .15 * 1.8) .* 500,
        dpi = 700=#
    )
    if show_plot == true
        display(plot!())
    end
    return plot!()
end

function plot_ring!(radius; color = RGBA(.5,.5,.5,.9))
    # Plots a ring in polar projection.
    circle = LinRange(0, 2π, 200)
    plot!(circle, (circle .* 0) .+ radius,
        linecolor = color,
        label = ""
    )
    return plot!()
end

function my_2d_histogram(x, y, x_bin_edges, y_bin_edges; weights = ones(length(x)))
    @assert length(x) == length(y) == length(weights) "x, y, and weight vectors must be same length"

    x_nbins = length(x_bin_edges) - 1
    y_nbins = length(y_bin_edges) - 1
    
    result = zeros(x_nbins, y_nbins)
    
    for x_idx = 1:x_nbins
        x_slice = x_bin_edges[x_idx] .<= x .< x_bin_edges[x_idx+1]
        for y_idx = 1:y_nbins
            y_slice = y_bin_edges[y_idx] .<= y .< y_bin_edges[y_idx+1]
            slice = x_slice .&& y_slice
            result[x_idx, y_idx] = sum(weights[slice])
        end
    end
    return result
end

function my_1d_histogram(x, bin_edges)
    weights = ones(length(x))
    return my_1d_histogram(x, weights, bin_edges)
end

function my_1d_histogram(x, weights, bin_edges)
    nbins = length(bin_edges)-1
    result = zeros(nbins)
    for i = 1:nbins
        slice = bin_edges[i] .<= x .< bin_edges[i+1]
        result[i] = sum(weights[slice])
    end
    return result
end