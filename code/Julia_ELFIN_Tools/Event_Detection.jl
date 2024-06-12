include("./Events.jl")
include("./Visualization.jl")
using Statistics
using LinearAlgebra
using Dates
using Plots
using Plots.PlotMeasures
using StatsBase
using NumericalIntegration
using NPZ
using HDF5

#=
Event_Detection.jl

This file provides functions for the detection of EMIC-induced precipitation in ELFIN data, including
functions that calculate precipitation scores for high and low energy channels, and an EIP detection
function.

Written by Julia Claxton (julia.claxton@colorado.edu)
Released under MIT License (see <project top level>/LICENSE.txt for full license)
=#

function check_for_bad_events(fullday_event::Event, slice)
    # Low energy channels hot
    e_flux, _ = integrate_flux(fullday_event, energy = true)
    slice_avg = mean(e_flux[slice,:])
    slice_std = std(e_flux[slice,:])

    if log10(slice_std) <= log10(slice_avg) - .5 # These are events with hot low channels and a very hot and uniform PAD.
        return true
    end
end

function high_energy_precipitation_score(event::Event; time_range_idx = 1:event.n_datapoints)
    Jprec_over_Jperp = event.Jprec_over_Jperp[time_range_idx, :]

    # Kill pixels with low amounts of precipitation
    Jprec_over_Jperp[Jprec_over_Jperp .< 10^(-.75)] .= NaN

    # Find pixels with few neighbors and kill them
    old_J_Jprec_over_Jperp = copy(Jprec_over_Jperp) # To avoid feedback loop of pixel killing
    for t = 1:length(time_range_idx)
        t_to_check = clamp.(t-1:t+1, 1, length(time_range_idx))
        t_to_check = unique(t_to_check) # Don't double count pixels on the edges

        for e = 1:16
            e_to_check = clamp.(e-1:e+1, 1, 16)
            e_to_check = unique(e_to_check) # Don't double count pixels on the edges

            n_neighbors = sum(isfinite.(old_J_Jprec_over_Jperp[t_to_check, e_to_check])) - 1 # Minus one because this includes current pixel. It can go below zero if current pixel is non-finite but that's fine.
            if n_neighbors < 3
                Jprec_over_Jperp[t,e] = NaN
            end
        end
    end
    nonfinite_idxs = findall(.!isfinite.(Jprec_over_Jperp))
    Jprec_over_Jperp[nonfinite_idxs] .= 0

    # Crop out low energy channels
    high_energy_idxs = findall(event.energy_bins_min .>= 500)
    Jprec_over_Jperp_cropped = Jprec_over_Jperp[:, high_energy_idxs]

    # Sum up for score
    score = dropdims(sum(Jprec_over_Jperp_cropped, dims = 2), dims = 2)

    # Return
    return score
end

function low_energy_precipitation_score(event::Event; time_range_idx = 1:event.n_datapoints)
    Jprec_over_Jperp = event.Jprec_over_Jperp[time_range_idx ,:]

    # Zero out low precipitation
    Jprec_over_Jperp[Jprec_over_Jperp .< 10^(-.75)] .= 0
    
    # Remove NaNs and Infs
    nonfinite_idxs = findall(.!isfinite.(Jprec_over_Jperp))
    Jprec_over_Jperp[nonfinite_idxs] .= 0

    # Crop out high energy channels
    low_energy_idxs = findall(event.energy_bins_min .<= 500)
    Jprec_over_Jperp = Jprec_over_Jperp[:, low_energy_idxs]

    # Sum up for score
    score = dropdims(sum(Jprec_over_Jperp, dims = 2), dims = 2)
end

function quiet_detector(event::Nothing; time_range_idx = nothing, show_plot = nothing)
    return nothing
end

function quiet_detector(event::Event; time_range_idx = 1:event.n_datapoints, show_plot = true)
    if length(time_range_idx) < 2; return [], []; end

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

    low_e_score = smoothed_low_e

    mask = (high_e_score .<= .25) .&& (low_e_score .<= .25)
    mask[end] = 0 # If the last datapoint is high, call it the end of an event. If it's low, no impact.
    event_edges = diff(mask)
    if event_edges[1] == -1
        event_edges[1] = 0 # If first datapoint is high and second low, discard detection
    else
        event_edges[1] = mask[1] # If the first datapoint is high, call it the start of an event
    end

    # Get info on detected events
    start_idxs = findall(event_edges .== 1)
    stop_idxs = findall(event_edges .== -1)

    # Plotting
    if show_plot == true
        # Heatmap
        yticks = Int.(round.(event.energy_bins_mean))
        yticks = yticks[1:3:16]
        p1 = heatmap(time_range_idx, 1:16, log10.(event.Jprec_over_Jperp[time_range_idx,:]'),
            background_color_inside = RGB(.75,.75,.75),
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
        )
        plot!(time_range_idx, high_e_score, label = "High-Energy Precip. Score")
        plot!(time_range_idx, low_e_score, label = "Low-Energy Precip. Score")
        p2 = plot!()

        # Detection area
        plot(time_range_idx, 1 .* mask,
            fill = true,
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
            suptitle = "$(event.name)",
            dpi = 350
        )
        display(plot!())
    end
    return time_range_idx[start_idxs], time_range_idx[stop_idxs]
end

function emic_detector(event::Nothing; time_range_idx = nothing, show_plot = nothing)
    return nothing
end
  
function emic_detector(event::Event; time_range_idx = 1:event.n_datapoints, show_plot = true)
    if length(time_range_idx) < 2; return [], []; end

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
    if show_plot == true
        # Heatmap
        yticks = Int.(round.(event.energy_bins_mean))
        yticks = yticks[1:3:16]
        p1 = heatmap(time_range_idx, 1:16, log10.(event.Jprec_over_Jperp[time_range_idx,:]'),
            background_color_inside = RGB(.75,.75,.75),
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
        )
        plot!(time_range_idx, high_e_score, label = "High-Energy Precip. Score")
        plot!(time_range_idx, smoothed_low_e, label = "Low-Energy Precip. Score")
        p2 = plot!()

        # Detection area
        plot(time_range_idx, 1 .* emic_mask,
            fill = true,
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
            suptitle = "$(event.name)",
            dpi = 350
        )
        display("image/png",plot!())
    end
    return time_range_idx[emic_start_idxs], time_range_idx[emic_stop_idxs]
end