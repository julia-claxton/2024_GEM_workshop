cd(@__DIR__)
println("Compiling libraries...")
include("../Julia_ELFIN_Tools/Events.jl")
include("../Julia_ELFIN_Tools/Visualization.jl")
include("../Julia_ELFIN_Tools/Event_Detection.jl")
using Statistics
using LinearAlgebra
using Plots

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
    ylims = (.01, 2),
    yscale = :log10,
    yminorticks = true,

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

    ylims = (.01, 2),
    yscale = :log10,
    yminorticks = true,
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

    ylims = (.01, 2),
    yscale = :log10,
    yminorticks = true,
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

png("../../results/figures/supplementary/energy_profile.png")
display(plot!())