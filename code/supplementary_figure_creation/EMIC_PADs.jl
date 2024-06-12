cd(@__DIR__)
println("Compiling libraries...")
include("../Julia_ELFIN_Tools/Events.jl")
include("../Julia_ELFIN_Tools/Visualization.jl")
include("../Julia_ELFIN_Tools/Event_Detection.jl")
using Statistics
using LinearAlgebra

#########################################
# Functions                             #
#########################################
function create_pad_base_plot()
    plot(
        xlabel = "NH Pitch Angle, deg",
        xlims = (0, 180),
        xticks = 0:30:180,
        xminorticks = 3,

        ylabel = "J/J90",
        ylims = (.01, 2),
        yscale = :log10,
        
        aspect_ratio = (180/(2-.01)),
        grid = true,
        minorgrid = true,
        topmargin = 5mm,
        bottommargin = 5mm,

        size = (1,1) .* 500,
        dpi = 350
    )
    plot!(twinx(),
        xlims = (0, 180),
        xticks = 0:30:180,
        xminorticks = 3,

        ylims = (.01, 2),
        yscale = :log10,
        grid = true,
        minorgrid = true,
        aspect_ratio = (180/(2-.01))
    )
    plot!(twiny(),
        xlims = (0, 180),
        xticks = 0:30:180,
        xminorticks = 3,

        ylims = (.01, 2),
        yscale = :log10,
        grid = true,
        minorgrid = true,
        aspect_ratio = (180/(2-.01))
    )
    return plot!()
end

function get_emic_pad_subplots()
    plots = [create_pad_base_plot() for i = 1:16]

    for i = 1:n_emic
        print("\r$(i)/$(n_emic)")
        event = create_event(start_times[i], stop_times[i], sat[i])
        e_fluence, n_fluence = integrate_flux(event, time = true)

        for e = 1:16
            # Normalize flux to 90º flux (J/J90º)
            e_flux = e_fluence[e, :] ./ event.duration
            _, normalizing_idx = findmin(abs.(90 .- event.avg_pitch_angles))
            J90 = e_flux[normalizing_idx]
            normalized_flux = e_flux ./ J90

            # Get pitch angles and LC/ALC angles
            α = event.avg_pitch_angles
            α_lc = event.avg_loss_cone_angle
            α_alc = event.avg_anti_loss_cone_angle
            if event.avg_loss_cone_angle > 90
                reverse!(normalized_flux)
                α_lc = 180 - α_lc
                α_alc = 180 - α_alc
            end

            # Plot PAD
            plot!(plots[e], α, normalized_flux, # yaxis is in log10 scale
                title = "$(Int(round(event.energy_bins_mean[e]))) keV",
                linecolor = RGBA(0,0,0,.1),
                label = false
            )
            vline!([α_lc, α_alc],
                linecolor = RGBA(0,0,0,.01),
                label = false
            )
        end
    end
    println()
    return plots
end

#########################################
# Script                                #
#########################################
# Get data from manually-pruned EMIC detection algorithm
data = readdlm("../../results/EMIC_detections.csv", ',')
start_times = DateTime.(data[2:end, 1])
stop_times = DateTime.(data[2:end, 2])
sat = data[2:end, 3]
MLT = data[2:end, 4]
L = data[2:end, 5]
n_emic = length(start_times)

#########################################
# Plots                                 #
#########################################
plots = get_emic_pad_subplots()

plots_we_want = [1, 6, 10, 13]
plot(plots[plots_we_want]...,
    suptitle = "\nEIP Pitch Angle Distributions",
    layout = (4,1),
    size = (1,3.5) .* 450,
    dpi = 350
)
png("../../results/figures/supplementary/EMIC_PADs.png")
display(plot!())


# Animation
animation = @animate for e = 1:16
    plot(plots[e])
end
display(gif(animation, "../../results/figures/supplementary/EMIC_PADs.gif", fps = 5, show_msg = false))