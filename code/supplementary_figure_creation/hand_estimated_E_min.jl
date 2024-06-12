cd(@__DIR__)
println("Compiling libraries...")
include("../Julia_ELFIN_Tools/Events.jl")
include("../Julia_ELFIN_Tools/Visualization.jl")
include("../Julia_ELFIN_Tools/Event_Detection.jl")
using Statistics
using LinearAlgebra
using Plots

precip_ratio = 10 .^(-1.6:.2:-.4)
E_min = [140, 170, 170, 180, 310, 530, 1020] 

plot(precip_ratio, E_min,
    marker = true,
    markercolor = RGB(0xb3a2de),
    markerstrokecolor = :black,
    linecolor = RGB(0xb3a2de),
    linewidth = 2.5,
    label = false,

    title = "Minimum Scattering Energy \nby Scattering Strength\n",
    titlefontcolor = :white,

    xlabel = "J/J90",
    xlims = (.01, 1),
    xscale = :log10,
    xminorticks = true,
    xflip = true,
    xguidefont = :white,
    xtickfontcolor = :white,

    ylabel = "Eₘᵢₙ, keV",
    ylims = (100, 1100),
    yscale = :log10,
    yminorticks = true,
    yguidefont = :white,
    ytickfontcolor = :white,

    grid = true,
    minorgrid = true,

    background_color_inside = :white,
    background_color_outside = BG_COLOR,

    size = (1, 1.1) .* 400,
    dpi = 350
)
plot!(twinx(),
    xlims = (.01, 1),
    xscale = :log10,
    xminorticks = true,
    xflip = true,
    xguidefont = :white,
    xtickfontcolor = :transparent,

    ylims = (100, 1100),
    yscale = :log10,
    yminorticks = true,
    yguidefont = :white,
    ytickfontcolor = :transparent,

    grid = true,
    minorgrid = true
)
plot!(twiny(),
    xlims = (.01, 1),
    xscale = :log10,
    xminorticks = true,
    xflip = true,
    xguidefont = :white,
    xtickfontcolor = :transparent,

    ylims = (100, 1100),
    yscale = :log10,
    yminorticks = true,
    yguidefont = :white,
    ytickfontcolor = :transparent,

    grid = true,
    minorgrid = true
)
png("../../results/figures/supplementary/E_min.png")
display(plot!())
