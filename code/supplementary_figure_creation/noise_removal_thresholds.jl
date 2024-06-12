cd(@__DIR__)
include("../Julia_ELFIN_Tools/Events.jl")
include("../Julia_ELFIN_Tools/Visualization.jl")
include("../Julia_ELFIN_Tools/Event_Detection.jl")
using Statistics
using LinearAlgebra
using Dates
using Glob
using Plots

#########################################
# Functions                             #
#########################################
function get_frequency(bin_edges, data)
    nbins = length(bin_edges)-1
    result = zeros(nbins)
    for i = 1:nbins
        result[i] = sum(bin_edges[i] .<= data .< bin_edges[i+1])
    end
    return result
end

function get_data(n_days)
    dirpath = "../../data/processed_position_data"
    files = glob("*.npz", dirpath)
    dates = copy(files)
    sats = copy(files)
    for i = 1:length(files)
        dates[i] = files[i][end-15:end-8]
        sats[i] = string(files[i][end-4])
    end
    dates = Date.(dates, dateformat"yyyymmdd")

    e_data = Array{Float64}(undef, 0, 16)
    n_data = Array{Float64}(undef, 0, 16)

    # We pick random days because holding all 1500 days of data in memory is 10+ GB
    days_to_get_data_from = sample(1:length(dates), n_days, replace = false)

    for i = 1:n_days
        print("\r\t$(i)/$(n_days)")
        datafile_index = days_to_get_data_from[i]
        event = create_event(dates[datafile_index], sats[datafile_index])
        if event == nothing; continue; end
    
        e_flux, n_flux = integrate_flux(event, pitch_angle = true)
        e_data = cat(e_data, e_flux, dims = 1)
        n_data = cat(n_data, n_flux, dims = 1)
    
        if i % 30 == 0; GC.gc(); end
    end
    return e_data, n_data
end

function calculate_and_plot_thresholds(data, flux_type)
    if flux_type == "n"
        max_threshold = 5
    elseif flux_type == "e"
        max_threshold = 7.75
    else
        error("Argument flux_type must == \"e\" or \"n\"")
    end

    nbins = 75
    bin_edges = LinRange(log10(min(data[data .> 0]...)), log10(max(data...)), nbins + 1)
    bin_means = [mean([bin_edges[i], bin_edges[i+1]]) for i = 1:nbins]

    distro = zeros(nbins, 16)
    peak_idxs = zeros(16)

    for i = 1:16
        channel = data[i,:]
        channel = channel[channel .> 0]
        distro[:,i] = get_frequency(bin_edges, log10.(channel))
        _, peak_idxs[i] = findmax(distro[bin_means .< max_threshold, i])
    end
    peak_idxs = Int.(peak_idxs)

    a = create_event(Date("2020-01-08"), "a")
    e = log10.(a.energy_bins_mean)
    x_range = max(e...) - min(e...)
    y_range = max(bin_means...) - min(bin_means...)

    heatmap(e, bin_means, log10.(distro),
        title = "$(uppercase(flux_type))-Flux Frequency",
        xlabel = "Energy, keV",
        ylabel = "Flux, log10(keV/[MeV cm^2 s])",
        colorbar_title = "log10(Frequency)",
        xlims = (min(e...), max(e...)),
        ylims = (min(bin_means...), max(bin_means...)),
        xticks = (e[1:2:end], Int.(round.(10 .^e[1:2:end]))),
        colormap = :bone_1,
        background_color_inside = :black,
        background_color_outside = :transparent,
        aspect_ratio = (x_range / y_range) * 1,
        dpi = 300
    )

    threshold = round.(bin_means[peak_idxs] .+ .7, sigdigits = 3)

    plot!(e, threshold,
        linestyle = :dash,
        linewidth = 2,
        linecolor = :black,
        label = ""
    )

    println("$(flux_type) thresholds:")
    println(threshold)

    png("../../results/figures/supplementary/$(flux_type)_thresholds.png")
    display("image/png", plot!())
end

#########################################
# Script                                #
#########################################

println("Gathering flux data...")
e_data, n_data = get_data(300)
println()

println("Calculating thresholds...\n")
calculate_and_plot_thresholds(collect(n_data'), "n")
calculate_and_plot_thresholds(collect(e_data'), "e")

println("\nDone")