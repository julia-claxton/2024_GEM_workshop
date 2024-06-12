cd(@__DIR__)
println("Compiling libraries...")
include("../Julia_ELFIN_Tools/Events.jl") # Provides type Event for ELFIN data
include("../Julia_ELFIN_Tools/Visualization.jl")
include("../Julia_ELFIN_Tools/Event_Detection.jl")
using Statistics
using LinearAlgebra
using Dates
using DelimitedFiles
using TickTock
using Glob


mlt = []
dst = []

dirpath = "../../data/processed_position_data"
files = glob("*.npz", dirpath)
dates = Array{Date}(undef, length(files))
sats = copy(files)
N = length(files)

# Get MLT and dst for all datapoints
for i = 1:N
    print("\r\t$(i)/$(N)")
    dates[i] = Date(files[i][end-15:end-8], dateformat"yyyymmdd")
    sats[i] = string(files[i][end-4])

    event = create_event(dates[i], sats[i])
    if event == nothing; continue; end

    append!(mlt, event.MLT)
    append!(dst, event.dst)

    if i % 30 == 0; GC.gc(); end
end

# Bin into 2D histogram
MLT_nbins = 24 * 6
dst_nbins   = 14 * 2

MLT_bin_edges = LinRange(0, 24, MLT_nbins + 1)
dst_bin_edges = LinRange(-100, 40, dst_nbins + 1)

result = zeros(dst_nbins, MLT_nbins)

for mlt_idx = 1:MLT_nbins
    print("\r$(mlt_idx)/$(MLT_nbins)")
    mlt_slice = MLT_bin_edges[mlt_idx] .<= mlt .< MLT_bin_edges[mlt_idx+1]
    for dst_idx = 1:dst_nbins
        dst_slice = dst_bin_edges[dst_idx] .<= dst .< dst_bin_edges[dst_idx+1]
        slice = mlt_slice .&& dst_slice

        result[dst_idx, mlt_idx] += sum(slice)
    end
    # Then normalize each MLT slice to be a pdf (i.e. integrates to 1)
    result[:, mlt_idx] ./= sum(result[:, mlt_idx]) # We can use sum instead of integration here because Δdst is uniform and would thus get divided out using rectangular integration, and it's computationally faster
end


# Create contour plot
contour(log10.(result),
    levels = 5,
    fill = true,

    title = "Geomagnetic Activity Coverage by MLT",
    xlabel = "MLT, hour (IGRF)",
    ylabel = "dst",
    colorbar_title = "Log10 Density",

    xlims = (0, MLT_nbins),
    xticks = (LinRange(0, MLT_nbins, length(0:3:24)), 0:3:24),
    ylims = (0, dst_nbins),
    yticks = (LinRange(0, dst_nbins, length(-100:10:40)), -100:10:40),
    clims = (-4.25, -.5),

    colormap = :bone_1,
    background = :transparent,
    framestyle = :box,
    grid = true,
    xminorgrid = true,
    xminorticks = 3,
    dpi = 300
)
display("image/png", plot!())


println("\nDone")



