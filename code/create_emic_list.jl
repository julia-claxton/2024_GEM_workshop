println("\n========= analyze.jl START =========\n")
println("Compiling libraries... ")
cd(@__DIR__)
include("./Julia_ELFIN_Tools/Events.jl") 
include("./Julia_ELFIN_Tools/Visualization.jl") 
include("./Julia_ELFIN_Tools/Event_Detection.jl")
using Statistics
using LinearAlgebra
using Dates
using Plots
using Glob
using DelimitedFiles

#=
create_emic_list.jl

This script is used to generate the list of EIP events for the GEM poster. The EIP detection
algorithm is run on the ELFIN lifetime and each detection is presented to the user for hand
verification. After the verification is complete, the start times, stop times, L, and MLT of
each EIP event is saved to <project top level>/results/EMIC_detections.csv.

Written by Julia Claxton (julia.claxton@colorado.edu)
Released under MIT License (see <project top level>/LICENSE.txt for full license)
=#


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ANALYSIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

print("If using VSCode, type a space and then enter twice because readline() is broken in VSCode. Otherwise, just press enter.")
readline()

# Find all EMIC-driven events
emic_MLT = []
emic_L = []
emic_start_times = []
emic_stop_times = []
emic_satellites = []

data_path = "../data/processed_scientific_data"
files = glob("*.npz", data_path)
N = length(files)

# For each day of data
for i = 1:N
    print("\r\t$(i)/$(N)")
    date = Date(files[i][end-15:end-8], dateformat"yyyymmdd")
    sat = string(files[i][end-4])
    event = create_event(date, sat)
    if event == nothing; continue; end

    # For each observation period in the day
    for obs_number = 1:event.n_observations
        start_idx = event.observation_edge_idxs[obs_number]
        stop_idx = event.observation_edge_idxs[obs_number+1]-1
        slice = start_idx:stop_idx

        emic_start_idxs, emic_stop_idxs = emic_detector(event, time_range_idx = slice, show_plot = false)
        n_emic = length(emic_start_idxs)
        
        if n_emic > 0
            emic_detector(event, time_range_idx = slice, show_plot = true)
            for candidate_number = 1:n_emic
                print("\n\t$(candidate_number): ")
                response = readline()

                if response != ""
                    println("Event rejected")
                    continue
                end
                println("Event accepted")
                append!(emic_start_times, event.time_datetime[emic_start_idxs])
                append!(emic_stop_times, event.time_datetime[emic_stop_idxs])
                append!(emic_MLT, event.MLT[emic_start_idxs])
                append!(emic_L, event.L[emic_start_idxs])
                append!(emic_satellites, repeat(event.satellite, n_emic))
            end
        end
    end
    if i % 30 == 0; GC.gc(); end
end
println()

header = ["start_time" "stop_time" "satellite" "MLT" "L"]
data = hcat(emic_start_times, emic_stop_times, emic_satellites, emic_MLT, emic_L)
to_write = vcat(header, data)
writedlm("../results/EMIC_detections.csv", to_write, ',')

println("\nAnalysis Complete")
println("\n========== analyze.jl END ==========")