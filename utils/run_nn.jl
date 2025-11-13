#!/usr/bin/env julia

"""
run_nn.jl

This script demonstrates how to use the compiled C library for nearest neighbor calculations.
It reads coordinates from a CSV file, finds the nearest neighbors, and prints the results.
"""

using CSV
using DataFrames
using DelimitedFiles
using PrettyTables

# Function to find the path to the compiled library
function get_lib_path()
    if Sys.isapple()
        return joinpath(@__DIR__, "julia-R-nn-ccall2", "libnn.dylib")
    elseif Sys.iswindows()
        return joinpath(@__DIR__, "julia-R-nn-ccall2", "libnn.dll")
    else
        return joinpath(@__DIR__, "julia-R-nn-ccall2", "libnn.so")
    end
end

# Check if the library exists
lib_path = get_lib_path()
if !isfile(lib_path)
    error("Library not found at $lib_path. Please run compile_c_lib.jl first.")
end

# Read coordinates from CSV file
coords_file = joinpath(@__DIR__, "julia-R-nn-ccall2", "coords.csv")
if !isfile(coords_file)
    error("Coordinates file not found at $coords_file")
end

# Read raw data, skipping the header
raw_data = readdlm(coords_file, ',', skipstart=1)
coords = raw_data[:, 1:2]

# Flatten the coordinates for C function
coords_flat = vcat(coords[:, 1], coords[:, 2])

# Parameters
n = size(coords, 1)  # Number of points
m = 5               # Number of nearest neighbors
nThreads = 1        # Force serial execution

# Calculate the size of the index array
nIndx = trunc(Int, (1+m)/2*m+(n-m-1)*m)

# Initialize arrays with safe initialization
nnIndx = zeros(Int32, nIndx)
nnDist = zeros(Float64, nIndx)
nnIndxLU = zeros(Int32, 2*n)

# Convert parameters to pointers for C call
nPtr = Ref(Int32(n))
mPtr = Ref(Int32(m))
nThreadsPtr = Ref(Int32(nThreads))

# Ensure coordinates are in correct format
coords_ptr = pointer(coords_flat)

# Call the C function
ccall(
    (:mkNNIndxCB, lib_path),
    Cvoid,
    (Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, Ref{Int32}),
    nPtr, mPtr, coords_ptr, nnIndx, nnDist, nnIndxLU, nThreadsPtr
)

# Print results
println("Nearest neighbor calculations complete!")
println("Number of points: $n")
println("Number of nearest neighbors: $m")
println("Number of threads used: $nThreads")

# Create a table with the results
results = DataFrame(
    Point = Int[],
    Neighbor = Int[],
    Distance = Float64[]
)

# Safe indexing with bounds checking
for i in 1:n
    start_idx = nnIndxLU[i]
    num_neighbors = nnIndxLU[n+i]
    
    # Ensure start_idx is valid
    if start_idx < 1 || start_idx > length(nnIndx)
        println("Warning: Invalid start index for point $i")
        continue
    end
    
    for j in 1:min(num_neighbors, m)
        # Ensure we don't go out of bounds
        idx = start_idx + j - 1
        if idx > 0 && idx <= length(nnIndx) && 
           nnIndx[idx] > 0 && 
           nnIndx[idx] <= n
            
            neighbor_idx = nnIndx[idx]
            distance = nnDist[idx]
            
            push!(results, (i, neighbor_idx, distance))
        end
    end
end

# Print the first 10 results
println("\nFirst 10 nearest neighbors:")
if !isempty(results)
    pretty_table(results[1:min(10, nrow(results)), :])
end

println("\nResults saved to 'nn_results.csv'")
CSV.write("nn_results.csv", results) 