#!/usr/bin/env julia

"""
test_installation.jl

This script tests the installation of the ST_SpatialFactor project.
It verifies that:
1. All required packages are installed
2. The C library can be compiled
3. The nearest neighbor calculations can be run
"""

using BenchmarkTools
using CSV
using DataFrames
using Distributions
using GeoStats
using IJulia
using IterativeSolvers
using JLD2
using Plots
using PrettyTables
using ProgressMeter
using RCall
using StatsPlots

println("All packages loaded successfully!")

# Test RCall
R"""
print("R is working!")
"""

# Test plotting
p = plot(randn(100), title="Test Plot")
display(p)

# Test DataFrame
df = DataFrame(A = 1:4, B = ["M", "F", "F", "M"])
pretty_table(df)

println("All functionality tests passed!")

# Function to test package installation
function test_packages()
    println("Testing package installation...")
    
    # List of required packages
    packages = [
        "CSV",
        "DataFrames",
        "Plots",
        "StatsPlots",
        "GeoStats",
        "JLD2",
        "PrettyTables",
        "Distributions",
        "Random",
        "Distances",
        "LinearAlgebra",
        "SparseArrays",
        "IterativeSolvers",
        "ProgressMeter",
        "BenchmarkTools",
        "IJulia",
        "RandomizedLinAlg",
        "RData",
        "RCall", 
        "MCMCChains",
        "MCMCDiagnosticTools",
        "StatsBase",
        "DelimitedFiles"
    ]
    
    # Check if all packages are installed
    missing_packages = String[]
    for pkg in packages
        try
            @eval using $(Symbol(pkg))
            println("✓ $pkg is installed")
        catch
            push!(missing_packages, pkg)
            println("✗ $pkg is not installed")
        end
    end
    
    if !isempty(missing_packages)
        println("\nThe following packages are missing:")
        for pkg in missing_packages
            println("  - $pkg")
        end
        println("\nPlease run 'julia install_packages.jl' to install them.")
        return false
    end
    
    println("All packages are installed!")
    return true
end

# Function to test C library compilation
function test_c_library()
    println("\nTesting C library compilation...")
    
    # Check if the C source files exist
    cpp_file = joinpath(@__DIR__, "utils", "julia-R-nn-ccall2", "nn.cpp")
    h_file = joinpath(@__DIR__, "utils", "julia-R-nn-ccall2", "nn.h")
    
    if !isfile(cpp_file) || !isfile(h_file)
        println("✗ C source files not found")
        return false
    end
    
    # Try to compile the C library
    try
        include(joinpath(@__DIR__, "utils", "compile_c_lib.jl"))
        println("✓ C library compiled successfully")
        return true
    catch e
        println("✗ Failed to compile C library: $e")
        return false
    end
end

# Function to test nearest neighbor calculations
function test_nn_calculations()
    println("\nTesting nearest neighbor calculations...")
    
    # Check if the coordinates file exists
    coords_file = joinpath(@__DIR__, "utils", "julia-R-nn-ccall2", "coords.csv")
    if !isfile(coords_file)
        println("✗ Coordinates file not found")
        return false
    end
    
    # Try to run the nearest neighbor calculations
    try
        include(joinpath(@__DIR__, "utils", "run_nn.jl"))
        println("✓ Nearest neighbor calculations completed successfully")
        return true
    catch e
        println("✗ Failed to run nearest neighbor calculations: $e")
        return false
    end
end

# Main function
function main()
    println("Testing ST_SpatialFactor installation...\n")
    
    # Test package installation
    packages_ok = test_packages()
    
    # Test C library compilation
    c_library_ok = test_c_library()
    
    # Test nearest neighbor calculations
    nn_calculations_ok = test_nn_calculations()
    
    # Print summary
    println("\nInstallation test summary:")
    println("✓ Package installation: $(packages_ok ? "OK" : "FAILED")")
    println("✓ C library compilation: $(c_library_ok ? "OK" : "FAILED")")
    println("✓ Nearest neighbor calculations: $(nn_calculations_ok ? "OK" : "FAILED")")
    
    if packages_ok && c_library_ok && nn_calculations_ok
        println("\nAll tests passed! The installation is complete.")
    else
        println("\nSome tests failed. Please fix the issues and run the tests again.")
    end
end

# Run the main function
main() 