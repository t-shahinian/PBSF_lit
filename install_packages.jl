#!/usr/bin/env julia

"""
install_packages.jl

This script installs all the required Julia packages for the ST_SpatialFactor project.
"""

using Pkg

# Note: For Python-R interoperability, install rpy2 in your Python environment
# conda install -c conda-forge rpy2 -y

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

# Add packages to the project
println("Installing required packages...")
for pkg in packages
    println("Installing $pkg...")
    Pkg.add(pkg)
end

# Special configuration for RCall
try
    println("\nSetting up RCall...")
    
    # Try to find R installation
    r_home = nothing
    try
        r_home = chomp(read(`R RHOME`, String))
        println("Found R installation at: $r_home")
    catch
        println("Could not automatically detect R installation.")
        println("You'll need to manually configure RCall with the deps/build.jl script.")
    end
    
    # If R was found, build RCall with the correct R_HOME
    if r_home !== nothing
        ENV["R_HOME"] = r_home
        println("Building RCall with R_HOME=$r_home")
        Pkg.build("RCall")
        println("RCall setup complete!")
    end
catch e
    println("Error setting up RCall: $e")
    println("You'll need to manually configure RCall with the deps/build.jl script.")
end

println("\nAll packages installed successfully!") 