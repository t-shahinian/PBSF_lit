#!/usr/bin/env julia

"""
build.jl

This script configures the R environment for RCall package.
It sets the correct path to the R installation and rebuilds RCall if needed.
"""

# Set R_HOME environment variable to point to the correct R installation
R_HOME = "/opt/homebrew/Cellar/r/4.4.3_1/lib/R"
ENV["R_HOME"] = R_HOME

println("Setting R_HOME to: $R_HOME")

# Check if RCall is installed and rebuild it with the new R_HOME
try
    using Pkg
    Pkg.build("RCall")
    println("RCall has been rebuilt with the new R_HOME setting.")
    
    # Verify RCall works by loading it
    println("Verifying RCall installation...")
    using RCall
    println("RCall loaded successfully!")
catch e
    println("Error configuring RCall: $e")
    println("Please manually rebuild RCall using: ")
    println("   julia -e 'using Pkg; ENV[\"R_HOME\"] = \"$R_HOME\"; Pkg.build(\"RCall\")'")
end

# Ensure we can find the correct libR.dylib
libR_path = joinpath(R_HOME, "lib", "libR.dylib")
if isfile(libR_path)
    println("Found libR.dylib at: $libR_path")
else
    println("Warning: libR.dylib not found at expected location: $libR_path")
    println("You may need to adjust the R_HOME path.")
end
