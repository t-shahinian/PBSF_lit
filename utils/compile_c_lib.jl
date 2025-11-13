#!/usr/bin/env julia

"""
compile_c_lib.jl

This script compiles the C++ library for nearest neighbor calculations.
It creates the necessary directory structure and compiles the C++ code
into a shared library usable by Julia.

Usage:
- Serial build: julia utils/compile_c_lib.jl
- OpenMP build: julia utils/compile_c_lib.jl --openmp
"""

# Parse command line arguments
use_openmp = "--openmp" in ARGS

# Get the absolute path of the script's directory
script_dir = @__DIR__
lib_dir = joinpath(script_dir, "julia-R-nn-ccall2")

# Ensure the directory exists
mkpath(lib_dir)

# Full paths to source and output files
cpp_file = joinpath(lib_dir, "nn.cpp")
obj_file = joinpath(lib_dir, "nn.o")

# Check if source file exists
if !isfile(cpp_file)
    error("C++ source file not found at $cpp_file")
end

# Compile the C++ code
try
    if Sys.isapple()
        # macOS: use Homebrew LLVM clang++
        clangxx = "/opt/homebrew/opt/llvm/bin/clang++"
        if !isfile(clangxx)
            error("LLVM clang++ not found at $clangxx. Install with: brew install llvm")
        end

        if use_openmp
            run(`$clangxx -c -fpic -fopenmp $cpp_file -o $obj_file`)
            lib_file = joinpath(lib_dir, "libnn.dylib")
            run(`$clangxx -shared -fopenmp -o $lib_file $obj_file`)
            println("Compiled libnn.dylib (macOS, OpenMP, LLVM clang++)")
        else
            run(`$clangxx -c -fpic $cpp_file -o $obj_file`)
            lib_file = joinpath(lib_dir, "libnn.dylib")
            run(`$clangxx -shared -o $lib_file $obj_file`)
            println("Compiled libnn.dylib (macOS, serial, LLVM clang++)")
        end

    elseif Sys.iswindows()
        # Windows
        if use_openmp
            run(`clang++ -c -fpic -fopenmp $cpp_file -o $obj_file`)
            lib_file = joinpath(lib_dir, "libnn.dll")
            run(`clang++ -shared -fopenmp -o $lib_file $obj_file`)
            println("Compiled libnn.dll (Windows, OpenMP)")
        else
            run(`clang++ -c -fpic $cpp_file -o $obj_file`)
            lib_file = joinpath(lib_dir, "libnn.dll")
            run(`clang++ -shared -o $lib_file $obj_file`)
            println("Compiled libnn.dll (Windows, serial)")
        end

    else
        # Linux and other Unix
        if use_openmp
            run(`clang++ -c -fpic -fopenmp $cpp_file -o $obj_file`)
            lib_file = joinpath(lib_dir, "libnn.so")
            run(`clang++ -shared -fopenmp -o $lib_file $obj_file`)
            println("Compiled libnn.so (Linux, OpenMP)")
        else
            run(`clang++ -c -fpic $cpp_file -o $obj_file`)
            lib_file = joinpath(lib_dir, "libnn.so")
            run(`clang++ -shared -o $lib_file $obj_file`)
            println("Compiled libnn.so (Linux, serial)")
        end
    end

    println("Compilation complete!")

catch e
    println("Compilation failed:")
    println(e)
    rethrow(e)
end

