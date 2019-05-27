using BenchmarkTools, TaylorModels

@static if VERSION > v"0.7.0"
    using LinearAlgebra, SparseArrays
end

SUITE = BenchmarkGroup()  # parent BenchmarkGroup to contain our suite

include("daisy.jl")
