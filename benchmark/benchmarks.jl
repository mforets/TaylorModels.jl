using BenchmarkTools, TaylorModels

SUITE = BenchmarkGroup()  # parent BenchmarkGroup to contain our suite

include("daisy.jl")
