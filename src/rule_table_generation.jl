using DoubleFloats
using Optim
using TaylorDiff
using PreallocationTools
using GaussTuranQuadrature
using JuliaFormatter

include("linear_algebra_fix.jl")

filename = "src/rules.jl"
T = Double64
N = 1:10
S = 0:4

open(filename, "w") do file
    intro1 = "# This file was automatically generated with rule_table_generation.jl"
    intro2 = "# (n, s) => (; X, W)"
    dict = "const rules = Dict{Type,Dict}(Float64 => Dict{Tuple{Int, Int}, @NamedTuple{X::Vector{Float64}, W::Matrix{Float64}}}())"
    write(file, join([intro1, intro2, dict], "\n"))
end

for n in N
    for s in S
        I, res = GaussTuranComputeRule(
            T,
            n,
            s,
            optimization_options = Optim.Options(;
                x_abstol = T(1e-250),
                g_tol = T(1e-250),
                show_trace = true,
                show_every = 1_000,
                iterations = 10_000
            )
        )
        open(filename, "a") do file
            write(file,
                "\n\n rules[Float64][($n,$s)] = (; X = $(Float64.(I.X)), W = $(Float64.(I.W)))")
        end
        format(filename)
        println("Completed rule for n = $n, s = $s")
    end
end
