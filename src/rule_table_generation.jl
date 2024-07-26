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
    dict = "rules = Dict{Tuple{Int, Int}, @NamedTuple{X::Vector{Float64}, W::Matrix{Float64}}}()"
    write(file, join([intro1, intro2, dict], "\n"))
end

for n in N
    for s in S
        rhs = 1 ./ collect(range(one(T), T(2 * (s + 1) * n + 1)))

        ϕ = if s == 0
            (x, j) -> (x^(j - 1),)
        else
            (x, j) -> derivatives(x -> x^(j - 1), x, one(T), Val(2s + 1)).value
        end

        I, res = GaussTuranComputeRule(ϕ,
            n,
            s,
            rhs,
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
                "\n\n rules[($n,$s)] = (; X = $(Float64.(I.X)), W = $(Float64.(I.W)))")
        end
        format(filename)
        println("Completed rule for n = $n, s = $s")
    end
end
