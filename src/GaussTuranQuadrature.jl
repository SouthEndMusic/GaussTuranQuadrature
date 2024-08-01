module GaussTuranQuadrature

# Define function name for GaussTuranQuadratureComputeRulesExt.jl
"""
    GaussTuranComputeRule(::Type{T}, n::Int, s::Int; kws...) where {T}

Compute a monomial Gauss-Turán rule of order `(n, s)` and precision `T`.
This function is implemented in a package extension requiring
`using GaussTuranQuadrature, Optim, PreallocationTools`
"""
function GaussTuranComputeRule(T, n, s; kws...)
    error("No tabulated rule for n = $n, s = $s, and $T precision. To compute a rule try `using Optim, PreallocationTools`")
end

include("rules.jl")
include("rule_evaluation.jl")

export GaussTuranRule, GaussTuranComputeRule
export TaylorDiffIntegrand

end # module GaussTuranQuadrature
