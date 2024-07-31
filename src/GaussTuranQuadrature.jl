module GaussTuranQuadrature

# Define function name for GaussTuranQuadratureComputeRulesExt.jl
"""
    GaussTuranComputeRule

This function is implemented in a package extension requiring
`using GaussTuranQuadrature, Optim, PreallocationTools`
"""
function GaussTuranComputeRule end

include("rules.jl")
include("rule_evaluation.jl")

export GaussTuranRule, GaussTuranComputeRule

end # module GaussTuranQuadrature
