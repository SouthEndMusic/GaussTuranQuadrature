module GaussTuranQuadrature

using TaylorDiff: derivatives

# Define function name for GaussTuranQuadratureComputeRulesExt.jl
function GaussTuranComputeRule end

include("rule_setup.jl")
include("rule_evaluation.jl")

export GaussTuranRule, GaussTuranComputeRule

end # module GaussTuranQuadrature