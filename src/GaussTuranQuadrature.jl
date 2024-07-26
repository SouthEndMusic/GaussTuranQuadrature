module GaussTuranQuadrature

# Define function name for GaussTuranQuadratureComputeRulesExt.jl
function GaussTuranComputeRule end

include("rules.jl")
include("rule_evaluation.jl")

export GaussTuranRule, GaussTuranComputeRule

end # module GaussTuranQuadrature
