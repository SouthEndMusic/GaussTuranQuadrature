using Documenter
using GaussTuranQuadrature

makedocs(;
    sitename = "GaussTuranQuadrature.jl",
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "Evaluating quadratures" => "evaluate_rule.md",
            "Computing quadratures" => "compute_rule.md"
        ],
        "API" => "api.md",
        "Theory" => "theory.md"
    ]
)
