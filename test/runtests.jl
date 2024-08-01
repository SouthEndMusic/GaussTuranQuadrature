using SafeTestsets

@safetestset "aqua" begin include("aqua.jl") end
@safetestset "accuracy" begin include("accuracy.jl") end
@safetestset "quadrature" begin include("quadrature.jl") end
@safetestset "TaylorDiff" begin include("taylordiff.jl") end
@safetestset "new rules" begin include("newrules.jl") end
