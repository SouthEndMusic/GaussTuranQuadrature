using Test
using GaussTuranQuadrature
using TaylorDiff

@testset "GaussTuranRule" begin

    n = 3
    s = 2

    I = GaussTuranRule(n, s)
    f = TaylorDiffIntegrand(sin)

    evaluation = I(f)
    @test evaluation ≈ (1 - cos(1)) atol=1e-15

    I = GaussTuranRule(n, s; domain = (0, π/2))
    evaluation = I(f)
    @test evaluation ≈ (1 - cos(π/2)) atol=1e-15
end
