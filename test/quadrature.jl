using GaussTuranQuadrature
using Test
# using QuadGK

@testset "GaussTuranRule" begin

    n = 3
    s = 2

    function f(x)
        val = sin(x)
        der = sqrt(1 - val^2) # cos
        (val, der, -val, -der, val)
    end
    for (a, b) in [(0.0, 1.0), (0.0, pi/2), (-0.1, pi/2)]
        I = GaussTuranRule(n, s; domain=(a, b))
        @test I(f) ≈ (-cos(b) - -cos(a)) atol=1e-15
    end
end
