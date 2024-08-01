using Test
using GaussTuranQuadrature

@testset "monomial rules ($n, $s)" for (n, s) in keys(GaussTuranQuadrature.rules)
    xs, ws = GaussTuranQuadrature.rules[(n, s)]
    for i in 1:2*(s+1)*n
        # integrate x^(i-1) on [0, 1]
        I_f = 1/i
        Q_f = zero(eltype(xs))
        for j in 1:n
            x = xs[j]
            df = x^(i-1)
            for k in 1:2s+1
                w = ws[j,k]
                Q_f += w*df
                df *= (i-k)/x # numerically compute derivative of monomial
            end
        end
        @test I_f ≈ Q_f atol=10eps(one(eltype(xs)))
    end
end
