module GaussTuranQuadratureTaylorDiffExt
using GaussTuranQuadrature: n_derivs, TaylorDiffIntegrand, evalrule
import GaussTuranQuadrature: GaussTuranRule
using TaylorDiff: derivatives, value

function (I::GaussTuranRule)(f::TaylorDiffIntegrand)
    (; a, b, W, X) = I
    _f = f.f
    order = f.order === nothing ? Val(n_derivs(I)) :
            f.order isa Integer ? Val(f.order) :
            f.order
    return evalrule(x -> value(derivatives(_f, x, 1, order)), W, X, a, b)
end
end
