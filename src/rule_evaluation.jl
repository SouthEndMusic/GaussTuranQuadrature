struct GaussTuranRule{WType <: AbstractMatrix, XType <: AbstractVector, xType <: AbstractFloat}
    a::xType
    b::xType
    W::WType # (n_nodes, n_derivs)
    X::XType
    function GaussTuranRule(W::WType, X::XType; domain::Tuple = (0,1)) where {WType, XType}
        new{WType, XType, eltype(X)}(domain[1], domain[2], W, X)
    end
end

n_derivs(I::GaussTuranRule) = size(I.W)[2]

function (I::GaussTuranRule)(f)
    (; a, b, W, X) = I
    l = b - a
    out = zero(T)
    for (i, x) in enumerate(X)
        derivs = derivatives(f, a + l * x, l, Val(n_derivs(I))).value
        for (m, deriv) in enumerate(derivs)
            out += W[i, m] * deriv
        end
    end
    out
end