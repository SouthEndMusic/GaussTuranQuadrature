struct GaussTuranRule{WType, XType} where {WType <: abstractMatrix, XType <: AbstractVector}
    W::WType # (n_nodes, n_derivs)
    X::XType
end

n_derivs(I::GaussTuranRule) = size(I.W)[2]

function (I::GaussTuranRule)(f)
    (; W, X) = I
    T = eltype(X)
    out = zero(T)
    for (i, x) in enumerate(X)
        derivs = derivatives(f, x, one(T), Val(n_derivs(I))).value
        for (m, deriv) in enumerate(derivs)
            out += W[i, m] * deriv
        end
    end
    out
end