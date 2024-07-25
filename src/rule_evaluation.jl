struct GaussTuranRule{
    WType <: AbstractMatrix, XType <: AbstractVector, xType <: AbstractFloat}
    a::xType
    b::xType
    W::WType # (n_nodes, n_derivs)
    X::XType
    function GaussTuranRule(W::WType, X::XType; domain::Tuple = (0, 1)) where {WType, XType}
        new{WType, XType, eltype(X)}(domain[1], domain[2], W, X)
    end
end

function Base.show(io::IO, I::GaussTuranRule)
    (; a, b, W, X) = I
    deriv_max = n_derivs(I) - 1
    n = length(X)
    println(io,
        "Gauss-Turán quadrature rule on the interval ($a, $b), with $n nodes and derivatives up to order $deriv_max.")
    b = IOBuffer()
    # Get formatted X as string
    show(b, "text/plain", X)
    X_string = String(take!(b))

    # Get formatted W as string
    b = IOBuffer()
    # Get formatted X as string
    show(b, "text/plain", W)
    W_string = String(take!(b))

    println(io, "* Nodes:\n$X_string")
    println(io, "* Weights:\n$W_string")
end

n_derivs(I::GaussTuranRule) = size(I.W)[2]

function (I::GaussTuranRule)(f)
    (; a, b, W, X) = I
    l = b - a
    out = zero(eltype(X))
    for (i, x) in enumerate(X)
        derivs = derivatives(f, a + l * x, l, Val(n_derivs(I))).value
        for (m, deriv) in enumerate(derivs)
            out += W[i, m] * deriv
        end
    end
    out
end
