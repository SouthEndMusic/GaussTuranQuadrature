"""
    GaussTuranRule(W::AbstractMatrix, X::AbstractVector; domain=(0,1))
    GaussTuranRule(n::Integer, s::Integer; domain=(0,1))

Construct a Gauss-Turán quadrature for the given interval either from known weights `X` and
nodes `X` or for a fixed number of nodes `n` and derivatives `s`.

    (::GaussTuranRule)(f)

Evaluate the integral of the function `f` using a computed Gauss-Turán rule. `f(x)` must return
a tuple or vector of the integrand and its first `2s+1` derivatives.
"""
struct GaussTuranRule{
    WType <: AbstractMatrix, XType <: AbstractVector, xType <: AbstractFloat}
    a::xType
    b::xType
    W::WType # (n_nodes, n_derivs)
    X::XType
    function GaussTuranRule(W::WType,
            X::XType;
            domain::Tuple{<:Number,<:Number} = (0,1)) where {WType <: AbstractMatrix, XType <: AbstractVector}
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

n_derivs(I::GaussTuranRule) = _n_derivs(I.W)
_n_derivs(W) = size(W, 2)

function GaussTuranRule(n::Integer, s::Integer; domain::Tuple{<:Number, <:Number} = (0, 1))
    key = (n, s)
    if haskey(rules, key)
        data = rules[(n, s)]
        Δx = domain[2] - domain[1]
        # Transformed nodes and weights
        X = @. domain[1] + Δx * data.X
        W = copy(data.W)
        for d in 1:(2s+1)
            W[:,d] *= Δx^d
        end
        return GaussTuranRule(W, X; domain)
    else
        error("No tabulated rule for n = $n, s = $s.")
    end
end

function (I::GaussTuranRule)(f::F) where {F}
    (; a, b, W, X) = I
    return evalrule(f, W, X, a, b)
end

function evalrule(f, W, X, a, b)
    # unroll first iteration to get right type of integrand
    x0, rest = Iterators.peel(X)
    derivs0 = f(x0)
    @assert length(derivs0) == _n_derivs(W)
    deriv0, drest = Iterators.peel(derivs0)
    out = W[1,1] * deriv0
    for (_m, d0) in enumerate(drest)
        out += W[1,_m+1] * d0
    end
    # remaining points
    for (_i, x) in enumerate(rest)
        i = _i + 1
        derivs = f(x)
        for (m, deriv) in enumerate(derivs)
            out += W[i, m] * deriv
        end
    end
    return out
end

"""
    TaylorDiffIntegrand(f, order=nothing)

Construct an integrand from `f` for use with [`GaussTuranRule`](@ref) having derivatives
computed automatically using TaylorDiff.jl.
The second argument `order` should be `Val(2s+1)`, where `s` is the derivative order of the Gauss-Turán rule.
This functionality is implemented in a package extension requiring `using TaylorDiff`
"""
struct TaylorDiffIntegrand{F,N}
    f::F
    order::N
end
TaylorDiffIntegrand(f) = TaylorDiffIntegrand(f, nothing)
