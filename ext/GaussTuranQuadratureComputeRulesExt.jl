module GaussTuranQuadratureComputeRulesExt
using GaussTuranQuadrature
using Base.Threads
if isdefined(Base, :get_extension)
    using Optim
    using PreallocationTools
else
    using ..Optim
    using ..PreallocationTools
end

"""
Cached data for the `GaussTuranLoss!` call.
"""
struct GaussTuranCache{T}
    n::Int
    s::Int
    N::Int
    ε::T
    rhs_upper::Vector{T}
    rhs_lower::Vector{T}
    M_upper_buffer::LazyBufferCache{typeof(identity)}
    M_lower_buffer::LazyBufferCache{typeof(identity)}
    X_buffer::LazyBufferCache{typeof(identity)}
    W_buffer::LazyBufferCache{typeof(identity)}
    function GaussTuranCache(
            n,
            s,
            N,
            ε::T,
            rhs_upper::Vector{T},
            rhs_lower::Vector{T}
    ) where {T}
        new{T}(
            n,
            s,
            N,
            ε,
            rhs_upper,
            rhs_lower,
            LazyBufferCache(),
            LazyBufferCache(),
            LazyBufferCache(),
            LazyBufferCache()
        )
    end
end

"""
Function whose root defines the quadrature rule.
"""
function GaussTuranLoss!(ϕ, ΔX::AbstractVector{T}, cache) where {T}
    (;
    n,
    s,
    N,
    rhs_upper,
    rhs_lower,
    M_upper_buffer,
    M_lower_buffer,
    W_buffer,
    X_buffer
) = cache
    M_upper = M_upper_buffer[ΔX, (N, N)]
    M_lower = M_lower_buffer[ΔX, (n, N)]
    W = W_buffer[ΔX, N]
    X = X_buffer[ΔX, n]

    # Compute X from ΔX
    cumsum!(X, ΔX)

    # Evaluating ϕ derivatives
    for (i, x) in enumerate(X)
        Threads.@threads for j in 1:N
            M_upper[j, i:n:N] .= ϕ(x, j)
        end
        Threads.@threads for j in (N + 1):(N + n)
            M_lower[j - N, i:n:N] .= ϕ(x, j)
        end
    end

    # Solving for W
    W .= M_upper \ rhs_upper

    # Computing output
    out = zero(eltype(ΔX))
    for i in eachindex(X)
        out_term = -rhs_lower[i]
        for j in eachindex(W)
            out_term += W[j] * M_lower[i, j]
        end
        out += out_term^2
    end
    sqrt(out)
end

function GaussTuranRule(res, cache::GaussTuranCache{T}, dϕ) where {T}
    (; W_buffer, s, n, N) = cache
    X = cumsum(res.minimizer)
    dϕ.f(res.minimizer)
    W = reshape(W_buffer[T[], N], (n, 2s + 1))
    GaussTuranQuadrature.GaussTuranRule(W, X)
end

"""
    GaussTuran(ϕ, n, s; ε = nothing, X₀ = nothing, optimization_options::Optim.Options = Optim.Options())

Compute a Gauss-Turán quadrature rule by solving a constrained non-linear optimization problem.
For details about the method see (_reference to theory_).

## Inputs

  - `ϕ`: Function with signature `ϕ(x::T, j)::T` that returns ϕⱼ(x), ∂¹ϕⱼ(x), …, ∂²ˢϕⱼ(x)
  - `n`: The number of nodes in the quadrature rule
  - `s`: Determines the highest order derivative required from the functions ϕⱼ, currently 2(s + 1)
  - `rhs`: The integrals of the functions ϕⱼ weighted by `w`. The element type of this vector (<: AbstractFloat)
    deterimnes the float type used for the computations. Use e.g. `DoubleFloats.Double64`` or `BigFloat` (slow!) for high accuracy results.

## Keyword Arguments

  - `ε`: the minimum distance between nodes. Defaults to 1e-3 * / (n + 1).
  - `X₀`: The initial guess for the nodes. Defaults to uniformly distributed over (0, 1).
  - `optimization_kwargs`: The key word arguments passed to `Optim.Options` for the minization problem
    for finding X.
"""
function GaussTuranQuadrature.GaussTuranComputeRule(
        ϕ,
        n::Integer,
        s::Integer,
        rhs::AbstractVector{T};
        ε = nothing,
        X₀ = nothing,
        optimization_options::Optim.Options = Optim.Options()
) where {T <: AbstractFloat}
    # Initial guess
    if isnothing(X₀)
        X₀ = collect(range(zero(T), one(T), length = n + 2)[2:(end - 1)])
    else
        @assert length(X₀) == n
    end
    ΔX₀ = diff(X₀)
    pushfirst!(ΔX₀, X₀[1])

    # Minimum distance between nodes
    if isnothing(ε)
        ε = 1e-3 / (n + 1)
    else
        @assert 0 < ε ≤ 1 / (n + 1)
    end
    ε = T(ε)
    N = (2s + 1) * n
    rhs_upper = rhs[1:N]
    rhs_lower = rhs[(N + 1):(N + n)]

    # Solving constrained non linear problem for ΔX, see
    # https://julianlsolvers.github.io/Optim.jl/stable/examples/generated/ipnewton_basics/

    # The cache for evaluating GaussTuranLoss
    cache = GaussTuranCache(n, s, N, ε, rhs_upper, rhs_lower)

    # The function whose root defines the quadrature rule
    # Note: the optimization method requires a Hessian, 
    # which brings the highest order derivative required to 2s + 2
    func(ΔX) = GaussTuranLoss!(ϕ, ΔX, cache)
    dϕ = TwiceDifferentiable(func, ΔX₀; autodiff = :forward)

    # The constraints on ΔX
    ΔX_lb = fill(ε, length(ΔX₀))
    ΔX_ub = fill(one(T) - 2 * ε, length(ΔX₀))

    # Defining the variable and constraints nε ≤ ∑ΔX ≤ 1 - ε
    sum_variable!(c, ΔX) = (c[1] = sum(ΔX); c)
    sum_jacobian!(J, ΔX) = (J[1, :] .= one(eltype(ΔX)); J)
    sum_hessian!(H, ΔX, λ) = nothing
    sum_lb = [n * ε]
    sum_ub = [one(T) - ε]
    constraints = TwiceDifferentiableConstraints(
        sum_variable!,
        sum_jacobian!,
        sum_hessian!,
        ΔX_lb,
        ΔX_ub,
        sum_lb,
        sum_ub
    )

    # Solve for the quadrature rule by minimizing the loss function
    res = Optim.optimize(dϕ, constraints, T.(ΔX₀), IPNewton(), optimization_options)

    # Return the computed rule for the interval (0, 1) and the optimizer stats
    GaussTuranRule(res, cache, dϕ), res
end

end # module IntegralsGaussTuranExt
