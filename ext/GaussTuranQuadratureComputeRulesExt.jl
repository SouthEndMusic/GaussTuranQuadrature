module GaussTuranQuadratureComputeRulesExt
using GaussTuranQuadrature
using Base.Threads
if isdefined(Base, :get_extension)
    using Optim
    using TaylorDiff
    using PreallocationTools
else
    using ..Optim
    using ..TaylorDiff
    using ..PreallocationTools
end

function DEFAULT_w(x::T)::T where {T}
    one(T)
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
    A_buffer::LazyBufferCache{typeof(identity)}
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
    A_buffer,
    X_buffer
) = cache
    M_upper = M_upper_buffer[ΔX, (N, N)]
    M_lower = M_lower_buffer[ΔX, (n, N)]
    A = A_buffer[ΔX, N]
    X = X_buffer[ΔX, n]

    # Compute X from ΔX
    cumsum!(X, ΔX)

    # Evaluating ϕ derivatives
    for (i, x) in enumerate(X)
        Threads.@threads for j in 1:N
            M_upper[j, i:n:N] .= derivatives(x -> ϕ(x, j), x, 1, Val(2s + 1)).value
        end
        Threads.@threads for j in (N + 1):(N + n)
            M_lower[j - N, i:n:N] .= derivatives(x -> ϕ(x, j), x, 1, Val(2s + 1)).value
        end
    end

    # Solving for A
    A .= M_upper \ rhs_upper

    # Computing output
    out = zero(eltype(ΔX))
    for i in eachindex(X)
        out_term = -rhs_lower[i]
        for j in eachindex(A)
            out_term += A[j] * M_lower[i, j]
        end
        out += out_term^2
    end
    sqrt(out)
end

"""
    Callable result object of the Gauss-Turán quadrature rule
    computation algorithm.
"""
struct GaussTuranResult{T, RType, dϕType}
    X::Vector{T}
    A::Matrix{T}
    res::RType
    cache::GaussTuranCache
    dϕ::dϕType

    function GaussTuranResult(res, cache::GaussTuranCache{T}, dϕ) where {T}
        (; A_buffer, s, n, N) = cache
        X = cumsum(res.minimizer)
        dϕ.f(res.minimizer)
        A = reshape(A_buffer[T[], N], (n, 2s + 1))
        new{T, typeof(res), typeof(dϕ)}(X, A, res, cache, dϕ)
    end
end

function GaussTuranRule(res, cache::GaussTuranCache{T}, dϕ) where {T}
    (; A_buffer, s, n, N) = cache
    X = cumsum(res.minimizer)
    dϕ.f(res.minimizer)
    A = reshape(A_buffer[T[], N], (n, 2s + 1))
    GaussTuranQuadrature.GaussTuranRule(A, X)
end

"""
    Input: function f(x, d) which gives the dth derivative of f
"""
function (I::GaussTuranResult{T} where {T})(integrand)
    (; X, A, cache) = I
    (; s) = cache
    out = zero(eltype(X))
    for (i, x) in enumerate(X)
        derivs = derivatives(integrand, x, 1.0, Val(2s + 1)).value
        for (m, deriv) in enumerate(derivs)
            out += A[i, m] * deriv
        end
    end
    out
end

"""
    GaussTuran(ϕ, n, s; w = DEFAULT_w, ε = nothing, X₀ = nothing, integration_kwargs::NamedTuple = (;), optimization_options::Optim.Options = Optim.Options(), T::Type{<:AbstractFloat} = Float64)

Compute a Gauss-Turán quadrature rule by solving a constrained non-linear optimization problem.
For details about the method see (_reference to theory_).

## Inputs

  - `ϕ`: Function with signature `ϕ(x::T, j)::T` that returns ϕⱼ at x
  - `n`: The number of nodes in the quadrature rule
  - `s`: Determines the highest order derivative required from the functions ϕⱼ, currently 2(s + 1)
  - `rhs`: The integrals of the functions ϕⱼ weighted by `w`. The element type of this vector (<: AbstractFloat)
    deterimnes the float type used for the computations. Use e.g. `DoubleFloats.Double64`` or `BigFloat` (slow!) for high accuracy results.

## Keyword Arguments

  - `w`: the integrand weighting function, must have signature w(x::Number)::Number. Defaults to `w(x) = 1`.
  - `ε`: the minimum distance between nodes. Defaults to 1e-3 * (b - a) / (n + 1).
  - `X₀`: The initial guess for the nodes. Defaults to uniformly distributed over (a, b).
  - `integration_kwargs`: The key word arguments passed to `solve` for integrating w * fⱼ
  - `optimization_kwargs`: The key word arguments passed to `Optim.Options` for the minization problem
    for finding X.
"""
function GaussTuranQuadrature.GaussTuranComputeRule(
        ϕ,
        n::Integer,
        s::Integer,
        rhs::AbstractVector{T};
        w = DEFAULT_w,
        ε = nothing,
        X₀ = nothing,
        integration_kwargs::NamedTuple = (;),
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

    # Integrate w * ϕ
    integrand = (out, x, j) -> out[] = w(x) * ϕ(x, j)
    function integrate(j)
        prob = IntegralProblem{true}(integrand, (zero(T), one(T)), j)
        res = solve(prob, QuadGKJL(); integration_kwargs...)
        res.u[]
    end
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
