# Computing Gauss-Turán qudrature rules

In this section we look at how to use the algorithm in `GaussTuranQuadrature.jl` to compute new quadrature rules. For the mathematical background of this algorithm see the [theory section](theory.md).
We use the example presented in [this article](https://www.sciencedirect.com/science/article/pii/S0898122100850014?via%3Dihub), with the functions

```math
\begin{equation}
    \varphi_j(x) = x^{p_j}, \quad j = 1, \ldots, 2(s + 1)n
\end{equation}
```

where

```math
\begin{equation}
    p_j =
    \begin{cases}
        \frac{j}{2} - 1 \text{ if $j$ is even} \\
        \frac{j - 1}{2} - \frac{1}{3} \text{ if $j$ is odd}
    \end{cases},
\end{equation}
```

and the trivial weight function $w(x) = 1$.

## Computing the rule

Below is some setup of the inputs for computing the rule.

```@example 1
using TaylorDiff: derivatives

p(n, s) = [j % 2 == 0 ? j / 2 - 1 : (j - 1) / 2 - 1 / 3 for j in 1:(2 * (s + 1) * n)]

n = 5
s = 1

# Powers of the input functions
P = p(n, s)

# Integrals over (0, 1) of the input functions
rhs = @. 1 / (P + 1)

# The functions themselves and their derivatives
ϕ = (x, j) -> derivatives(x -> x^P[j], x, 1.0, Val(2s + 1)).value;
```

Now we can compute the rule.

```@example 1
using Optim
using TaylorDiff
using PreallocationTools
using GaussTuranQuadrature

I, res = GaussTuranComputeRule(ϕ, n, s, rhs)

I
```

And evaluate the rule.

```@example 1
# Create a function which computes exp(x) once and fills a tuple of length 2s+1 with that value
f = x -> (val = Base.exp(x); ntuple(_ -> val, 2s + 1))

evaluation = I(f)
error = abs(evaluation - (Base.exp(1) - 1))
error
```

We can also have a look at the optimization statistics.

```@example 1
res
```

## Obtaining higher accuracy
