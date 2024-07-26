# Evaluating Gauss-Turán quadrature rules

`GaussTuranQuadrature` contains precomputed quadrature rules, which can be used as follows.

```@example 1
using GaussTuranQuadrature

n = 3
s = 2

I = GaussTuranRule(n, s)
```

```@example 1
# Define a function with cheap computation of the derivatives
function f(x)
    val = sin(x)
    der = sqrt(1 - val^2) # cos
    (val, der, -val, -der, val)
end

evaluation = I(f)
error = abs(evaluation - (1 - cos(1)))
```

Other finite domains $(a,b)$ can also be used.

```@example 1
I = GaussTuranRule(n, s; domain = (0, π/2))
I(f)
```