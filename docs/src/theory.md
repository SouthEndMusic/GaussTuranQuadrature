# Gauss-Turán quadrature rules

A Gauss-Turán quadrature rule is an operation on a function $f : (a,b) \rightarrow \mathbb{R}$ of the form

```math
\begin{equation}
    I(f) = \sum_{m=1}^{2s + 1} \sum_{i = 1}^n w_{m,i}f^{(m-1)}(x_i),
\end{equation}
```

defined by the weights $w_{m,i} \in \mathbb{R}$ and the $n$ nodes $x_i \in (a, b)$, where for the rest of this section we will assume that $a = 0$ and $b = 1$. These variables are chosen such that for a set of linearly independent functions

```math
\begin{equation}
    \varphi_1, \ldots \varphi_{2(s+1)n} \in C^{2s+1}(0,1)
\end{equation}
```

(i.e. $(2s+1)$ times continuously differentiable on $(0,1)$) the operator $I$ gives an exact value for the integral:

```math
\begin{equation}
    I(\varphi_j) = \int_0^1\varphi_j(x)\text{d}x, \quad j=1,\ldots,2(s+1)n.
\end{equation}
```

By linearity, $I$ will also given exact integrals for linear combinations of the functions $\varphi_j$. For functions $f$ that are well approximated by linear combinations of the $\varphi_j$, $I$ will give a good approximation of the integral:

```math
\begin{equation}
    I(f) \approx \int_0^1f(x)\text{d}x.
\end{equation}
```

The value $s \in \mathbb{N}$ determines the highest order derivative ($2s$) that is required to evaluate the quadrature rule. [Gaussian quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature) uses $s = 0$ and orthogonal polynomials $\varphi_j$.

## Generalizing to any interval $(a,b)$

Say we have a quadrature rule as defined in the previous section, but now we want to use it to integrate a function over the interval $(a,b)$. By integral transformation rules we have that


```math
\begin{align}
    \int_a^b f(t)\text{d}t &=& (b - a) \int_0^1 f(a + (b-a)x)\text{d}x \\
    &\approx& (b-a)I(\tilde{f})
\end{align}
```

where $\tilde{f}(x) = f(a + (b-a)x)$. We can interpret $f \mapsto (b-a)I(\tilde{f})$ as a quadrature rule with 
- weights multiplied by $b-a$;
- transformed nodes $a + (b-a)x_i$;
- derivatives in the "direction" $b - a$.

## Deriving Gauss-Turán quadrature rules

The algorithm presented here is not published anywhere as far as the author is aware. For an example of an older algorithm that uses an analytically derived Jacobian see [here](https://www.sciencedirect.com/science/article/pii/S0898122100850014?via%3Dihub).

### Deriving a loss function

For given functions $\varphi_j$, (3) defines $2(n+1)n$ equations in the weights and nodes of the quadrature rule. These give $(2s + 1)n + n = 2(s + 1)n$ degrees of freedom, and thus we have a system of equations which in principle we can solve. Due to the $\varphi_j$ this is a nonlinear problem, and so we derive a numerical method which will converge to the quadrature rule.

We reformulate (3) as a linear system

```math
\begin{equation}
    M(\mathbf{x})\mathbf{w} = \mathbf{b},
\end{equation}
```

where

- the matrix $M(\mathbf{x})$ is of size $2(s+1)n \times (2s + 1)n$ (note that it is not square!) and consisting of the evaluations of the derivatives of the $\varphi_j$:

```math
\begin{equation}
M(\mathbf{x})_{j,k(m, i)} = \varphi_j^{(m-1)}(x_i);
\end{equation}
```

- the column vector $\mathbf{w} = (w_{k(m, i)})$ is of length $(2s+1)n$;
- the column vector $\mathbf{b}$ is of length $2(s+1)n$ and consisting of the integrals of the $\varphi_j$:

```math
\begin{equation}
    b_j = \int_a^b \varphi_j(x)\text{d}x.
\end{equation}
```

here $k(m,i)$ is an unrolling of the double sum in (1) into a single dimension.

For given nodes $\mathbf{x}$ we can solve the first $(2s+1)n$ equations of this linear system for the weights $\mathbf{x}$. To this end we split the matrix $M(\mathbf{x})$ into an upper part $M_\text{upper}(\mathbf{x})$ of (square) size $(2s+1)n \times (2s+1)n$ and a lower part of $M_\text{lower}$ of size $n \times (2s + 1)n$. The right hand side vector $\mathbf{b}$ is split analogously into $\mathbf{b}_\text{upper}$ and $\mathbf{b}_\text{lower}$. This gives

```math
\begin{equation}
    \mathbf{w} = M_\text{upper}(\mathbf{x})^{-1}\mathbf{b}_\text{upper},
\end{equation}
```

where in the implementation this will be solved as a linear system instead of computing the inverse $M_\text{upper}(\mathbf{x})^{-1}$ explicitly.

Given these weights, we can check the validity of the last $n$ equations, which yields the loss function

```math
\begin{equation}
    \text{loss}(\mathbf{x}) = \|M_\text{lower}(\mathbf{x}) - \mathbf{b}_\text{lower}\|_2 
    = 
    \|M_\text{lower}(\mathbf{x})M_\text{upper}(\mathbf{x})^{-1}\mathbf{b}_\text{upper} - \mathbf{b}_\text{lower} \|_2,
\end{equation}
```

fully expressed in terms of the nodes $\mathbf{x}$.

### Constraining the optimization problem

We want to solve for the quadrature rule by approximating the root of the above loss function with a Newton solver. However, a naive implementation of this can lead to numerical problems, because the different nodes can converge towards eachother which makes the matrix $M_\text{lower}(\mathbf{x})$ singular. Furthermore, we want the nodes to stay within the interval $(0,1)$.

To this end we define the new variables 

```math
\begin{equation}
    \Delta \mathbf{x} = (\Delta x_1, \Delta x_2, \ldots, \Delta x_n) = (x_1, x_2 - x_1, \ldots, x_n - x_{n-1})
\end{equation}
```

which implies

```math
\begin{equation}
    x_i = \sum_{i' = 1}^i \Delta x_{i'}
\end{equation}
```

and a new loss function $\text{Loss}(\Delta\mathbf{x})$ in terms of these variables. Given a minimal distance between the nodes $0 <\varepsilon \le \frac{1}{n+1}$, we define the constraints

```math
\begin{equation}
    \varepsilon \le \Delta x_{i'} \le 1-2\varepsilon, \quad i' = 1, \ldots, n,
\end{equation}
```

and the final constraint

```math
\begin{equation}
    n\varepsilon\le \sum_{i'=1}^n \Delta x_{i'} \le 1-\varepsilon.
\end{equation}
```

!!! note
    Since we solve this problem with `Optim.IPNewton` which requires a hessian of the loss function, we actually need the functions $\varphi_j$ to be $2s + 2$ times continuously differentiable.