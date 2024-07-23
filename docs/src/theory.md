# Gauss-Turán quadrature rules

A Gauss-Turán quadrature rule is an operation on a function $f : (a,b) \rightarrow \mathbb{R}$ of the form

$$
    I(f) = \sum_{m=1}^{2s + 1} \sum_{i = 1}^n w_{m,i}f^{(m-1)}(x_i),
$$

defined by the weights $w_{m,i} \in \mathbb{R}$ and the $n$ nodes $x_i \in (a, b)$, where for the rest of this section we will assume that $a = 0$ and $b = 1$. These variables are chosen such that for a set of linearly independent functions

$$
    \varphi_1, \ldots \varphi_{2(s+1)n} \in C^{2s+1}(0,1)
$$

(i.e. $(2s+1)$ times continuously differentiable on $(a,b)$) the operator $I$ gives an exact value for the integral over the interval $(a,b)$:

$$
    I(\varphi_j) = \int_0^1\varphi_j(x)\text{d}x, \quad j=1,\ldots,2(s+1)n.
$$

By linearity, $I$ will also given exact integrals for linear combinations of the functions $\varphi_j$. For functions $f$ that are well approximated by linear combinations of the $\varphi_j$, $I$ will give a good approximation of the integral:

$$
    I(f) \approx \int_0^1f(x)\text{d}x.
$$

The value $s \in \mathbb{N}$ determines the highest order derivative ($2s$) that is required to evaluate the quadrature rule. Gaussian quadrature uses $s = 0$ and orthogonal polynomials $\varphi_j$.

## Generalizing to any interval $(a,b)$

Say we have a quadrature rule as defined in the previous section, but now we want to use it to integrate a function over the interval $(a,b)$. By integral transformation rules we have that


$$
\begin{align*}
    \int_a^b f(t)\text{d}t &=& (b - a) \int_0^1 f(a + (b-a)x)\text{d}x \\
    &\approx& (b-a)I(\tilde{f})
\end{align*}
$$

where $\tilde{f}(x) = f(a + (b-a)x)$. We can interpret $f \mapsto (b-a)I(\tilde{f})$ as a quadrature rule with 
- weights multiplied by $b-a$;
- transformed nodes $a + (b-a)x_i$;
- derivatives in the "direction" $b - a$.

## Deriving Gauss-Turán quadrature rules

### Deriving a loss function

For given functions $\varphi_j$, (_reference to exact integrals_) defines $2(n+1)n$ equations in the weights and nodes of the quadrature rule. These give $(2s + 1)n + n = 2(s + 1)n$ degrees of freedom, and thus we have a system of equations which in principle we can solve. Due to the $\varphi_j$ this is a nonlinear problem, and so we derive a numerical method which will converge to the quadrature rule.

We reformulate (_reference to exact integrals_) as a linear system

$$
    M(\mathbf{x})\mathbf{w} = \mathbf{b},
$$

where

- $M(\mathbf{x})$ is a $2(s+1)n \times (2s + 1)n$ matrix (note that it is not square!) consisting of the evaluations of the derivatives of the $\varphi_j$:

$$
M(\mathbf{x})_{j,k(m, i)} = \varphi_j^{(m-1)}(x_i);
$$

- $\mathbf{w} = (w_{k(m, i)})$ is a column vector of length $(2s+1)n$;
- $\mathbf{b}$ is a column vector of length $2(s+1)n$ consisting of the integrals of the $\varphi_j$:

$$
    b_j = \int_a^b \varphi_j(x)\text{d}x.
$$

here $k(m,i)$ is an unrolling of the double sum in (_reference to quadrature rule definition_) into a single dimension.

For given nodes $\mathbf{x}$ we can solve the first $(2s+1)n$ equations of this linear system for the weights $\mathbf{x}$. To this end we split the matrix $M(\mathbf{x})$ into an upper part $M_\text{upper}(\mathbf{x})$ of (square) size $(2s+1)n \times (2s+1)n$ and a lower part of $M_\text{lower}$ of size $n \times (2s + 1)n$. The right hand side vector $\mathbf{b}$ is split analogously into $\mathbf{b}_\text{upper}$ and $\mathbf{b}_\text{lower}$. This gives

$$
    \mathbf{w} = M_\text{upper}(\mathbf{x})^{-1}\mathbf{b}_\text{upper},
$$

where in the implementation this will be solved as a linear system instead of computing the inverse $M_\text{upper}(\mathbf{x})^{-1}$ explicitly.

Given these weights, we can check the validity of the last $n$ equations, which yields the loss function

$$
    \text{loss}(\mathbf{x}) = \|M_\text{lower}(\mathbf{x}) - \mathbf{b}_\text{lower}\|_2 
    = 
    \|M_\text{lower}(\mathbf{x})M_\text{upper}(\mathbf{x})^{-1}\mathbf{b}_\text{upper} - \mathbf{b}_\text{lower} \|_2,
$$

fully expressed in terms of the nodes $\mathbf{x}$.

### Constraining the optimization problem

We want to solve for the quadrature rule by approximating the root of the above loss function with a Newton solver. However, a naive implementation of this can lead to numerical problems, because the different nodes can converge towards eachother which makes the matrix $M_\text{lower}(\mathbf{x})$ singular. Furthermore, we want the nodes to stay within the interval $(0,1)$.

To this end we define the new variables 

$$
    \Delta \mathbf{x} = (\Delta x_1, \Delta x_2, \ldots, \Delta x_n) = (x_1, x_2 - x_1, \ldots, x_n - x_{n-1})
$$

which implies

$$
    x_i = \sum_{i' = 1}^i \Delta x_{i'}
$$

and a new loss function $\text{Loss}(\Delta\mathbf{x})$ in terms of these variables. Given a minimal distance between the nodes $0 <\varepsilon \le \frac{1}{n+1}$, we define the constraints

$$
    \varepsilon \le \Delta x_{i'} \le 1-2\varepsilon, \quad i' = 1, \ldots, n,
$$

and the final constraint

$$
    n\varepsilon\le \sum_{i'=1}^n \Delta x_{i'} \le 1-\varepsilon.
$$
