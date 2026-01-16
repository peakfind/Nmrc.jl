# Finite Element Method with Transparent Boundary Condition

To handle unbounded computational domains, we introduce Dirichlet-to-Neumann 
map (DtN) which maps the Dirichlet value to the Neumann value of the solution 
of a boundary value problem. In our problem, we construct the DtN map based on 
the Rayleigh expansion radiation conditon. Then the DtN map give rise to a transparent 
boundary conditon on the artificial boundary which truncates the unbounded domain.

## DtN map and transparent boundary condition

The DtN map is defined on the exterior domain ``G_{b} = \\{(x_{1}, x_{2}) \\in \\mathbb{R}^{2} : x_{1} \\in (0, p) \\text{ and } x_{2} \\in (b, +\\infty) \\}`` 
with boundary ``\\Gamma_{b} = \\{(x_{1}, x_{2}) \\in \\mathbb{R}^{2} : x_{2} = b\\}``. The mapping ``T`` 
takes Dirichlet value ``f \\in C^{\\infty}`` to the corresponding Neumann value ``Tf = \\frac{\\partial v}{\\partial n}\\vert_{\\Gamma_{b}}``, 
where ``v`` is the solution to the Helmholtz equation in ``G_{b}`` satisfying the Dirichlet 
boundary condition ``v = f`` on ``\\Gamma_{b}`` and the Rayleigh expansion radiation condition. We 
conclude the following boundary value problem
```math
\begin{equation*}
    \begin{cases}
        \Delta v + k^{2}v = 0, \\
        v = f, \\
    \end{cases}
\end{equation*}
```

## Finite element discretization