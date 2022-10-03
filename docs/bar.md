---
commits:
    20221002@WuJC 
---

# Bar problem

## Galerkin formulation
Consider a 1D bar problem with a strong form

$$
\left \{
\begin{array}{ll}
(EAu_{,x})_{,x} + b = 0 & x \in (0,L) \\
u = g & x = 0 \\
u_{,x} = t & x = L
\end{array}
\right .
$$

where $EA$ is the stiffness coefficient, $b$ is the body force in the whole bar. The displacement $u$ at left side of bar has been enforced, and the right side of bar is applied by a traction $t$. The corresponding total potential energy functional of bar problem is given by
$$
\Pi(u) = \int_0^L \frac{1}{2} EA u_{,x}^2 dx - ut\vert_{x=L} - \int_0^L ub dx
$$

With a variation with potential energy functional, the Galerkin weak form can be expressed by
$$
\delta \Pi(u) = \int_0^L \delta u_{,x} EA u_{,x} dx - \delta ut\vert_{x=L} - \int_0^L \delta ub dx = 0
$$

Hence, the bar problem becomes a **Ritz** problem as
$$
\text{find} \; u \in H^1, \quad a(v,u) = f(v) \quad \forall v \in H^1
$$
where
$$
a(v,u) = \int_0^L v_{,x} EA u_{,x} dx
$$
$$
f(v) = vt\vert_{x=L} + \int_0^L vb dx
$$

where $a$ is a bilinear form, and $f$ is a linear form.

## Approximation
Herein, we introduce **finite element method (FEM)** to approximate the displacement $u$ and its virtual one $v$. Based upon the methodology of FEM, the whole domain of bar problem is divided into a set of cells $(x_I,x_{I+1}), I = 1,2,\dots, n_e$, where $n_e$ is the total number of cells. And the shape functions of FEM are evaluated by parametric coordinate $\xi$ in each cell. For instance, as $x\in[x_I,x_{I+1}]$, the Cartesian coordinate $x$ and the approximant of displacement $u^h$ can be expressed as
$$
x = N_I(\xi) x_I + N_{I+1}(\xi) x_{I+1}, \quad \xi \in [-1,1]
$$
$$
u^h(\xi) = N_I(\xi)d_I + N_{I+1}(\xi)d_{I+1}, \quad \xi \in [-1,1]
$$
in which
$$
N_I(\xi) = \frac{1-\xi}{2}, \quad N_{I+1}(\xi) = \frac{1+\xi}{2}
$$
$N_I$'s and $d_I$'s are FEM shape functions and coefficients at the node $x_I$, and with this definition, we can find the so-call Kronecker delta property,
$$
N_I(x_J) = \delta_{IJ} = \left \{ \begin{array}{ll}
1, & I=J \\
0, & I\ne J
\end{array} \right .
$$
and $u^h(x_I) = d_I$. On the other sight, the approximants of $u^h$ and $v^h$ can be rewritten by a global form
$$
v^h(x) = \sum_{I=1}^{n_p} N_I(x) \delta d_I, \quad
u^h(x) = \sum_{I=1}^{n_p} N_I(x) d_I
$$
where $n_p$ denotes the total number of discrete points. Substituting them into weak form and eliminating $\delta d$ can obtain the discrete equilibrium equations
$$
\boldsymbol{Kd} = \boldsymbol{f}
$$
where $\boldsymbol{K}$ is the stiffness matrix and $\boldsymbol{f}$ is the force vector, and their components are listed as follows
$$
K_{IJ} = \int_0^L N_{I,x} EA N_{J,x} dx
$$
$$
f_I = N_I(L)t + \int_0^L N_I b dx
$$

## Numerical implementation
The integrations in stiffness matrix and force vector are carried out by traditional Gauss quadrature rule in each cells $(x_K,x_{K+1})$.
$$
\begin{split}
K_{IJ} &= \int_{x_K}^{x_{K+1}} N_{I,x} EA N_{J,x} dx \\
&= \int_{-1}^{1} N_{I,x} EA N_{J,x} \frac{dx}{d\xi}d\xi \\
&= \sum_{G=1}^{n_g} N_{I,x}(\xi_G) EA N_{J,x}(\xi_G) J(\xi_G) w_G
\end{split}
$$
where
$$
J(x) = \frac{dx}{d\xi} = \frac{x_{K+1}-x_K}{2} = \frac{L}{2}
$$

In `ApproxOperator.jl`, the above operator is implemented by the method of type `Operator{:∫vₓuₓdx}`

```julia
function (op::Operator{:∫vₓuₓdx})(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    EA = op.EA
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += B[i]*EA*B[j]*𝑤
            end
        end
    end
end
```
in which `𝑤` stands for $J(\xi_G)w_G$. In this operator, the stiffness $EA$ is embedded in operator. Moreover, the linear form can be implemented by standard linear operator listed as follow
```julia
function (op::Operator{:𝑓𝑣})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        u = ξ.u
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] += N[i]*u*𝑤
        end
    end
end
```
the `u` is a general variable that can stands for body force $b$ or traction $t$.
