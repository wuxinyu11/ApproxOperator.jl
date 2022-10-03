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

Herein, we introduce finite element method to approximate the displacement $u$ and its virtual one $v$, yields
$$
v(x) = \sum_{I=1}^{n_p} N_I(x) \delta d_I, \quad
u(x) = \sum_{I=1}^{n_p} N_I(x) d_I
$$
substituting them into weak form and eliminating $\delta d$ can obtain the discrete equilibrium equations
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
