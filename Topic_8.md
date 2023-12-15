# Numerical Solutions of Partial Differential Equations

## Types of Partial Differential Equations

In this topic, we will focus on *linear second order partial differential equations with constant coefficients*:

> In general, any **linear second order partial differential equation with constant coefficients** can be expressed as:
> $$A\frac{\partial^2 u}{\partial x^2} + B \frac{\partial^2u}{\partial x\partial y} + C \frac{\partial^2u}{\partial y^2} + D\frac{\partial u}{\partial x} + E \frac{\partial u}{\partial y} + Fu = G(x, y)$$
> where $u(x, y)$ is the unknown function.



There are three fundamental kinds of partial differential equations: *parabolic*, *hyperbolic* and *elliptic*, depending on the sign of a parameter we define as: $\Delta = B^2 - 4AC$.

> A **Hyperbolic PDE** is a linear second order partial differential equation such that $\Delta > 0$. An example of one such PDE is the *Wave Equation*:
> $$u_{tt} - c^2 u_{xx} = 0$$

> A **Parabolic PDE** is a linear second order partial differential equation such that $\Delta = 0$. An example of one such PDE is the *Heat Equation*:
> $$u_{t} - c u_{xx} = 0,\ c>0$$

> An **Elliptic PDE** is a linear second order partial differential equation such that $\Delta < 0$. An example of one such PDE is the *Laplace Equation*:
> $$u_{xx} + u_{yy} = 0$$

## Finite Differences

The simplest method of solving PDEs numerically is by the use of finite differences. We can use either backward, centered or forward finite differences. Typically, the kind of finite difference we use for the first derivatives depends on the problem. If we denote $u_{ij} \equiv u(x_i, x_j)$, where $x_{i+1} = x_i + h$ and $y_{j+1} = y_j + k$:

> The *Forward Finite Differences* for the first derivative in the $x$ and $y$ coordinates, respectively, are:
> $$\begin{align*}
    u_{x}(x_i, y_i)&\approx\frac{u_{i+1,j}-u_{i,j}}{h} & u_{y}(x_i, y_i)&\approx\frac{u_{i,j+1}-u_{i,j}}{k}
\end{align*}$$

> The *Backward Finite Differences* for the first derivative in the $x$ and $y$ coordinates, respectively, are:
> $$\begin{align*}
    u_{x}(x_i, y_i)&\approx\frac{u_{i,j}-u_{i-1,j}}{h} & u_{y}(x_i, y_i)&\approx\frac{u_{i,j}-u_{i,j-1}}{k}
\end{align*}$$

> The *Centered Finite Differences* for the first derivative in the $x$ and $y$ coordinates, respectively, are the average of the backward and forward formulas:
> $$\begin{align*}
    u_{x}(x_i, y_i)&\approx\frac{u_{i+1,j}-u_{i-1,j}}{2h} & u_{y}(x_i, y_i)&\approx\frac{u_{i,j+1}-u_{i,j-1}}{2k}
\end{align*}$$

For the second $x$ and $y$ derivatives, we normally use centered differences:

> The *Centered Finite Differences* for the second derivative in the $x$ and $y$ coordinates, respectively, are:
> $$\begin{align*}
    u_{xx}(x_i, y_i)&\approx\frac{u_{i+1,j}-2u_{i,j}+u_{i-1,j}}{h^2}\\ u_{yy}(x_i, y_i)&\approx\frac{u_{i,j+1} -2u_{i,j}+u_{i,j-1}}{k^2}
\end{align*}$$

## The Hyperbolic Problem

> The **wave equation** is defined as
> $$u_{tt} = c^2u_{xx}$$
> in an interval $a<x<b$, for $t>0$. In order to obtain the solution, we need initial conditions for both the position and the velocity of the wave, which are given by functions of $x$:
> $$\begin{align*}
    u(x, 0) &= f(x) & u_t(x, 0) &= g(x) & a\leq x\leq b
\end{align*}$$
> and we also need some boundary conditions at each end, which are given functions of $t$:
> $$\begin{align*}
    u(a, t) &= l(t) & u_t(b, t) &= r(t) & t>0
\end{align*}$$

### Explicit Method
In order to solve this problem, use the centered finite difference formulas for the second derivative:

$$
\frac{u_{i,j+1} -2u_{i,j}+u_{i,j-1}}{k^2} = c^2\frac{u_{i+1,j}-2u_{i,j}+u_{i-1,j}}{h^2}
$$

If we solve for the term at time $t_{j + 1}$, we obtain:

$$
\boxed{u_{i, j + 1} = \alpha^2(u_{i+1,j} + u_{i-1,j})+2(1-\alpha^2)u_{i,j}-u_{i,j-1}} \tag{1}
$$

where $\alpha = \frac{ck}{h}$ is the so called *Courant Parameter*, which is important for studying the convergence of the method using the ***Courant–Friedrichs–Lewy condition***.

> The ***Courant–Friedrichs–Lewy condition*** (*CFL condition*) states that the numerical solution to a partial differential equation converges if and only if the following relation holds:
> $$\alpha = \frac{c\Delta t}{\Delta x}<C_{max}$$
> for some $C_{max}$ that depends on the method used. $\alpha$ is called the *Courant Parameter*.

In this particular case, we have $C_{max} =1$. In order to start the method, we need to find the solution at $t_0$ and $t_1 = t_0+k$. From the initial condition $u(x, 0) = f(x)$, we can get the solution at $t_0$:

$$
\boxed{u_{i,0} = f_i}
$$

On the other hand, from the initial condition $u_t(x, 0) = g(x)$, we obtain (using centered finite diferences for the time derivative):

$$
u_{t}(x_i, t_0)=\frac{u_{i,1}-u_{i,-1}}{2k} = g(x_i)\to u_{i,-1} = u_{i,1}-2kg_i
$$

Then, if we evaluate the expression in $(1)$ at $j = 0$:

$$
u_{i, 1} = \alpha^2(u_{i+1,0} + u_{i-1,0})+2(1-\alpha^2)u_{i,0}-u_{i,-1}=
$$
$$
= \alpha^2(u_{i+1,0} + u_{i-1,0})+2(1-\alpha^2)u_{i,0}-(u_{i,1}-2kg_i)
$$

then:
$$
u_{i, 1}= \alpha^2\frac{(u_{i+1,0} + u_{i-1,0})}{2}+(1-\alpha^2)u_{i,0}+kg_i
$$

and finally:
$$
\boxed{u_{i, 1}= \alpha^2\ \frac{f_{i+1} + f_{i-1}}{2}+(1-\alpha^2)f_{i}+kg_i}
$$


Once we have $u_{i, 0}$ and $u_{i, 1}$, we can keep using the formula in $(1)$ to find the solution at the next time steps.


<center><img src="https://i.imgur.com/bFHhIGX.png" width="300px"></center>

## The Parabolic Problem

> The **heat equation** is defined as
> $$u_{t} = cu_{xx}$$
> in an interval $a<x<b$, for $t>0$. In order to obtain the solution, we need the initial condition:
> $$\begin{align*}
    u(x, 0) &= f(x) & a< x<b
\end{align*}$$
> and we also need some boundary conditions at each end, which are given as:
> $$\begin{align*}
    u(a, t) &= T_a & u_t(b, t) &= T_b & t>0
\end{align*}$$

### Explicit Method

In order to solve this problem, use the forward finite difference formula for the time derivative and the central finite difference formula for the spatial second derivative:

$$
\frac{u_{i,j+1} -u_{i,j}}{k} = c\frac{u_{i+1,j}-2u_{i,j}+u_{i-1,j}}{h^2}
$$

If we solve for the term at time $t_{j + 1}$, we obtain:

$$
\boxed{u_{i, j + 1} = \alpha(u_{i+1,j} + u_{i-1,j})+(1-2\alpha)u_{i,j}} \tag{2}
$$

where $\alpha = \frac{ck}{h^2}$. In this case, the method converges for $\alpha \leq \frac{1}{2}$, and the fastest speed of convergence is achieved for $\alpha = \frac{1}{6}$.

In order to start the method, we need to find the solution at $t_0$. From the initial condition $u(x, 0) = f(x)$:

$$
\boxed{u_{i,0} = f_i}
$$

On the other hand, from the boundary conditions:

$$
\boxed{u_{0,j} = T_a}
$$
$$
\boxed{u_{n,j} = T_b}
$$

where $n+1$ is the total number of spatial points in the mesh.



Once we have $u_{i, 0}$ and $u_{i, 1}$, we can keep using the formula in $(1)$ to find the solution at the next time steps.


<center><img src="https://i.imgur.com/HavykZx.png" width="300px"></center>


## The Elliptic Problem

> The **Laplace equation** is defined as
> $$u_{xx}+u_{yy} =0$$
> in the rectangle $R = \{ (x, y): a<x<b, c<y<d\}$, where we know the value of $u(x, y)$ at the boundary of the rectangle $R$.

In order to solve this problem, use centered finite difference formulas for both derivatives:

$$
\frac{u_{i+1,j}-2u_{i,j}+u_{i-1,j}}{h^2}+\frac{u_{i,j+1} -2u_{i,j}+u_{i,j-1}}{k^2} = 0
$$

This gives the following square linear system:

$$
\boxed{\alpha^2 (u_{i-1,j} + u_{i+1,j}) + u_{i,j-1} + u_{i,j+1}-2(\alpha^2+1)u_{i,j} = 0}\tag{3}
$$

where $\alpha = k / h$. 



<center><img src="https://i.imgur.com/Ez2nZTg.png" width="300px"></center>

The system in $(3)$ can often be very large, so solving it by traditional methods is not generally an option. Instead, we use iterative methods such as *Jacobi*, *Gauss-Seidel* and the *Successive Over-Relaxation Method*.



