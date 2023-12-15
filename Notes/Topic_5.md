# Numerical Solutions of Ordinary Differential Equations (I): One-Step Methods

## The Initial Value Problem (IVP)

The type of problems we will solve in this topic are of the form:

> Find a real-valued function $y\in C^1(I)$ ( where $I$ is an interval in $\mathbb{R}$) such that:
> $$y'(t) = f(t, y(t)) , \ t\in I$$
> $$y(t_0) = y_0$$
> where $f$ is  a given real-valued function defined on the strip $S = I \times (−\infty, +\infty)$, which is continuous
with respect to both variables.

The initial condition enables us to select one of the infinite solutions that are possible, and which differ in only an additive constant. 

If $f$ is continuous, then the *Cauchy* problem is equivalent to:

$$y(t) - y_0 = \int_{t_0}^t f(\tau, y(\tau))d\tau$$


## Well Posedness

> The IVP
> $$y'(t) = f(t, y),\ a\leq t\leq b,\ y(a) = y_0$$
> is **well posed** if:
> * There exists a unique solution of the problem $y(t)$.
> * There exist constants $\epsilon_0 > 0$ and $k > 0$ such that for any $\epsilon$, with $\epsilon_0 > \epsilon > 0$, whenever $\delta (t)$ is continuous with $|\delta (t)| < \epsilon$ for all $t$ in $[a, b]$, and when $|\delta_0 | < \epsilon$, the initial-value problem
>   $$z'(t) = f(t, z) + \delta(t),\ a \leq t \leq b, \ z(a) = y_0 + \delta_0,$$
>   has a unique solution $z(t)$ that satisfies
> $$|y(t) -z(t)|<k\epsilon, \ \text{for all }t\in [a, b].$$

In order to develop this definition of *well posedness*, we need define the concept of *Lipschitz Continuity*:

> **Lipschitz Continuity**\
> We say that $f(t, y)$ is locally **Lipschitz continuous** at $(t_0 , y_0 )$ with respect to $y$, if there exist a neighborhood $J \subseteq I$ of $t_0$ with width $r_J$, a neighbourhood $\Sigma$ of $y_0$ with width $r_\Sigma$, and a constant $L > 0$, such that
> $$|f(t, y_1) -f(t, y_2)|\leq L|y_1-y_2|, \ \forall t\in J \text{ and }\forall y_1, y_2 \in \Sigma.$$
> Note: Being *Lipschitz continuous* is slightly more than being continuous but slightly less than being differentiable.

Then, the definition of well posedness can be rewritten:

> **Well Posedness**\
> Suppose that $D = \{ (t, y): a \leq t\leq b\text{ and }-\infty < y < \infty \}$ and that $f(t, y)$ is continuous on $D$. If $f$ satisfies a Lipschitz continuity condition on $D$ in the variable $y$, then the IVP
> $$y'(t) = f(t, y),\ a\leq t\leq b,\ y(a) = y_0$$
> is well posed in $D$.

## One-step Numerical Methods

In this section, we will find some methods to help us in solving the following IVP:
$$y'(t) = f(t, y(t)) , \ t\in [t_0, t_0 + T]$$
$$y(t_0) = y_0$$

Specifically, we will focus on *one-step numerical methods*:

> A numerical method for the approximation of the Cauchy problem is called a **one-step method** if $\forall i \geq 0$, $u_{i+1}$ depends only on $u_i$. Otherwise, the scheme is called a multistep method.



### Euler's Methods

There are two variations of **Euler's Method**, the *implicit* and the *explicit* Euler methods. 

> A method is called *explicit* (or *forward*) if $u_{i+1}$ can be computed directly in terms of (some of) the previous values $u_k$, $k \leq i$. A method is said to be *implicit* (or *backward*) if $u_{i+1}$ depends implicitly on itself through $f$.

Thus, in *implicit* methods, we often need to solve a root-finding problem at each step, for a certain equation with $u_{i+1}$ as the unknown.

For the *explicit* method:

> **Euler's Explicit Method** constructs the solution of the IVP with the following iterative formula:
> $$u_{i+1} = u_i  + h f_i, \text{ for }i = 0, 1, ..., N$$
> where $f_i = f(t_i, u_i)$, $h$ is the step size and the initial condition is $u_0 = y(t_0) = y_0$.

For the *implicit* method:

> **Euler's Implicit Method** constructs the solution of the IVP with the following iterative formula:
> $$u_{i+1} = u_i  + h f_{i+1}, \text{ for }i = 0, 1, ..., N$$
> where $f_{i+1} = f(t_{i+1}, u_{i+1})$, $h$ is the step size and the initial condition is $u_0 = y(t_0) = y_0$.

### Crank-Nicolson Method

The **Crank-Nicolson Method** (or *Trapezoidal Method*)

> The **Crank-Nicolson Method** constructs the solution of the IVP with the following iterative formula:
> $$u_{i+1} = u_i  + \frac{h}{2}(f_i + f_{i+1}), \text{ for }i = 0, 1, ..., N$$
> where $f_k = f(t_k, u_k)$, $h$ is the step size and the initial condition is $u_0 = y(t_0) = y_0$.

As the expression for $u_{i+1}$ includes $f_{i+1}$, we can see that the *Crank Nicolson* method is *implicit*.

### Heun's Method

**Heun's Method** comes from attempting to make the *Cranck-Nicolson* method explicit by replacing $u_{i+1}$ by $u_i+hf_i$:

> **Heun's Method** constructs the solution of the IVP with the following iterative formula:
> $$u_{i+1} = u_i  + \frac{h}{2}(f_i + f(t_{i+1}, u_i+hf_i)), \text{ for }i = 0, 1, ..., N$$
> where $f_i = f(t_i, u_i)$, $h$ is the step size and the initial condition is $u_0 = y(t_0) = y_0$.

## Analysis of One-Step Methods

> The general expression of all one-step methods is:
> $$u_{n + 1} = u_n +h\Phi(t_n, u_n, f_n; h), \ 0\leq n\leq N_h-1, \ u_0 = y_0 $$
> where $\Phi(t_n, u_n, f_n; h)$ is called the increment function, which fully characterizes the method.

When we use the analytical solution of the IVP in the general expression of a one-step method in order to obtain the next numerical value of in the solution, we obtain:

$$y_{n+1} = y_n +h\Phi(t_n, y_n, f(t_n, y_n); h) + \varepsilon_{n+1}$$

where $\varepsilon_{n+1}$ is the residual (the difference between the real value of the solution at the next point and the one obtained from the given one-step method), which we can rewrite as:

$$
\varepsilon_{n+1} = h\tau_{n +1}(h)
$$

and $\tau_{n +1}(h)$ is the *local truncation error* at $t_{n+1}$

> The **Global Truncation Error** is defined as:
> $$ \tau(h) = \max_{0\leq n\leq N_h-1} |\tau_{n +1}(h)|$$

### Consistency

> A method is **consistent** if its increment function fulfills:
> $$\lim_{h \to 0}\Phi(t_n, y_n, f(t_n, y_n); h) = f(t_n, y_n), \ \forall t_n \geq t_0,$$
> which means that
> $$\lim_{h \to 0}\tau (h) = 0.$$
> Moreover, the method has *order of consistency* $p$ if the solution fulfills $\tau(h) = O(h^p)$ for $h\to0$.

### Zero-Stability

> A method is **zero-stable** if its increment function $\Phi$ is Lipschitz Continuous with respect to the second argument, with constant $\Lambda$ independent of $h$ and of the nodes $t_i \in[t_0, t_0 + T]$.

### Convergence

> A method is **convergent** if it is *zero-stable* and *consistent*.

### Truncation Error

When studying the error of a method, it may seem that in order to get better solutions we should just make $h$ smaller and smaller. However, this is not true, because there comes a point when the computer's round-off errors start to be significant. This should be taken into account when choosing $h$ for a method.

### Absolute Stability

Informally, **absolute stability** in a numerical method means that, for a fixed $h$, $u_n$ remains bounded as $t_n \to +\infty$. This is different from **zero-stability**, as zero-stability means that for a finite interval, $u_n$ remains bounded when $h \to 0$.

> The **test problem** is the following IVP:
> $$\begin{cases}
    y'(t) = \lambda y(t), t > 0,\\
    y(0) = 1,
\end{cases}$$
> where $\lambda\in \mathbb{C}$. The solution to this IVP is $e^{\lambda t}$. Notice that $\lim_{t\to\infty}y(t) = 0$ if $Re(\lambda) < 0$.

The formal definition of absolute stability is the following:

> A numerical method is **absolutely stable** if on approximating test problem (with $Re(\lambda) < 0$ and for a fixed $h\lambda$) it verifies:
> $$|u_n| \to 0 \text{ as }t_n \to +\infty$$

According to this, we can define the region of absolute stability:

>  The **region of absolute stability** of a numerical method for the test problem is the subset of the complex plane
> $$A = \{z = h\lambda \in C \text{ for which the test problem is absolutely stable}\}$$

The regions of absolute stability for the methods we have studied are:

* Explicit Euler: $|1 + h\lambda| < 1$
* Implicit Euler: $|1 + h\lambda| > 1$
* Crank-Nicolson: $Re(\lambda) < 0$
* Heun: $∣1 + h\lambda +h  \lambda^2∣ < 1$

### $A$-Stability

>  We say that a method is **$A$-stable** if $A ∩ \mathbb{C}^- =\mathbb{C}^-$, i.e. if for $Re(\lambda) < 0$, the condition of absolute stability is satisfied for all values of $h$ (the region of absolute stability includes all the negative complex plane). Otherwise, the method is said to be conditionally stable.

**Note:** there are no explicit unconditionally absolutely stable schemes.

### Summary

In summary, the analysis of these methods includes:

$$
\text{Analysis of Methods}
\begin{cases}
\text{Convergence: }|u_n - y_n|\to 0\text{ as }h\to 0
\begin{cases}
\text{Consistency: }\tau_{n + 1}(h)\to 0\text{ as }h\to 0\\
\text{Zero Stability: }\Phi\text{ Lipschitz Continuous in }u_n
\end{cases}\\
\text{Absolute Stability: }|u_n|\to 0\text{ as }t_n\to\infty\text{ for the test problem}
\end{cases}
$$