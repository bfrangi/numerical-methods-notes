# Numerical Solutions of Ordinary Differential Equations (II): Multi-Step Methods


## Multi-Step Methods
A multi-step method is defined as:

> A $q$-step method ($q \geq 1$) is one which, $\forall n \geq q - 1$,  $u_{n+1}$ depends on $u_{n+1}-q$, but not on the values $u_k$ with $k < n + 1 - q$.

Some examples include the *Midpoint Method* and the *Simpson Method*.

> A general $(p+1)$-step method can be expressed as:
> $$u_{n + 1} = \sum_{j = 0}^p a_j u_{n - j} + h\sum_{j = 0}^p b_jf_{n - j} + hb_{-1}f_{n + 1}\ \text{ for }n = p, p+1, ...$$
> which means the we need to know the values of the first $p + 1$ steps $(u_0, u_1 ,..., u_p)$ to start with, which we can get using a one-step method.

The multistep method is completely described by the real coefficients $a_j$ and $b_j$ with $a_p$ or $b_p$ different from 0. Also, if $b_{-1} = 0$, the method is *explicit*, and *implicit* otherwise.

## Adams Methods

Adams methods are derived from the integral form of the Cauchy problem (or IVP):

$$
y(t) - y_0 = \int_{t_0}^tf(\tau, y(\tau))d\tau
$$

using a numerical approximation of the integral between $t_n$ and $t_{n+1}$, using also $u_{n-1}$, $u_{n-2}$, ... . This approximation is obtained supposing equispaced nodes, and integrating integrate the interpolating polynomial of $f$ on $p + 1$ distinct nodes. The resulting schemes are consistent by construction and have the following form:

$$
u_{n + 1} = u_n + h \sum_{j = - 1}^p b_j f_{n - j}, \ n \geq p
$$

The interpolation nodes can be either:

* $t_n$, $t_{n-1}$, ... , $t_{n-p}$ (in this case $b_{-1} = 0$ and the resulting method is **explicit**).
* $t_{n+1}$, $t_n$, ..., $t_{n-p+1}$ (in this case $b_{-1} \neq 0$ and the scheme is **implicit**).

The **implicit** schemes are called **Adams-Moulton methods**, while the **explicit** ones are called **Adams-Bashforth methods**.

### Adams-Bashford Methods

For this subfamily, $b_{-1}= 0$. Depending on the value of $p$, we obtain different methods:

* With $p = 0$, we obtain the **Explicit Euler Method**:

$$
u_{n + 1} = u_n +hf(t_n, u_n)
$$

* With $p = 1$, we obtain the **Two-Step Adams-Bashford Method**:

$$
u_{n + 1} = u_n + \frac{h}{2}(3f_n - f_{n - 1})
$$

* With $p = 2$, we obtain the **Three-Step Adams-Bashford Method**:

$$
u_{n + 1} = u_n + \frac{h}{12}(23f_n - 16f_{n - 1} + 5f_{n - 2})
$$

* With $p = 3$, we obtain the **Four-Step Adams-Bashford Method**:

$$
u_{n + 1} = u_n + \frac{h}{24}(55f_n - 59f_{n - 1} + 37f_{n - 2}-9f_{n - 3})
$$


### Adams-Moulton Methods

For this subfamily, $b_{-1}\neq 0$. Depending on the value of $p$, we obtain different methods:

* With $p = -1$, we obtain the **Implicit Euler Method**:

$$
u_{n + 1} = u_n +hf(t_{n+1}, u_{n+1})
$$

* With $p = 0$, we obtain the **Cranck-Nicolson Method**:

$$
u_{n + 1} = u_n + \frac{h}{2}(f_n + f_{n + 1})
$$

* With $p = 1$, we obtain the **Two-Step Adams-Moulton Method**:

$$
u_{n + 1} = u_n + \frac{h}{12}(5f_{n+1} + 8f_n - f_{n - 1})
$$

* With $p = 2$, we obtain the **Three-Step Adams-Moulton Method**:

$$
u_{n + 1} = u_n + \frac{h}{24}(9f_{n+1} + 19f_n - 5f_{n - 1}+f_{n - 2})
$$

* With $p = 2$, we obtain the **Three-Step Adams-Moulton Method**:

$$
u_{n + 1} = u_n + \frac{h}{720}(251f_{n+1} + 646f_n - 264f_{n - 1}+106f_{n - 2}-19f_{n - 3})
$$

The $q$-step *Adams-Moulton* methods have order $q + 1$.

## BFD Methods

The **Backwards Differentiation Formulas** (BDFs) are *implicit multi-step methods* derived using the interpolation polynomial for $y$ for the points $t_{n+1}$, $t_n$ , ..., $t_{n-k}$, differentiating it and asserting that it is equal to $f_{n+1}$.

* For $k = 0$, we obtain:

    $$
    y(t) =  y_{n + 1} + (t - t_{n + 1})\frac{y_{n + 1}- y_n}{t_{n + 1}- t_n}
    $$

    and

    $$
    y'(t) \approx  \frac{y_{n + 1}- y_n}{t_{n + 1}- t_n} = f(t_{n + 1}, y_{n + 1})
    $$

    so, we obtain the **Implicit Euler Method**:

    $$
    u_{n + 1} = u_n +hf(t_{n+1}, u_{n+1})
    $$

* For $k = 1$, we obtain:

    $$
    u_{n + 1} = \frac43u_n - \frac13u_{n - 1} + \frac23 h f_{n + 1}
    $$

Only the first 6 BDF formulas are zero-stable, so they are the only ones used. The first two are *A-stable*, which means that they are a very good option for stiff problems.

## Properties of Multi-Step Methods

> **First Dahlquist Barrier**: There isn’t any zero-stable, $p$-step linear multistep method with order greater than $p + 1$ if $p$ is odd, $p + 2$ if $p$ is even.

> **Second Dahlquist barrier**: A linear explicit multistep method can be neither *A-stable*, nor $\theta$-stable. Moreover, there is no *A-stable* linear multistep method with order greater than $2$. Finally, for any
$\theta \in (0, \pi/2)$, there only exist $\theta$-stable $p$-step linear multistep methods of order $p$ for $p = 3$ and $p = 4$. The theorem states that the maximum order of a multi-step method that is *A-stable* is $2$, and that among all order $2$ multi-step methods, the one with the smallest local error is the Crank-Nicolson.

Why not using always Crank-Nicolson?
* Sometimes order 2 is not enough.
* It has the draw back that it doesn't dump errors for some (stiff) problems.

