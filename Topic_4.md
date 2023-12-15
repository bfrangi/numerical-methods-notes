# Numerical Integration

Numerical integration formulas, also called **quadrature formulas**, are formulas that help us to approximate the value of integrals at certain points. They are useful because some functions cannot be integrated analytically.

Most quadrature formulas use interpolation to substitute functions that are difficult to integrate by simple polynomials.

> In order to find the approximation of the integral of a function, we can substitute the function by an interpolation polynomial of order $n$:
> $$I (f) = \int_{x_0}^{x_n}f(x)dx \approx \int_{x_0}^{x_n}f_n(x)dx = I_n(f)$$
> Then, if $f\in C^0([x_0,x_n])$ then the quadrature error $E_n (f) = I(f) -I_n(f)$ satisfies:
> $$|E_n(f)| \leq \int_{x_0}^{x_n}|f(x)-f_n(x)|dx\leq\\ \leq (x_n-x_0)||f(x)-f_n(x)||_\infty$$

From this, we can learn two things:
1. This error estimation is very conservative, as it considers the worst case scenarios.
2. Using high order polynomial interpolation is not a good idea to calculate an integral, as it generates regions (especially at the interval ends) where $||f(x)-f_n(x)||$ is very large.

## Lagrange Quadrature

A popular interpolatory quadrature formula uses **Lagrange's Polynomial Interpolation**:

> Using *Lagrange's Polynomial Interpolation*, we can find the following **Quadrature Formula**:
> $$I_n(f) = \int_{x_0}^{x_n}P_n(x)dx = \int_{x_0}^{x_n}\sum_{k = 0}^nf(x_k)l_k(x)dx = \\ =\sum_{k = 0}^nf(x_k)\int_{x_0}^{x_n}l_k(x)dx = \sum_{k = 0}^n\omega_kf(x_k)$$
> where $x_k$ are the nodes of the quadrature and $\omega_k$ are their weights, defined by *Lagrange's Quadrature Formula*:
> $$\omega_k = \int_{x_0}^{x_n}l_k(x)dx$$

This kind of quadrature is called the *Newton-Cotes Quadrature* if the nodes are equispaced.

## Properties of Quadrature Formulae

Quadrature formulae can be *open* or *closed*. *Open* quadrature formulae do not include the ends of the integration interval in their formula, whereas *closed* quadrature fromulae do include them.

We can evaluate interpolatory quadrature formulae using the degree of exactness:

> The **degree of exactness** of a quadrature formula is the maximum integer $r\geq 0$ for which $I_n(f) = I(f)$, $\forall f \in \mathbb{P}_r$. Any interpolatory quadrature formula that makes use of $n + 1$ distinct nodes has degree of exactness equal to at least $n$.

In other words, the degree of exactness of an interpolatory quadrature formula measures what is the maximum degree $r$ of polynomials whos integral is approximated exactly by the quadrature formula.

Another way to qualify an interpolatory quadrature formula is by its infinitessimal order:

> The Infinitessimal Order of an interpolatory quadrature formula is the maximum integer $p$ such that $|I(f) - I_n(f)| = O(h^p)$.

If we define the following quantities:

> For $n$ even:
> $$M_n = \begin{cases}
\int_0^n t^2 (t-1)(t-2)\cdots(t-n)dt\ \text{for closed formulas},\\
\int_{-1}^{n+1} t^2 (t-1)(t-2)\cdots(t-n)dt\ \text{for open formulas}.        
\end{cases}$$
> For $n$ odd:
> $$K_n = \begin{cases}
\int_0^n t (t-1)(t-2)\cdots(t-n)dt\ \text{for closed formulas},\\
\int_{-1}^{n+1} t (t-1)(t-2)\cdots(t-n)dt\ \text{for open formulas}.       
\end{cases}$$


Then, we can find the order of the error:

> Suppose the following quadrature formula for $n+1$ points:
> $$\int_{x_0}^{x_n}f(x)dx \approx \sum _{i = 0}^n \omega_if(x_i)$$	
> Then, the error involved is:\
> For $n$ even:
> $$\text{Error}= \frac{h^{n+3}f^{(n+2)}(\xi)}{(n+2)!}M_n = O(h^{n+3})$$
> For $n$ odd:
> $$\text{Error}= \frac{h^{n+2}f^{(n+1)}(\xi)}{(n+1)!}K_n = O(h^{n+2})$$

## Interpolatory Quadrature Rules

### Midpoint Rule

> **Midpoint Rule**\
> The *Midpoint Rule* is defined as
> $$\int_a^bf(x)dx \approx (b-a) f\left(\frac{a+b}{2}\right)$$
> It has degree of exactness $1$ and is $O(h^2)$.

<!-- <details>
<summary>Proof (finish later)</summary>

From the previous topic on interpolation, we know that the error of the *Lagrange Interpolation* is:

> $$E_n(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!}\prod_{i = 0}^n(x-x_i)\ x_0 < \xi < x_n$$

We can find the error of the quadrature by integrating this expression. As $n = 1$ (we use two points to do this interpolation):

$$E_1(f) = \int_a^b\left(\frac{f^{(1+1)}(\xi)}{(1+1)!}\prod_{i = 0}^1(x-x_i)\right)dx =\\= \frac{f''(\xi)}{2}\int_a^b(x-a)(x-b)dx$$


If $f\in C^2([a, b])$, its error is:
$$E_0 (f) = \frac{h^3}{24}f''(\xi)$$

</details> -->


### Trapezoid Rule

> **Trapezoid Rule**\
> The *Trapezoid Rule* is obtained by substituting the function by its *Lagrange interpolaion of degree 1*.
> $$\int_a^bf(x)dx \approx \frac{h}{2}\left(f(a)+f(b)\right) = I_1(f)$$
> where $h = b-a$. It has degree of exactness $1$ and is $O(h^2)$.

### Simpson's $1/3$ Rule

> **Simpson's $1/3$ Rule**\
> The *Simpson's $1/3$ Rule* is obtained by substituting the function by its *Lagrange interpolaion of degree 2*.
> $$\int_a^bf(x)dx \approx \frac{h}{6}\left(f(a)+4f\left(\frac{a+b}{2}\right) + f(b)\right) = I_2(f)$$
> where $h = b-a$. It has degree of exactness $3$ and is $O(h^5)$.

## Composite Interpolatory Quadratures

As we have seen, low order interpolations do not give all the precision that we would want; and high order interpolations are not a good idea. The other option we have is to divide the intergration interval into subintervals in which we can apply low order interpolations and obtain better results.

### Composite Midpoint Rule

> **Composite Midpoint Rule**\
> If we divide the integration interval into $m$ subintervals (in composite midpoint, the total number of points is $m + 2$, one in the middle of each interval plus the two end points), we obtain the following expression for the *Composite Midpoint Rule*:
> $$I_{0, m}(f) = \frac{b-a}{m}\sum_{k =0}^{m-1}f(x_k)$$
> It has degree of exactness $1$ and is $O(h^2)$. Its error is:
> $$E_{0, m}(f) = \frac{b-a}{24}H^2f''(\xi) $$
> where $H = \frac{b-a}{m}$ and $\xi$ is some value $\xi \in(a, b)$.

### Composite Trapezoid Rule

> **Composite Trapezoid Rule**\
> If we divide the integration interval into $m$ subintervals, we obtain the following expression for the *Composite Trapezoid Rule*:
> $$I_{1, m}(f) = \frac{b-a}{2m}\sum_{k =0}^{m-1}\left(f(x_k)+f(x_{k+1})\right)$$
> It has degree of exactness $1$ and is $O(h^2)$. Its error is:
> $$E_{1, m}(f) = -\frac{b-a}{12}H^2f''(\xi) $$
> where $H = \frac{b-a}{m}$ and $\xi$ is some value $\xi \in(a, b)$.


### Composite Simpson's $1/3$ Rule

> **Composite Simpson's $1/3$ Rule**\
> If we divide the integration interval into an *even* number $m$ of subintervals (we have $m+1$ points, where $x_0 = a$ and $x_m = b$), we obtain the following expression for the *Composite Simpson's $1/3$ Rule*:
> $$I_{2, m}(f) = \frac{b-a}{3m}\left[
    f(a)+2\sum_{k = 2, 4, ...}^{m - 2}f(x_k) + 4\sum_{k =1,3,...}^{m-1}f(x_k) +f(b)
    \right]$$
> It has degree of exactness $3$ and is $O(h^4)$. Its error is:
> $$E_{0, m}(f) = -\frac{b-a}{180}H^4f^{(4)}(\xi) $$
> where $H = \frac{b-a}{m}$ and $\xi$ is some value $\xi \in(a, b)$.

