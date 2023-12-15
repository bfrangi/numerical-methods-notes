# Finite Difference Methods: Differentiation, Interpolation and Integration (I)

## Differentiation: Finite Differences

The definition of derivative is:

>$$\left.\frac{df}{dx}\right|_{x = x_0} = f'(x_0) = \lim_{h \rightarrow 0} \frac{f(x_0+h)-f(x_0)}{h}$$

Therefore, it is reasonable to think that:

> For small enough $h$:
> $$f'(x_0) \approx \frac{f(x_0+h)-f(x_0)}{h}$$

In this way, we approximate the differential operator $D = \frac{d}{dx}$ by an expression involving the values of the function at a discrete set of points. This kind of formula is called a **Finite Difference Formula**. The discrete set of points is called the **Stencil** of the formula.

All *finite difference formulas* are based on the concept of the *Taylor Series*:

>A function $f(x)$ continuous and differentiable in the interval $[x_0, x]$ can be approximated by:
>$$f(x) = \sum_{k = 0}^{\infty} \frac{f^{(k)}(x_0)}{k!}(x - x_0)^k = \sum_{k = 0}^{n} \frac{f^{(k)}(x_0)}{k!}(x - x_0)^k + R_n$$
>Where $R_n = \frac{1}{(n+1)!}(x - x_0)^{n+1}f^{(n+1)}(\xi)$, usually abbreviated as $O(h^{n+1})$.

The exponent of $O(h^{n+1})$ is the order of the finite difference formula. It means that the error made when approximating the derivative with a certain formula with $O(h^p)$ will have an error of the order of $h^p$. As we normally use $h<<1$, the error will be smaller for larger order formulas.

The choice of stencil is very important. Normally we use equispaced stencils, which means that all consecutive points are spaced by a constant amount $(x_{i+n} = x_i + nh)$. A general stencil would be:

| $x$       | $x_{i-n}$ | $\cdots$ | $x_{i-1}$ | $x_{i}$ | $x_{i+1}$ | $\cdots$ |$x_{i+m}$ |
| ---------:| :-------: | :------: | :-------: | :-----: | :-------: | :------: | :------: |
| $y =f(x)$ | $y_{i-n}$ | $\cdots$ | $y_{i-1}$ | $y_{i}$ | $y_{i+1}$ | $\cdots$ |$y_{i+m}$ |

In order to derive the dinite difference formulas, it is useful to define a set of linear differential operators:

>**Difference Operators:**\
>Forward Difference: $$\Delta u = u_{i+1}-u_i$$
>Backward Difference: $$\nabla u = u_i-u_{i-1}$$
>Central Difference: $$\delta u = u_{i+1/2}-u_{i-1/2}$$

### Finite Forward Difference Formulas

Forward difference formulas are those that enable us to compute the value $y_i'$ of the derivative of a function $y$ at a certain point $x_i$ based on the value $y_i$ of the function at that same point $x_i$ and a certain number of points larger than $x_i$ ($y_{i+1}, y_{i+2}, ...$).

These formulas derive from the taylor expansion of $y$ around $x_i$ up to order $N$:

$$y_{i+n} = \sum_{k = 0}^N \frac{(x_{i+n}-x_i)^k}{k!}y_i^{(k)} + O(h^{N+1})$$

The order up to which the expansion is used, $N$, is the order of the finite difference formula. If we consider $x_{i+1} - x_i = h$, then $x_{i+n} - x_i = nh$, and the formula is simplified to:

$$y_{i+n} = \sum_{k = 0}^N \frac{(nh)^k}{k!}y_i^{(k)} + O(h^{N+1})$$

If we want to find the expression of the derivative $y_i^{(r)}$, then we need to find the previous expansion up to $N \geq r$ and evaluate it at $n = 0, 1, ..., N$. This will allow us to build a system of equations with enough equations to find the coefficients $a_n$ of:

$$y_i^{(r)} = a_{0}y_{i} + a_{1}y_{i+1} + \cdots+ a_Ny_{i+N}+ O(h^N)$$

The order of forward finite difference formulas is usually equal to $k - 1$, where $k$ is the number of points of the function $y$ on which the formula depends.

<details>
<summary>Example</summary>

**Example: Second Order Finite Forward Difference Formula for the First Derivative**

We want to find the second order finite difference formula, so $N = 2$. We know this is possible because, as we want the first derivative ($r = 1$), we can verify that $N = 2 > 1 = r$. 

As $N =2$, we need to expand $y_{i+n}$ up to the second order and evaluate for $n = 0, 1, 2:$

$$y_{i+n} = \sum_{k = 0}^2 \frac{(nh)^k}{k!}y_i^{(k)} + O(h^{3}) = y_i + nhy_i' + \frac12 (nh)^2y_i'' + O(h^3)$$

Then:


$$y_i = y_i$$

$$y_{i+1} = y_i + hy_i'+\frac12h^2y_i''+ O(h^3)$$

$$y_{i+2} = y_i + 2hy_i'+2h^2y_i''+ O(h^3)$$

Now we have a system of equations to solve:

$$y_i' = a_{0}y_{i} + a_{1}y_{i+1} + a_{2}y_{i+2} + O(h^2)$$

Substituting:

$$y_i' = a_{0}y_{i} + a_{1}\left(y_i + hy_i'+\frac12h^2y_i''+ O(h^3)\right) +\\ a_{2}\left( y_i + 2hy_i'+2h^2y_i''+ O(h^3) \right) + O(h^2)=\\
=( a_0+a_1+a_2 )y_i + h(a_1+2a_2)y_i' + \frac{h^2}{2}(a_1 + 4a_2)y_i''$$

Finally, $y_i' =( a_0+a_1+a_2 )y_i + h(a_1+2a_2)y_i' + \frac{h^2}{2}(a_1 + 4a_2)y_i''$ requires:

$$a_0+a_1+a_2 = 0$$
$$h(a_1+2a_2) = 1$$
$$\frac{h^2}{2}(a_1+4a_2) = 0$$

yielding $a_0=-\frac{3}{2h}$, $a_1 = \frac{2}{h}$ and $a_2 = -\frac{1}{2h}$:

>The *Second Order Finite Forward Difference Formula for the First Derivative* is:
>$$y_i' = \frac{-3y_i+4y_{i+1}-y_{i+2}}{2h} + O(h^2)$$
</details>

### Finite Backward Difference Formulas

Backward difference formulas are those that enable us to compute the value $y_i'$ of the derivative of a function $y$ at a certain point $x_i$ based on the value $y_i$ of the function at that same point $x_i$ and a certain number of points smaller than $x_i$ ($y_{i-1}, y_{i-2}, ...$).

Again, these formulas derive from the taylor expansion of $y$ around $x_i$ up to order $N$ that we simplified before to:

$$y_{i+n} = \sum_{k = 0}^N \frac{(nh)^k}{k!}y_i^{(k)} + O(h^{N+1})$$

If we want to find the expression of the derivative $y_i^{(r)}$, then we need to find the previous expansion up to $N \geq r$ and evaluate it at $n = 0, -1, ..., -N$. This will allow us to build a system of equations with enough equations to find the coefficients $a_n$ of:

$$y_i^{(r)} = a_{0}y_{i} + a_{1}y_{i-1} + \cdots+ a_Ny_{i-N}+ O(h^N)$$

The order of backward finite difference formulas is usually equal to $k - 1$, where $k$ is the number of points of the function $y$ on which the formula depends. 

<details>
<summary>Example</summary>

**Example: Second Order Finite Backward Difference Formula for the First Derivative**

We want to find the second order finite difference formula, so $N = 2$. We know this is possible because, as we want the first derivative ($r = 1$), we can verify that $N = 2 > 1 = r$. 

As $N =2$, we need to expand $y_{i+n}$ up to the second order and evaluate for $n = 0, 1, 2:$

$$y_{i+n} = \sum_{k = 0}^2 \frac{(nh)^k}{k!}y_i^{(k)} + O(h^{3}) = y_i + nhy_i' + \frac12 (nh)^2y_i'' + O(h^3)$$

Then:


$$y_i = y_i$$

$$y_{i-1} = y_i - hy_i'+\frac12h^2y_i''+ O(h^3)$$

$$y_{i-2} = y_i - 2hy_i'+2h^2y_i''+ O(h^3)$$

Now we have a system of equations to solve:

$$y_i' = a_{0}y_{i} + a_{1}y_{i+1} + a_{2}y_{i+2} + O(h^2)$$

Substituting:

$$y_i' = a_{0}y_{i} + a_{1}\left(y_i - hy_i'+\frac12h^2y_i''+ O(h^3)\right) +\\ a_{2}\left( y_i - 2hy_i'+2h^2y_i''+ O(h^3) \right) + O(h^2)=\\
=( a_0+a_1+a_2 )y_i - h(a_1+2a_2)y_i' + \frac{h^2}{2}(a_1 + 4a_2)y_i''$$

Finally, $y_i' =( a_0+a_1+a_2 )y_i - h(a_1+2a_2)y_i' + \frac{h^2}{2}(a_1 + 4a_2)y_i''$ requires:

$$a_0+a_1+a_2 = 0$$
$$h(a_1+2a_2) = -1$$
$$\frac{h^2}{2}(a_1+4a_2) = 0$$

yielding $a_0=\frac{3}{2h}$, $a_1 = -\frac{2}{h}$ and $a_2 = \frac{1}{2h}$:

>The *Second Order Finite Backward Difference Formula for the First Derivative* is:
>$$y_i' = \frac{3y_i-4y_{i-1}+y_{i-2}}{2h}+O(h^2)$$

Notice that this is just the formula for the forward difference but substituting $h$ by $-h$ and substituting $i$ by $i- 2$.
</details>

### Finite Centered Difference Formulas

Centered difference formulas are those that enable us to compute the value $y_i'$ of the derivative of a function $y$ at a certain point $x_i$ based on the value $y_i$ of the function at that same point $x_i$ and a certain number of points of which half are smaller than $x_i$ ($y_{i-1}, y_{i-2}, ...$) and the other half are larger than $x_i$ ($y_{i+1}, y_{i+2}, ...$).

Again, these formulas derive from the taylor expansion of $y$ around $x_i$ up to order $N$ that we simplified before to:

$$y_{i+n} = \sum_{k = 0}^N \frac{(nh)^k}{k!}y_i^{(k)} + O(h^{N+1})$$

If we want to find the expression of the derivative $y_i^{(r)}$, then we need to find the previous expansion up to even $N \geq r$ and evaluate it at $n = - N/2, ..., -1 0, 1, ..., N/2$. This will allow us to build a system of equations with enough equations to find the coefficients $a_n$ of:

$$y_i^{(r)} = a_{-N/2}y_{i-N/2} +...+ a_{-1}y_{i-1} +\\ + a_{0}y_{i} + a_{1}y_{i+1} + \cdots+ a_{N/2}y_{i+N/2}+ O(h^{N+1})$$

The order of centered finite difference formulas is usually equal to the number of points of the function $y$ on which the formula depends.

<details>
<summary>Example</summary>

**Example: Second Order Finite Centered Difference Formula for the First Derivative**

We want to find the second order finite difference formula, so $N = 2$. We know this is possible because, as we want the first derivative ($r = 1$), we can verify that $N = 2 > 1 = r$. 

As $N =2$, we need to expand $y_{i+n}$ up to the second order and evaluate for $n = 0, 1, 2:$

$$y_{i+n} = \sum_{k = 0}^2 \frac{(nh)^k}{k!}y_i^{(k)} + O(h^{3}) = y_i + nhy_i' + \frac12 (nh)^2y_i'' + O(h^3)$$

Then:


$$y_i = y_i$$

$$y_{i+1} = y_i + hy_i'+\frac12h^2y_i''+ O(h^3)$$

$$y_{i-1} = y_i - hy_i'+\frac12h^2y_i''+ O(h^3)$$

Now we have a system of equations to solve:

$$y_i' = a_{-1}y_{i-1} + a_{0}y_{i} + a_{1}y_{i+1} + O(h^3)$$

Substituting:

$$y_i' = a_{-1}\left(y_i - hy_i'+\frac12h^2y_i''+ O(h^3)\right) + a_{0}y_{i}+\\ +\ a_{1}\left(y_i + hy_i'+\frac12h^2y_i''+ O(h^3)\right) + O(h^3)=$$
$$=(a_{-1}+a_0+a_1)y_i + h(-a_{-1} + a_0 + a_1)y_i' + \frac{h^2}{2} (a_{-1} + a_1)y_i''+O(h^3)$$

Finally, $y_i' =(a_{-1}+a_0+a_1)y_i + h(-a_{-1} + a_0 + a_1)y_i' + \frac{h^2}{2} (a_{-1} + a_1)y_i''$ requires:

$$a_{-1}+a_0+a_1 = 0$$
$$h(-a_{-1} + a_0 + a_1) = 1$$
$$\frac{h^2}{2} (a_{-1} + a_1) = 0$$

yielding $a_{-1} = -\frac{1}{2h}$, $a_0=0$ and $a_1 = \frac{1}{2h}$:

>The *Second Order Finite Centered Difference Formula for the First Derivative* is:
>$$y_i' = \frac{y_{i+1}-y_{i-1}}{2h} + O(h^3)$$
</details>

Online finite differences calculator [here](https://web.media.mit.edu/~crtaylor/calculator.html).

