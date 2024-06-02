---
title: Fourier Transforming Polar Domain - An Anatomy to Fourier Bessel Transformation
date: 2024-6-1 21:11:05 +/-TTTT
categories: [Mathematics, Analysis]
tags: [Fourier, Bessel, Analysis]     # TAG names should always be lowercase
math: true
---

![/pics/BesselBasis/BesselBasis.png](/pics/BesselBasis/BesselBasis.png)

## Preface

I started on this topic around 4 days ago, since some of my research algorithms requires Fourier transforming some signals defined over a 2D polar system. I thought this is a simple problem with some change of variables, but it turns out to be a massive rabbit hole with many concepts beyond my knowledge base. 

Therefore, I documented my entire research process and explained in as much detail as necessary to understand the Fourier Bessel Transformation. Trust me, it is nontrivial but it is not that hard.
Fourier Bessel Transformation is really closely related to the buzz word of "Spherical Harmonics". It essentially defined a "Disk Harmonics" following almost the same principles.

To best understanding this page, you should be farmiliar with at minimum:

- The continuous and discrete Fourier transformation
- The idea of functions being decomposed into a combination of basis functions (The idea of Harmonic Functions)

And, here we go.
## 1: Recap on Fourier transforms

#### 1.1 Continuous Fourier Transform

The complex exponential forms an orthogonal basis for the space of absolutely integrable functions:

$$e^{i\omega x} = \cos(\omega x) + i \sin(\omega x)$$

We can show the orthognality property by:

$$\langle  e^{i \omega_{1} x}, e^{i \omega_{2} x} \rangle_{x} = \int_{\mathbb{R} } e^{i \omega_{1}x} [e^{i \omega_{2}x}] \, dx  = 2\pi\delta (\omega_{1}-\omega_{2})$$

Which means:

$$\omega_{1}\neq \omega_{2}\implies \langle e^{i\omega_{1}x} , e^{i \omega_{2}x} \rangle_{x}=0 $$

To form an **orthnormal basis**, we can simply add a normalization factor of:

$$\frac{1}{\sqrt{ 2\pi } } e^{i\omega x}$$

which has the norm of 1 over the domain.

The fact that complex exponential basis covers all integrable functions came directly after the Fourier transform. intuitively, every integrable function is an an "weighted sum":

$$f(x) = \lim_{ n \to \infty } \sum_{-n} ^n W(\omega_{n}) e^{i \omega_{n} x} = \lim_{ n \to \infty } \sum_{-n} ^n |W(\omega_{n})|e^{i\varphi(\omega_{n})} e^{i \omega_{n} x}$$

Note: When we use countably infinitely many basis functions, our weights forms the **Fourier Series**. It **requires our function to be periodic**, and can only model the function within some period interval of $T$.

The weight $W$ is complex, it controls the magnitude (or amplitude) and the phase (the start point) for each basis function.

- These complex exponential essentially are "phase shiftable cosine waves". You can think the $i$ component as "below the ground", anything on this dimension will not contribute to the function.

- When the function is completely composed of frequencies without shifting, it is essentially an weighted sum of cos waves. The phase shift $e^{i\varphi(\omega_{n})}$ adjusts the amount of "sine" in this composition of "cosines".

When we want to sum over every possible frequency, we have an integral that can represent the function over the entire possible domain:

$$f(x) =\int _{-\infty} ^{\infty} F(\omega) e^{i \omega x } d\omega$$

- $F$ assigns weight for every $\omega$

The forward Fourier transform is an **"inner product"** that gives the "alignment" of a function against a particular angular frequency basis:

$$F(\omega)= \langle f,  e^{i\omega\_{} } \rangle_{x}$$

The 1 dimensional Fourier transformation (and its inverse transformation) is given by the following formula:

$$\begin{align}
(\mathcal Ff)(x) &=F(\omega) = \int _{-\infty} ^\infty f(x) e^{- i \omega x} \, dx \\
\mathcal (\mathcal F^{-1} F)(\omega) &= f(x) = \frac{1}{2\pi}\int _{-\infty} ^{\infty} F(\omega) e^{i \omega x } \, d\omega
\end{align} \quad (1.0) $$

- $\omega$ is the frequency parameter
- $F(\omega)$ tells you how much of the function is composed of this frequency. Its the $W$ we are looking for.
- $\mathcal{F}$ is the fourier function operator.
- $-$ the negative sign is just applied by convention. All Fourier transforms of real valued functions are symmetric w.r.t origin in the frequency domain, known as **Hermitian symmetry**.
- $\frac{1}{2\pi}$ is the normalization factor. It came from the dirac delta function 

$$\int _{-\infty} ^\infty e^{i\omega(x-x_{0})} \, d\omega = 2\pi \delta(x-x_{0})$$

during the proof of the inverse transform

Fourier transform is a function operator, in which it takes in a function that is **Absolutely Lebesgue Integrable** (which means it has an finite Lebesgue integral (or intuitively speaking have a finite "area under curve" over its domain)):

$$\int_{D} \lvert f(x) \rvert  \, dx < \infty$$

- where $D \subset (-\infty, \infty)$ is the function domain

Fourier transforms has the following properties:

- (linearity) The **addition / subtraction** of 2 functions in time domain is the same as the addition and subtraction of 2 functions in frequency domain:

$$f(x) \pm g(x) \iff F(\omega) \pm G(\omega)\quad(1.1.0)$$

- (linearity) The **multiplication** by a constant in the time domain is the multiplication by the constant in the frequency domain:

$$kf(x) \iff kF(\omega)\quad(1.1.1)$$

- The **convolution** of 2 function in time domain is their point-wise multiplication in the frequency domain.

$$(f * g)(x) \iff F(\omega)G(\omega)\quad(1.1.2)$$

- The **point-wise multiplication** of 2 functions in time domain is their weighted convolution in the frequency space (weighhted by the normalizing factor):

$$f(x)\cdot g(x) \iff  \frac{(F* G)(\omega)}{2\pi} \quad(1.1.4)$$

- The **Ideal Sampling / Shannon-Nyquist Theorem:** to idealy reconstruct any function, the sampling interval must be 1 over 2 times the period of its maximum non-zero frequency component. (Derived from Shah function property)

- Here is a list of commonly used Fourier pairs:

| Name     | Time Domain Function                                 | Frequency Domain Function                                         |
| -------- | ---------------------------------------------------- | ----------------------------------------------------------------- |
| Box      | $f(x)=1:\|x\| \leq \frac{1}{2} \text{ otherwise } 0$ | $sinc(\omega) : \frac{\sin(\pi \omega)}{\pi \omega}$              |
| Gaussian | $f(x)=e^{-\pi x^2}$                                  | $f(\omega)=e^{-\pi \omega^2}$                                     |
| Constant | $f(x)=1$                                             | $f(\omega)=\delta (\omega)$                                       |
| Sinusoid | $f(x)=\cos x$                                        | $f(\omega)=\pi (\delta (1-2 \pi \omega) +\delta (1+2\pi \omega))$ |
| Shah     | $f(x)=III_{T}(x)$                                    | $f(\omega)=III_{\frac{1}{T} }(\omega)$                             |

#### 1.2 Discrete Fourier Transform:

The discrete Fourier transform operates similarly as an inner product. For each angular frequency component, we take the dot product:

$$F(\omega)= \langle f,  e^{i\omega\_{} } \rangle$$

Since we do not have the analytic representation of $F$, we take $k$ samples along its domain, and do the sum:

$$F(\omega_{n}) = \sum_{i=1}^k f(x_{k}) \cdot e^{i\omega x_{k} }= \begin{bmatrix}
f(x_{1}) \\
\dots \\
f(x_{k})
\end{bmatrix}^T \cdot \begin{bmatrix}
 e^{i\omega x_{1} } \\
\dots \\
e^{i\omega x_{k} }
\end{bmatrix}$$

When we are looking for process multiple frequencies, we can represent the transform as a matrix vector product:

$$F= \begin{bmatrix}
\omega_{1}(x_{1}) & \dots & \omega_{1}(x_{k}) \\
&\dots \\
\omega_{n}(x_{1}) & \dots & \omega_{n}(x_{k})
\end{bmatrix} \begin{bmatrix}
f(x_{1}) \\
\dots \\
f(x_{k})
\end{bmatrix}$$

A really interesting fact is, the matrix on the left is orthogonal, which mean, if $n=k$, we have:

$$\begin{bmatrix}
 f(x_{1}) \\
\dots \\
f(x_{k})
\end{bmatrix}= C \begin{bmatrix}
\omega_{1}(x_{1}) & \dots & \omega_{1}(x_{k}) \\
&\dots \\
\omega_{n}(x_{1}) & \dots & \omega_{n}(x_{k})
\end{bmatrix} ^T F$$

where $C=\frac{1}{2\pi}$ is the normalizing factor, as we have seen before in continuous Fourier transform.

For this algorithm to properly function, we need $k$ to be double the maximum frequency $\omega_{n}$ we are looking for. So the matrix will be $2n \times 2n$, and we extend frequencies to $\omega_{2n}$, but only below $\omega_{n}$ will the transform be accurate.

## 2: Angular and Radial Basis Functions

#### 2.1 Cartesian case

In the Cartesian 2D Fourier transform, the basis functions are the same along each dimension, as:

$$e^{i \mathbf{k} \cdot \mathbf{r} } = e^{ik_{x} x} e^{ik_{y}y}$$

Where $\mathbf{r}=[x,y]$ is the position on $\mathbb{R}^2$, and $\mathbf{k}=[k_{x},k_{y}]$ is the "wave vector" at this point, representing wave strength on both axis.

To check orthogonal condition, we can do the inner product:

$$\langle e^{i \mathbf{k}_{1} \mathbf{r} }, e^{i\mathbf{k}_{2}\mathbf{r} } \rangle_{\mathbf{ r} }  = \int _{\mathbb{R} } e^{ik_{1,x}x} e^{ik_{2}x} \left[ \int _{\mathbb{R} } e^{ik_{1}y}e^{ik_{2}y} \, dx  \right]\, dy $$

Simplify:

$$= 4\pi^2 \delta (k_{1,x}-k_{2,x})\delta (k_{1,y}-k_{2,y})$$

Thus, orthogonality holds, and the 2D Fourier transform is thus given by:

$$\begin{align}
\mathcal{F} f(x,y) &= F(k_{x},k_{y}) = \int _{-\infty}^\infty  \int _{-\infty}^\infty f(x,y)   e^{ik_{x} x} e^{ik_{y}y}\, dx dy \\
\mathcal{F}^{-1} F(k_{x},k_{y}) &= f(x,y) =\frac{1}{(2\pi)^2} \int _{-\infty}^\infty  \int _{-\infty}^\infty F(k_{x},k_{y})   e^{ik_{x} x} e^{ik_{y}y}\, dk_{x} dk_{y}
\end{align}$$

Now, consider the change of basis for a function from Cartesian to Polar, as:

$$x = r \cos \varphi, y=r \sin \varphi$$

which implies:

$$\begin{align}
 dx &= \cos \varphi dr - r\sin \varphi d\varphi \\
dy &= \sin \varphi dr + r\cos \varphi d\varphi \\
dx\land dy &= r(\cos^2 \varphi + \sin^2 \varphi) dr \land d\varphi =r dr\land d\varphi
\end{align}$$

Therefore, we are looking for some new separable basis function:

$$\Psi(r,\varphi ; \mathbf{k}_{r,\varphi})=R(r;k_{r})\Phi(\varphi;k_{\varphi})$$

such that the Fourier premise works and each polar functions can be represented via:

$$f(r, \varphi) = \int _{\mathbb{R} } \int _{0}^{2\pi} W(\mathbf{k}) R(r;k_{r}) \Phi (\varphi;k_{\varphi}) r\, d\varphi  \, dr $$

- where $W$ is some complex weight assigned to basis at $n,m$

To actually solve for these basis functions, we can take advantage for the Laplacian representation of Fourier transfrom.
#### 2.2 Helmholtz Equation and Fourier Transform

For any function operator $\mathcal{F}$, we say $f$ is an **Eigenfunction** if and only if:

$$\exists \lambda \in \mathbb{R}: \mathcal F f =-\lambda f, \mathcal F f +\lambda f =0$$

The **Laplacian Operator** $\Delta$ or $\nabla^2$, is the function operator defined as follows:

$$\nabla^2 \mathcal f = \sum_{i=1}^n \frac{ {\partial^2 f} }{\partial x_{n}^2}$$

The **Helmholtz Equation** is a special Eigenfunction problem in the context of Laplacian operator, where:

$$\nabla^2 f =-k^2 f $$

where $k^2$ is some constant.

This equation is mostly used in physics context (which for now I am not particularly good at yet). However, we bring it up for the following property:

**Theorem 2.2.1** Let $\mathcal{F}$ be the Fourier operator and $f$ be a 1D integrable function, we have:

$$\mathcal{F}(\nabla^2f)(\omega) = -\lvert \omega \rvert^2 (\mathcal{F}f)(\omega)$$

**Proof:**

We start with the laplacian of our basis function:

$$\nabla^2_{x} e^{i\omega x}= \frac{\partial^2}{\partial x^2} e^{i\omega x}= (i\omega)^2 e^{i\omega x} = -\omega^2 e^{i\omega x}$$

This looks interesting, so we can try swap it in

$$\begin{align}
\mathcal{F}(\nabla^2f)(\omega) &= \int _{\mathbb{R} } [\nabla^2 f(x)] e^{-i\omega x} \, dx \\ &= \int _{\mathbb{R} } \left[ \frac{\partial^2}{\partial x^2} f(x) \right] e^{-i\omega x} \, dx \\ &= \int _{\mathbb{R} } \frac{\partial}{\partial x} \left( \frac{\partial f(x)}{\partial x} \right) e^{-i\omega x} \, dx
\end{align}$$

Using **integration by parts**, we focus on the first derivative term:

$$\begin{align} &= \left. \frac{\partial f(x)}{\partial x} e^{-i\omega x} \right|_{-\infty}^{\infty} - \int _{\mathbb{R} } \frac{\partial f(x)}{\partial x} \left( \frac{\partial}{\partial x} e^{-i\omega x} \right) \, dx \\ &= 0 - \int _{\mathbb{R} } \frac{\partial f(x)}{\partial x} (-i\omega) e^{-i\omega x} \, dx \\ &= i\omega \int _{\mathbb{R} } \frac{\partial f(x)}{\partial x} e^{-i\omega x} \, dx \end{align}$$

The term 

$$\left. \frac{\partial f(x)}{\partial x} e^{-i\omega x} \right|_{-\infty}^{\infty}$$ 

goes to zero from the integrable function assumption, which requires function $f(x)$ and its derivatives go to 0 as $$x\to \pm \infty$$

Repeat:

$$\begin{align} i\omega \int _{\mathbb{R} } \frac{\partial f(x)}{\partial x} e^{-i\omega x} \, dx &= i\omega \left( \left. f(x) e^{-i\omega x} \right|_{-\infty}^{\infty} - \int _{\mathbb{R} } f(x) \left( \frac{\partial}{\partial x} e^{-i\omega x} \right) \, dx \right) \\ &= i\omega \left( 0 - (-i\omega) \int _{\mathbb{R} } f(x) e^{-i\omega x} \, dx \right) \\ &= (-1) \omega^2 \int _{\mathbb{R} } f(x) e^{-i\omega x} \, dx \\ &= -\omega^2 (\mathcal{F}f)(\omega) \end{align}$$

Q.E.D-Theorem 2.2.1

Note, this theorem also applies to wave vector $\mathbf{k}$ for arbitrary dimensional Cartesian frame, as:

$$\mathcal{F}(\nabla^2f)(\mathbf{k}) = -\lvert \mathbf{k} \rvert^2 (\mathcal{F}f)(\mathbf{k})$$

Now, remember that this Theorem works only if the Fourier basis has this property

$$\nabla^2_{x_{1},\dots,x_{n} } e^{i\mathbf{k}\cdot \mathbf{r} }= -\lvert \mathbf{k} \rvert ^2 e^{i \mathbf{k} \cdot \mathbf{r} }$$

Which means: **Fourier basis functions are Eigenfunctions for the Laplacian Operator of eigenvalue the size of the frequency vector, or solutions to the particular Helmholtz Equation.**

if we want our **"Polar Fourier Transform"** to work, we need the polar basis functions $\Psi (r,\varphi;\mathbf{k}_{r,\varphi}) =R(r;k_{r})\Phi(\varphi;k_{\varphi})$ to be such solutions:

$$\nabla^2_{x,y} \Psi(r,\varphi;\mathbf{k}_{r,\varphi}) = -\lvert \mathbf{k}_{x,y} \rvert^2 \Psi(r,\varphi;\mathbf{k}_{r,\varphi})$$

Using the change of variable formula, we notice:

$$\nabla^2_{x,y} = \nabla^{2}_{r} + \frac{1}{r^2} \nabla^2_{\varphi}$$

Where:

$$\nabla ^2_{r} = \frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{\partial}{\partial r} \right), \nabla_{\varphi}^2 = \frac{\partial^2}{\partial \varphi^2}$$

- Note: this can be directly derived using change of variables

This turns our objective equation into "Helmholtz Equation in Polar coordinates":

$$\Phi (\varphi)\nabla^2_{r}R(r)+ \frac{R(r)}{r^2} \nabla^2_{\varphi} \Phi(\varphi)  = - k_{r} ^2 R(r)\Phi(\varphi) $$

note: $|\mathbf{k}|$ is sames as $\lvert k_{r} \rvert$ under polar coordinates

Rewrite, we have:

$$\frac{ {\nabla_{r}^2 R(r)} }{R(r)} + \frac{1}{r^2} \frac{ {\nabla_{\varphi}^2 \Phi(\varphi)} }{\Phi(\varphi)} = -\lvert  k_{r} \rvert^2 $$

And we can use a trick called **Separation of Variables.**

Because $r^2$ is a function of $r$, we organize by:

$$\frac{ {r^2\cdot\nabla_{r}^2 R(r;k_{r})} }{R(r;k_{r})} + k_{r} ^2 =-\frac{ {\nabla_{\varphi}^2 \Phi(\varphi;k_{\varphi})} }{\Phi(\varphi;k_{\varphi})} =\lambda $$

- where $\lambda$ can be any constant for now

Which will give us 2 ordinary differential equations:

$$\begin{align} \\
 \nabla_{r}^2 R(r;k_{r}) &= \left( \frac{\lambda}{r^2}-k_{r}^2 \right) R(r;k_{r}) \\
\nabla_{\varphi}^2 \Phi(\varphi;k_{\varphi}) &= - \lambda \Phi(\varphi;k_{\varphi})
\end{align}$$

Wait, the second equation is just another Laplacian eigenvalue problem! let me rewrite:

$$\nabla_{\varphi}^2 \Phi(\varphi;k_{\varphi}) +\lambda \Phi(\varphi;k_{\varphi}) = 0$$

And here comes our original Fourier basis which exactly satisfies this property:

$$\Phi(\varphi; k_{\varphi}) = e^{i \sqrt{ \lambda } \varphi}$$

since we need $k_{\varphi}$ as a parameter, we can deduce:

$$\sqrt{ \lambda} = k_{\varphi} \implies \lambda =k_{\varphi}^2$$

Thus, we can separate this equation into 2 equations for each basis function:

$$\begin{align}
&\nabla^2_{\varphi} \Phi(\varphi) + k_{\varphi}^2 \Phi(\varphi) = 0 \\
&\nabla^2_{r}R(r) + \left( k_{r}^2 -\frac{k_{\varphi}^2}{r^2} \right)R(r)=0
\end{align}$$

Solving these equations separately, we will get our desider basis functions.
#### 2.3 Angular Basis Function

The angular basis function is given by:

$$\Phi_{k_{\varphi} }(\varphi) = \frac{1}{\sqrt{ 2\pi } }e^{ik_{\varphi} \varphi}$$

which is the same as the canonical Fourier basis.
## 3. Bessel Functions and Radial Basis

Our second equation looks like this:

$$\nabla^2_{r}R(r) + \left( k_{r}^2 -\frac{k_{\varphi}^2}{r^2} \right)R(r)=0 $$

And after some rewrite:

$$\frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{\partial R}{\partial r} \right)  + \left( k_{r}^2 - \frac{k_{\varphi}^2}{r^2} \right) R = 0$$
- Note: this $\frac{1}{r}$ and $r$ split came from the Laplacian in cylindrical coordinates

Then, we can use the product rule:

$$\frac{1}{r}\left( r \frac{ {\partial^2 R} }{\partial r^2}  + \frac{ {\partial R} }{\partial r}\right)  + \left( k_{r}^2 - \frac{k_{\varphi}^2}{r^2} \right) R = 0$$

multiply both side by $r^2$, we haveL

$$r^2 \frac{ {\partial^2 R} }{\partial r^2} + r \frac{ {\partial R} }{\partial r} + (r^2k_{r}^2 - k_{\varphi}^2)=0$$

Finally, we preform a change of variable:

$$\psi = k_{r}r, \frac{\partial}{\partial r}= k_{r}\frac{d}{d\psi}$$

Put this in, we have:

$$(r^2k_{r}^2) \frac{ {\partial^2R} }{\partial r^2} + (rk_{r}) \frac{ {\partial R} }{\partial r} + ((rk_{r})^2 - k_{\varphi}^2)=0$$

What we have is a Bessel's Equation.

The **Bessel's Differential Equation** is defined animals of the following shape:

$$x^2 \frac{ {\partial^2 f} }{\partial x^2} +x \frac{ {\partial f} }{\partial x} +(x^2-\alpha^2) f =0$$

The Bessel functions are know to have the following sets of solutions:

$$f= AJ_{\alpha}(r) + BY_{\alpha}(r)$$

Where:

$$J_{\alpha}(x) = \left( \frac{x}{2} \right)^{\alpha} \sum_{m=0}^\infty  \frac{ { {(-1)^m} } }{m!
\Gamma(m+\alpha+1)} \left( \frac{x}{2}\right) ^{2m }$$

Is called the **Bessel Function of the first kind.**

And

$$Y_{\alpha}(x)= \frac{ {J_{\alpha}(x) \cos(\alpha \pi)-J_{-\alpha}(x)} }{\sin (\alpha \pi)}$$

is called the **Bessel Function of the second kind** or **Neumann Functions**
$A$ and $B$ are just some constant quantities.
### 3.1 Deriving Bessel Functions

#### 3.1.1 Bessel Function of the First Kind

![/pics/BesselBasis/BesselFirstKind.png](/pics/BesselBasis/BesselFirstKind.png)

The **Frobenius Method** says that, for any differential equations of shape:

$$y'' +\frac{p(x)}{x}y' + \frac{q(x)}{x^2}y =0$$

we can try something like this:

$$y= x^\alpha \sum_{m=1}^\infty A_{m}x^m$$

Take the partial derivative we have:

$$\begin{align}
y'(x) &= \alpha x^{\alpha-1} \sum_{m=0}^\infty A_{m}x^m + x^\alpha \sum_{m=0}^\infty m A_{m}x^{m-1} \\
y''(x) &= \alpha(\alpha-1) x^{\alpha-2}  \sum_{m=0}^\infty A_{m}x^m  + 2\cdot \left(\alpha x^{\alpha-1} \sum_{m=0}^\infty m A_{m}x^{m-1} \right) \\
&+ x^{\alpha} \sum_{m=0}^\infty m(m-1)A_{m}x^{m-2}\end{align}$$

as a result for the multiplication rule.

Now to recover the equaiton, we multiply $y'$ by $x$ and $y''$ by $x^2$, we end up with:

$$\begin{align}
xy'(x) &= \alpha x^{\alpha} \sum_{m=0}^\infty A_{m}x^m + x^\alpha \sum_{m=0}^\infty m A_{m}x^{m} \\
x^2y''(x) &= \alpha(\alpha-1) x^{\alpha}  \sum_{m=0}^\infty A_{m}x^m  + 2\cdot \left(\alpha x^{\alpha} \sum_{m=0}^\infty m A_{m}x^{m} \right) \\
&+ x^{\alpha} \sum_{m=0}^\infty m(m-1)A_{m}x^{m} \\
x^2y(x)&= x^\alpha \sum_{m=2}^\infty A_{m-2}x^m\end{align}$$

The last term $x^2y$ is specially for Bessel Equation

And you start to see where this is going, don't you?

$$x^2 \frac{ {\partial^2 f} }{\partial x^2} +x \frac{ {\partial f} }{\partial x} +(x^2-\alpha^2) f =0$$

Doing some plugins, and we have our equation rewritten as:

$$\begin{align}
&\left[\alpha(\alpha-1) x^{\alpha}  \sum_{m=0}^\infty A_{m}x^m  + 2\cdot \left(\alpha x^{\alpha} \sum_{m=0}^\infty m A_{m}x^{m} \right) + x^{\alpha} \sum_{m=0}^\infty m(m-1)A_{m}x^{m}\right] \\
+& \left[ \alpha x^{\alpha} \sum_{m=0}^\infty A_{m}x^m + x^\alpha \sum_{m=0}^\infty m A_{m}x^{m} \right] + \left[ x^\alpha \sum_{m=2}^\infty A_{m-2}x^m \right] - \left[\alpha^2 x^\alpha \sum_{m=0}^\infty A_{m}x^m \right] = 0
\end{align}$$

Now merge summation terms:

$$x^\alpha \sum_{m=0}^\infty [((\alpha^2-\alpha) +2\alpha m+(m^2-m) +\alpha + m -\alpha^2)A_{m } + A_{m-2}]x^{m}=0$$
- Note, we define $A_{-1}=A_{-2}=0$ to merge the sums

Cancel all these terms out, we have:

$$x^\alpha \sum_{m=0}^\infty [(m^2 +2\alpha m)A_{m} + A_{m-2}]x^m =0$$

And thus, we get an **recursive relationship**:

$$(m^2 +2\alpha m)A_{m} + A_{m-2}=0\implies A_{m}=-\frac{A_{m-2} }{(m^2 +2 \alpha m)}$$

From the initial condition that $A_{-1}=0$, we know all odd terms:

$$A_{1} =- \frac{0}{(m^2 +2\alpha m)}= A_{2m-1} = 0$$

However, $A_{-2}=0$ did not give us any condition, since:

$$A_{0}= \frac{0}{0} = ?$$

is not determined. Therefore, we can choose any conveinent $A_{0}$ to start with.

Look back to the equation, with all odd terms killed, we rewrite:

$$A_{2m} = -\frac{ {A_{2(m-1)} }}{((2m)^2 +4\alpha m)}=-\frac{A_{2(m-1)} }{4m(m+\alpha)}$$

Now, we can plug in some numbers to try this out:

$$\begin{align}
m=1: &&A_{2} &= (-1)\frac{A_{0} }{4(1+\alpha)} \\
m=2:&&A_{4}&= (-1)\frac{A_{2} }{2\cdot 4 (2+\alpha)} = (-1)^2\frac{A_{0} }{4^2\cdot(1\cdot 2) (1+\alpha)(2+\alpha)}  \\
\dots \\
m=k: && A_{2k}&= (-1)^k \frac{A_{0} }{4^k\cdot (k!) \prod_{j=\alpha+1}^{\alpha+k}j}
\end{align}$$

We know that the summation in the final power series equation needs to be convergent in order for it to be a solution to the Bessel equation. In addition, it would be advantageous to use the factorial $(\alpha-m)!$ instead of the ugly prod notation.

Therefore, we here by leveraging our right of choosing the initial term as:

$$A_{0}= \frac{1}{2^\alpha \alpha!}$$

- $2^\alpha$ is just for easier to combine the terms, later you will see.
- $\alpha!$ is to provide a support to remove the ugly prod term

Well, this only works if $\alpha$ is an integer, so we use our old friend the **Gamma function**, and say:

$$A_{0}= \frac{1}{2^\alpha \Gamma(\alpha+1)}$$

Where recall that:

$$\Gamma(\alpha):=\int _{0}^\infty e^{-t} t^{\alpha-1}  \, dt \implies \Gamma(\alpha+1) = \alpha \Gamma(\alpha)$$

- For people interested, Gamma function is derived from **Mellin transform** of **negative exponential**, search them up!

which means Gamma is an extension of factorials to real numbers:

$$\Gamma(\alpha+1)=\alpha!$$

- for integer $\alpha$

Then by putting in our carefully constructed first term.

$$A_{2m} = \frac{(-1)^m}{2^\alpha 4^m m!\Gamma(m+\alpha+1)} $$

Good, we are suspiciously close to the Bessel Function of first kind.

By doing another plugin, we have:

$$f(x) = x^\alpha \sum_{m=0}^\infty A_{m}x^m = x^\alpha \sum_{m=0}^\infty \frac{(-1)^m}{2^\alpha 2^{2m} m!\Gamma(m+\alpha+1)}x^{2m}$$

Regroup the terms a little bit:

$$f(x) = \left( \frac{x}{2} \right)^{\alpha} \sum_{m=0}^\infty \frac{(-1)^m}{m! \Gamma(m+\alpha+1)}  \left( \frac{x}{2} \right)^{2m}=J_{\alpha}(x)$$

Here we go

**Q.E.D-Bessel Function of the First Kind**

#### 3.1.2 Linearity of solutions to Bessel Equations

Suppose $y_1(x)$ and $y_2(x)$ are two solutions to this equation. define: 

$$ y(x) = c_1 y_1(x) + c_2 y_2(x) $$ 

where $c_1$ and $c_2$ are constants.

Substitute $y(x)$ into the Bessel equation:

$$ x^2 \frac{d^2}{dx^2} \left(c_1 y_1(x) + c_2 y_2(x)\right) + x \frac{d}{dx} \left(c_1 y_1(x) + c_2 y_2(x)\right) + \left(x^2 - \nu^2\right) \left(c_1 y_1(x) + c_2 y_2(x)\right) = 0 $$ 

Since differentiation is a linear operator, we can distribute the derivatives:

$$ c_1 \left(x^2 \frac{d^2 y_1}{dx^2} + x \frac{dy_1}{dx} + (x^2 - \alpha^2) y_1\right) + c_2 \left(x^2 \frac{d^2 y_2}{dx^2} + x \frac{dy_2}{dx} + (x^2 - \alpha^2) y_2\right) = 0 $$ 

Given that $y_1(x)$ and $y_2(x)$ are solutions to the Bessel equation, the original equation simplifies to:

$$ c_1 \cdot 0 + c_2 \cdot 0 = 0 $$


Therefore, linear combinations of solutions to the same equation is still a solution
#### 3.1.3 Bessel Function of the Second Kind

![/pics/BesselBasis/BesselSecondKind.png](/pics/BesselBasis/BesselSecondKind.png)

The need for the Bessel Function of Second Kind raised from the the following concerns.

In the differential equation, we noticed that:

$$x^2 \frac{ {\partial^2 f} }{\partial x^2} +x \frac{ {\partial f} }{\partial x} +(x^2-\alpha^2) f =0$$

$\alpha$ is given as $\alpha^2$, which means that:

$$J_{\alpha}(x) , J_{-\alpha}(x)$$

will all be valid solutions to our differential equation. By the linearity of Bessel Equation, we know that:

$$aJ_{\alpha}(x) + bJ_{-\alpha}(x)$$

for any real $a,b$ is still a solution to the equation.

$J_{\alpha}(x)$ and $J_{-\alpha}(x)$ are 2 **linearly independent**. We can examine the behavior for $x\to 0$ under different non-integer $\alpha$, and we find:

$$\begin{align}
\alpha > 0: && &\lim_{ x \to 0 }  \left( \frac{x}{2} \right)^{\alpha} \sum_{m=0}^\infty \frac{(-1)^m}{\Gamma(m+1) \Gamma(m+\alpha+1)}  \left( \frac{x}{2} \right)^{2m} = 0 \\
\alpha <0: &&  &\left( \frac{x}{2} \right)^\alpha = \left( \frac{2}{x} \right)^{-\alpha}\implies\lim_{ x \to 0 }  \left( \frac{x}{2} \right)^{\alpha} \sum_{m=0}^\infty \frac{(-1)^m}{\Gamma(m+1) \Gamma(m+\alpha+1)}  \left( \frac{x}{2} \right)^{2m} = \pm \infty \\
\alpha =0: &&  &\lim_{ x \to 0 }  \left( \frac{x}{2} \right)^{0} \sum_{m=0}^\infty \frac{(-1)^m}{\Gamma(m+1) \Gamma(m+0+1)}  \left( \frac{x}{2} \right)^{2m} = 1+0+0+\dots =1 \\

\end{align}$$

Thus $J_{-\alpha}$ cannot be a scaled multiple of $J_{\alpha}$

However, for $\alpha>0$, this guy has a massive problem:

$$J_{-\alpha}(x) = \left( \frac{x}{2} \right)^{-\alpha} \sum_{m=0}^\infty \frac{(-1)^m}{m! \Gamma(m-\alpha+1)} \left( \frac{x}{2} \right)^{2m}$$

1**Gamma Function Singularity**: The Gamma function, $\Gamma(z)$, has poles at non-positive integers. This means $\Gamma(m - \alpha + 1)$ will be undefined for certain values of $m$ when $\alpha$ is a non-integer positive number. Specifically, for $m = \alpha - 1, \alpha, \alpha + 1, \ldots$, the term $\Gamma(m - \alpha + 1)$ will be undefined, leading to singularities in the series.

2**Asymptotic Behavior**: The series for $J_{\alpha}(x)$ and $J_{-\alpha}(x)$ should converge to finite values for all $x$. The expression above, due to the singularities in $\Gamma(m - \alpha + 1)$, does not ensure this convergence. This is particularly problematic as the terms can become infinitely large or undefined, breaking the validity of the series.

#### 3.1.4 Hankel and Neumann Solutions

The **Hankel Solution** provides the general solution to Bessel function under integer degrees.

For positive integer $\alpha$ we start summing from $\alpha$ and ignore all prior undefined terms:

$$\begin{align}
J_{-\alpha}(x) &= \left( \frac{x}{2} \right)^{-\alpha} \sum_{m=\alpha}^\infty \frac{(-1)^m}{\Gamma(m+1)\Gamma(-\alpha+m+1)} \left( \frac{x}{2} \right)^{2m} \\
&= \left( \frac{x}{2} \right)^{-\alpha} \sum_{m=0}^\infty \frac{(-1)^{m+\alpha} }{\Gamma(\alpha+m+1)\Gamma(-\alpha+\alpha+m+1)} \left( \frac{x}{2} \right)^{2m+2\alpha}  \\
&= (-1)^\alpha\left( \frac{x}{2} \right)^{-\alpha} \left( \frac{x}{2} \right)^{2\alpha} \sum_{m=0}^\infty \frac{(-1)^m}{\Gamma(\alpha+m+1)\Gamma(m+1)}\left( \frac{x}{2} \right)^{2m} \\
&= (-1)^\alpha J_{\alpha}(x)
\end{align}$$

Therefore, we fined that in fact for any integer $\alpha$:

$$J_{\alpha}(x) = (-1)^{\alpha}J_{-\alpha}(x)$$

Holds.

By the linearity condition:

$$(-1)^{\alpha}J_{\alpha}(x) - J_{-\alpha}(x)=0$$

Is also a solution to the Bessel Equations.

Therefore, **Hankel** constructed a following structure:

$$\mathbf{Y}_{n}=\lim_{ \alpha \to n } \frac{ {(-1)^{\alpha}J_{\alpha}(x) -J_{-\alpha}(x)} }{\alpha-n}$$

This limit, when evaluated by L'Hospital's rule, gives the following result:

$$\begin{align}
\mathbf{Y}_{n}=\lim_{ \alpha \to n } \mathbf{Y}_{\alpha} &= \left[(-1)^{\alpha} \frac{ {\partial J_{\alpha}(x)} }{\partial \alpha} - \frac{ {\partial J_{-\alpha}(x)} }{\partial \alpha} \right]_{\alpha=n}
\end{align}$$

In the book *Treatise on the Theory of Bessel Functions (Watson G.A.)* Page 57, this limit is proven to exist by further differentiating. But for us, the process is unnecessarily long. So we take this result for granted in this note.

If we want to extend this result into all $\alpha$, we observe that:

$$\cos(\alpha \pi) = (-1)^\alpha : \alpha\in \mathbb{Z}$$

Therefore, we can replace our linear combination to:

$$\cos (\alpha\pi)J_{\alpha}(x) - J_{-\alpha}(x)$$

For the normalization term below, we know:

$$\sin(\alpha \pi)\to 0:\alpha-n\to 0$$

The **Hankel Solution** for non-integer values is thus:

$$\mathbf{Y}_{\alpha}(x) = 2\pi e^{i\pi \alpha} \frac{ {J_{\alpha}(x)\cos(\alpha \pi) - J_{-\alpha}(x)} }{\sin (2\alpha \pi)}$$

Taking the limit, we see:

$$\begin{align}
\lim_{ \alpha \to n } \mathbf{Y}_{\alpha}(x) &= \lim_{ \alpha \to n } \left[ \frac{ {\alpha-n} }{\alpha-n}  \pi e^{i\pi \alpha}  \frac{ {J_{\alpha}(x) \cos(\alpha \pi) - J_{-\alpha}(x)} }{\sin(\alpha \pi) \cos (\alpha \pi)} \right]  \\
&= \lim_{ \alpha \to n} \left[  \frac{\pi e^{i\pi \alpha} }{\cos (\alpha \pi)} \frac{ {\alpha-n} }{\sin(\alpha \pi)} \frac{ {J_{\alpha}(x)\cos(\alpha \pi) - J_{-\alpha}(x)} }{\alpha-n} \right] \\
&= (-1)^n\lim_{ \alpha \to n } \left[ \frac{ {J_{\alpha}(x)\cos(\alpha \pi) - J_{-\alpha}(x)} }{\alpha-n}\right] \\
&=  \mathbf{Y}_{n}(x) + \lim_{ \alpha \to n} \left[ \frac{(-1)^\alpha \cos(\alpha \pi)-1}{\alpha-n} J_{\alpha}(x) \right] \\
&= \mathbf{Y}_{n}(x)
\end{align}$$

- The third line is because: $\pi e^{i\pi n}=(-1)^n\pi,\cos(n\pi)=(-1)^n$
- Also, $\sin(\alpha \pi) \approx \pi (\alpha-n)\implies (\alpha-n) / \pi(\alpha-n)=1 / \pi$
- The last 2 lines came from the relations between $J_{\alpha}$ and $-J_{\alpha}$

Thus we get a general formed solution working on all $\alpha$

Later, this equation is simplified by **Carl Neumann** in the 19th century as:

$$Y_{\alpha}(x) = \frac{ {\cos(\alpha \pi) J_{\alpha}(x) -J_{-\alpha}(x)} }{\sin(\alpha \pi)}=\frac{\cos(\alpha \pi)}{\pi e^{i\pi \alpha} } \mathbf{Y_{\alpha} }(x)$$

This is the modern **Bessel Function of Second Kind**

And by doing the limit in the same way, we see:

$$Y_{n}(x) = \lim_{ \alpha \to n }\frac{ {\cos(\alpha \pi) J_{\alpha}(x) -J_{-\alpha}(x)} }{\sin(\alpha \pi)} = \frac{1}{\pi} \mathbf{Y}_{n}(x) $$

Therefore, any solutions to Bessel function can be represented as:

$$f(x) = AJ_{\alpha}(x) + BY_{\alpha}(x)$$

#### 3.1.5 Properties of Bessel functions:

Reference: [link](https://math.libretexts.org/Bookshelves/Differential_Equations/Introduction_to_Partial_Differential_Equations_(Herman)/05%3A_Non-sinusoidal_Harmonics_and_Special_Functions/5.05%3A_Fourier-Bessel_Series)

1: **Derivative Identities:**

$$\begin{align}
& \frac{\partial}{\partial x} [x^\alpha J_{\alpha}(x)] = x^\alpha J_{\alpha-1}(x) \\
& \frac{\partial}{\partial x}[x^{-\alpha}J_{\alpha}(x)] = -x^{-\alpha} J_{\alpha+1}(x0)
\end{align}$$

2: **Recursive Formulation:**

$$\begin{align}
&J_{\alpha-1}(x) +J_{\alpha+1}(x) = \frac{2\alpha}{x}J_{\alpha}(x) \\
&J_{\alpha-1}(x) - J_{\alpha+1}(x) = 2J'_{\alpha}(x)
\end{align}$$

3: **Sturm-Liouville Orthogonality:**

$$\int _{0}^a J_{\alpha}\left( j_{\alpha n} \frac{x}{a} \right) J_{\alpha} \left( j_{\alpha m} \frac{x}{a} \right) x\, dx  = \lVert J_{\alpha} \rVert^2 \delta_{m,n}$$

- where $j_{\alpha n}$ is the n-th root for the equation $J_{\alpha}(x)=0$ or $J'_{\alpha}(x)=0$. This property is really important and we will derive it in the next section
### 3.2 Sturm-Liouville Theory of Bessel Functions:
#### 3.2.1 Sturm-Liouville Theory

To make this thing usable, we need something called the **Sturm-Liouville (S-L) theory**, which says that for a problem defined as:

$$L[y] =\frac{d}{dx}\left[ p(x) \frac{dy}{dx} \right] +q(x)y = -\lambda \omega(x)y$$

- $y$ is a unknown function of $x$
- $L$ is the Sturm-Liouville (S-L) differential operator
- $\lambda$ is a eigenvalue, and $y(x)$ is the corresponding eigenfunction
- $p(x),q(x)$ are coefficient functions greater than 0 over the interval of interest.
- $\omega(x)$ is the weight function used in the inner product of the solutions under this problem. $\omega(x)>0$. It is intrpduced from a physics context.
- $y$ is defined on some closed interval $[a,b]$

The **Lagrangeís identity** gives this problem the following boundary conditions ($\alpha$ and $\beta$ are any of our choice from some range defined by the problem):

$$\alpha_{1}y(a) +\alpha_{2}y'(a)=0, \beta_{1}y(b)+\beta_{2}y'(b)=0$$

This can form a set of **Dirichlet** or **Neumann** boundary conditions given a differnt combination of $\alpha$ and $\beta$.

These conditions makes sure that the S-L differential operator is **self-adjoint**, meaning:

$$\langle L[y_{m}],y_{n} \rangle_{\omega} -\langle y_{m}, L[y_{n}] \rangle_{\omega} =  \int _{a}^b (L[y_{m}]y_{n}-y_{m}L[y_{n}]) \, \omega(x) dx = 0$$

These boundary conditions will give a limit on available $\lambda$s

A detailed walkthroughs and proofs can be found here: [link](https://math.libretexts.org/Bookshelves/Differential_Equations/Introduction_to_Partial_Differential_Equations_(Herman)/04%3A_Sturm-Liouville_Boundary_Value_Problems/4.02%3A_Properties_of_Sturm-Liouville_Eigenvalue_Problems)

If all the above condition mets, the **S-L theory** says that:

- The eigenvalues are nonnegative real numbers and can be numbered to form an increasing sequence $\lambda_{1},\dots,\lambda_{n}$
- The corresponding eigenfunctions can be uniquely determined up to a costant multiplier
- The eigenfunctions are mutually orthogonal with respect to the weight function:

- $$\langle y_{m},y_{n} \rangle = \int _{a}^b y_{n}(x) y_{m}(x) \omega(x) \, dx =0 :m\neq n $$

- The n-th eigenfunction has exactly $n-1$ zeros on the interval $[0,a]$
- The complete set of eigenfunctions forms a complete orthogonal set of functions defined on the interval $[0,a]$

And this is a property we need. We need to convert our radial basis into a S-L problem.

#### 3.2.2 Bessel Equation Orthogonality

For easier understanding, we use the standard Bessel Equation as an example:

$$x^2y'' + xy' +(k^2x^2-\alpha^2) = 0\implies y(x)=J_{\alpha}(kx)$$

defined for $x$ between $[0,a]$

rewrite:

$$\frac{d}{dx}\left( x \frac{dy}{dx} \right) + \left( k^2x-\frac{\alpha^2}{x} \right)y=0$$

Now, for 2 different candidates $k_{1}$ and $k_{2}$, we have:

$$y_{1}(x) = J_{\alpha}(k_{1}x)\quad y_{2}(x) = J_{\alpha}(k_{2}x)$$

which solves:

$$\begin{align}
\frac{d}{dx}\left( x \frac{dy_{1} }{dx} \right) + \left( k_{1}^2x-\frac{\alpha^2}{x} \right)y_{1}&=0 & &(1) \\
\frac{d}{dx}\left( x \frac{dy_{2} }{dx} \right) + \left( k_{2}^2x-\frac{\alpha^2}{x} \right)y_{2}&=0 & &(2)
\end{align}$$

Multiply $(1)$ by $y_{2}$ and $(2)$ by $y_{1}$, and subtract $(2)$ from $(1)$, we have:

$$(k_{1}^2 -k_{2}^2)y_{1}y_{2}x = -\frac{d}{dx}\left( x \frac{dy_{1} }{dx} \right)y_{2} + \frac{d}{dx}\left( x \frac{dy_{2} }{dx} \right)y_{1}$$

Integrating both side with respect to $x$ from $[0,a]$, we have:

$$(k_{1}^2 -k_{2}^2) \int _{0}^a  y_{1}y_{2}x \,dx = -\int _{0}^a  \frac{d}{dx}\left( x \frac{dy_{1} }{dx} \right)y_{2}\, dx + \int _{0}^a   \frac{d}{dx}\left( x \frac{dy_{2} }{dx} \right)y_{1} \, dx$$

Apply integration by part on the right for the inner component, we have:

$$\begin{align}
-\int _{0}^a  \frac{d}{dx}\left( x \frac{dy_{1} }{dx} \right)y_{2}\, dx  &= - \left[ x \frac{dy_{1} }{dx} y_{2} \right]_{0}^a + \int _{0}^a x y_{1}' y_{2}'\, dx \\
 \int _{0}^a  \frac{d}{dx}\left( x \frac{dy_{2} }{dx} \right)y_{1}\, dx  &=  \left[ x \frac{dy_{2} }{dx} y_{1} \right]_{0}^a -\int _{0}^a x y_{1}' y_{2}'\, dx
\end{align}$$

The last part cancels out, thus:

$$(k_{1}^2 -k_{2}^2) \int _{0}^a  y_{1}y_{2}x \,dx = - \left[ x \frac{dy_{1} }{dx} y_{2} \right]_{0}^a + \left[ x \frac{dy_{2} }{dx} y_{1} \right]_{0}^a $$

Now, since $x=0$ makes the first term zero, we have:

$$(k_{1}^2 -k_{2}^2) \int _{0}^a  y_{1}y_{2}x \,dx = - \left[ x
\frac{dy_{1} }{dx} y_{2} \right]_{x=a} + \left[ x \frac{dy_{2} }{dx} y_{1} \right]_{x=a}$$

Then we substitute Bessel functions into $y_{1}$ and $y_{2}$, we see:

$$(k_{1}^2 -k_{2}^2) \int _{0}^a  J_{\alpha}(k_{1}x) J_{\alpha}(k_{2}x)x \,dx =   a [k_{2}J_{\alpha}(k_{1}a)J_{\alpha}'(k_{2}a)-k_{1} J_{\alpha}(k_{2}a)J_{\alpha}'(k_{1}a) ]$$

Finally, if we substitute $k_{1}=\frac{j_{1} }{a}$ and $k_{2}=\frac{j_{ 2} }{a}$ as roots for $J_{\alpha}$ or $J_{\alpha}'$, we will see:

$$\frac{j_{1}^2 -j_{2}^2}{a^2}\int _{0}^a  J_{\alpha}\left( \frac{j_{1}x}{a} \right) J_{\alpha}\left( \frac{j_{2}x}{a} \right)x \,dx  = j_{2} J_{\alpha}(j_{1})  J_{\alpha}'(j_{2}) -  j_{1} J_{\alpha}(j_{2})  J_{\alpha}'(j_{1})$$

Note: the specific boundary conditions of S-L gives the kind of roots that we need
#### 3.2.3 Dirichlet Boundary Conditions (Roots of J)

The Dirichlet boundary condition applies on the solution directly. In our case:

$$J_{\alpha}(j_{n})=0$$

Then we will have:

$$\lVert J_{\alpha} \rVert^2 =\lim_{ p \to j_{n} }  \frac{a^2}{p^2 -j_{n}^2} [ j_{n} J_{\alpha}(p)  J_{\alpha}'(j_{n}) -  p J_{\alpha}(j_{n})  J_{\alpha}'(p)]$$

Now, we use L'Hopital's rule and:

$$\lVert J_{\alpha} \rVert^2 =\lim_{ p \to j_{n} }  \frac{a^2}{2p} \frac{d}{dp}[ j_{n} J_{\alpha}(p)  J_{\alpha}'(j_{n}) -  0]=\frac{a^2}{2} [J_{\alpha}'(j_{n})]^2 = \frac{a^2}{2} [J_{\alpha+1}(j_{n})]^2$$

Thus:

$$\int _{0}^a  J_{\alpha}\left( \frac{j_{n}x}{a} \right) J_{\alpha}\left( \frac{j_{m}x}{a} \right)x \,dx = \begin{cases}
\frac{a^2}{2} [J_{\alpha}'(j_{n})]^2 = \frac{a^2}{2} [J_{\alpha+1}(j_{n})]^2 & \text{if m=n} \\
0 & \text{ otherwise}
\end{cases}$$

#### 3.2.4 Neumann Boundary Conditions (Roots of J')

SImilarly:

$$J_{\alpha}'(j_{n})=0$$

following same steps as above, we have:

$$\int _{0}^a  J_{\alpha}\left( \frac{j_{n}x}{a} \right) J_{\alpha}\left( \frac{j_{m}x}{a} \right)x \,dx = \begin{cases}
\frac{a^2}{2} [J_{\alpha}(j_{n})]^2& \text{if m=n} \\
0 & \text{ otherwise}
\end{cases}$$

#### 3.2.5 Combined Boundary Condition

Now, for the roots of the general S-L boundary condition, we have:

$$\beta_{0} a J_{\alpha}(j_{n}) + \beta_{1} j_{n} J_{\alpha}'(j_{n})=0$$

where $\beta_{0}$ and $\beta_{1}$ are real values, $a$ is the range bound, and $j_{n}$ is the root such that this equation holds.

If $n \neq m$, the colinear property of $J_{\alpha}$ and $J'_{\alpha}$ makes the inner product zero.
If $n=m$, after some L'Hopitals we will get:

$$\lVert J_{\alpha} \rVert^2 = \frac{a^2}{2} (\left[ J_{\alpha}'(j_{n})^2 + \left( 1-\frac{\alpha^2}{j_{n}^2}  \right) J^2_{\alpha} (j_{n}) \right]$$

### 3.3 Solving Radial Basis

Now, as we have:

$$(r^2k_{r}^2) \frac{ {\partial^2R} }{\partial r^2} + (rk_{r}) \frac{ {\partial R} }{\partial r} + ((rk_{r})^2 - k_{\varphi}^2)=0$$

The Bessel functions give us the following sets of solutions:

$$R(r;k_{r})= AJ_{k_{\varphi} }(k_{r}r) + BY_{k_{\varphi} }(k_{r}r)$$

Now, because we need the basis to be non-singular at origin, and for polar coordinates $r\geq 0$. Thus, we can entirely ignore the second term, and write:

$$R(r;k_{r})= J_{k_{\varphi} }(rk_{r})$$

To prove this is a basis, we can check its orthogonality, as:

$$\int _{0}^\infty \, J_{k_{\varphi} }(rk_{1}) J_{k_{\varphi} } (rk_{2})r\, dr = \frac{1}{k_{1} }\delta (k_{1}-k_{2})$$

- The $r$ in the inner integral came from the change of variable from $dxdy$ to $dr d\varphi$

However, this seemingly proper function has a massive landmine. Unlike the beautiful complex exponential that holds its orthogonality over any domain, if we take a finite subdomain of Bessel Functions, we see:

$$\begin{align}
\int _{0}^a J_{k_{\varphi} }(rk_{1}) J_{k_{\varphi} } (rk_{2})\, dx = \frac{a}{k_{1}^2 -k_{2}^2} [k_{2}J_{k_\varphi}(k_{1}a)J_{k_\varphi}'(k_{2}a)-k_{1} J_{k_\varphi}(k_{2}a)J_{k_\varphi}'(k_{1}a) ]
\end{align}$$

Therefore, we will use **S-L Theory on Bessel Functions**

Retrieve our initial differential equation:

$$\frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{\partial R}{\partial r} \right)  + \left( k_{r}^2 - \frac{k_{\varphi}^2}{r^2} \right) R = 0$$

By regroup our terms a little bit, we have:

$$-\frac{\partial}{\partial r} \left( r \frac{ {\partial R} }{\partial r} \right) + \frac{k_{\varphi^2} }{r}R = k_{r}^2 r R$$

Writing in the S-L form, we have:

$$\begin{cases}
p(r)=r \\
q(r)=\frac{k_{\varphi}^2}{r} \\
\omega(r)=r \\
\lambda=k_{r}^2
\end{cases}$$

And thus:

$$-\frac{\partial}{\partial r}\left( p(r) \frac{ {\partial R} }{\partial r} \right) + q(r)R=\lambda \omega(r)R$$

with $r\in[0,a]$

now we bring in the boundary conditions:

$$\begin{cases}
R(0) \cos \alpha - p(0)R'(0) \sin \alpha = 0 \\
R(a) \cos \beta - p(a)R'(a) \sin \beta =0
\end{cases}$$

-  $\alpha,\beta \in[0,\pi)$
- The reason for using trigs is just for the ease of isolating one of the terms from the boundary condition to form a desired outcome

This form a complete S-L probelm.

Since we already know that $R(x)=J_{k_{\varphi} }(k_{r}r)$, is a general non-singular solution to the differential equation, we know that $\lambda=k_{r}^2$ as well. This gives us 2 conclusions:

- 1: with $\alpha=\frac{\pi}{2}$ and $r=0$, the first boundary condition is gone
- 2: all our $k_{r}$ selections needs to satisfy the following constraint when $r=a$: 

$$J_{k_{\varphi} } (k_{r}a)\cos \beta -k_{r}aJ'_{k_{\varphi} }(k_{r}a) \sin \beta =0$$

Using a change of variable with $j=k_{r}a$, we have

$$J_{k_{\varphi} }(j)\cos \beta - j J_{k_{\varphi} } (j) \sin \beta=0$$

Now this looks like **mixed boundary condition** problem as we discussed before.

Suppose $j_{k_{\varphi},1} < \dots < j_{k_{\varphi},n}<\dots$ are non negative solutions to this equation with $J_{k_{\varphi} }\left( \frac{j_{k_{\varphi},n}r}{a} \right)\neq 0$, then $k_{r}$ can take values from:

$$\left\{  \frac{j_{k_{\varphi},1} }{a},\dots, \frac{j_{k\varphi,n} }
{a}\dots  \right\}$$

Thus, we can define:

$$k_{r,n,k_\varphi}= \frac{ {j_{k_{\varphi},n} }}{a}$$

Now using our orthogonality results from earlier, we have:

$$\int _{0}^a J_{k_{\varphi} }(k_{r,n,k_\varphi}r) J_{k_{\varphi} }(k_{r,m,k_\varphi}r) r\, dr = \lVert  J_{k_{\varphi} } \rVert^2 \delta_{mn}$$

With:

$$\lVert J_{k_{\varphi} } \rVert^2_{n} = \frac{a^2}{2} (\left[ J_{k_{\varphi} }'(j_{n, k_\varphi})^2 + \left( 1-\frac{\alpha^2}{j_{n,k_\varphi}^2}  \right) J^2_{\alpha} (j_{n,k_\varphi}) \right]$$

And therefore, our **final n-th basis function** looks like this:

$$R_{n,k_{\varphi} }(r;k_{r,n,k_{\varphi} })=J_{k_{\varphi} }(k_{r,n,k_{\varphi} }r)$$

The normalized version will be:

$$R_{n,k_{\varphi} }(r;k_{r,n,k_{\varphi} })= \frac{1}{\sqrt{ \lVert J_{k_{\varphi} } s\rVert _{n}^2 } }J_{k_{\varphi} }(k_{r,n,k_{\varphi} }r)$$

dependent on the specific $\beta$ or boundary condition of our choice.

Now, for each radial function component:

$$f(r) = \sum_{n=1}^\infty \left[ \int _{0}^a f(r_{0}) R_{n,k_{\varphi} }(r_{0})r_{0} \, dr_{0}  \right] R_{n.k_{\varphi} } (r)$$

just like Fourier transform, its a inner product.
#### 3.3.2 Dirichlet (Zero value) Boundary Conditions and Fourier Bessel Series

Like earlier, the first boundary condition give:

$$J_{k_{\varphi} } (k_{r}a)=0$$

Thus:

$$\lVert J_{k_{\varphi} } \rVert^2_{n} = \frac{a^2}{2} J^2_{k_{\varphi}+1} (k_{r,n,k_{\varphi} }a)$$

Under this senario, we can define:

$$c_{n,k_{\varphi} } = \sqrt{  \frac{a^2  J^2_{k_{\varphi}+1} (k_{r,n,k_{\varphi} }a)}{2} }  \left[ \int _{0}^a f(r_{0}) R_{n,k_{\varphi} }(r_{0})r_{0} \, dr_{0}  \right]$$

and write:

$$f(r) = \sum_{n=1}^\infty c_{n,k_{\varphi} } J_{k_{\varphi} } (k_{r,n,k_{\varphi} }r)$$

This is called the **Fourier Bessel Series** of a radial function at degree $k_{\varphi}$

There is also the version of basis for **Neumann Boundary Condition**, but we will not discuss that in here to save some space.

## 4: The Fourier Bessel Transform

### 4.1: Complete Basis Structure

![/pics/BesselBasis/BesselBasis.png](/pics/BesselBasis/BesselBasis.png)

The basis function for the polar Fourier transform is composed of the radial and the angular parts.

For transforming the entire space, the basis is given by:

$$\Psi _{k_{r},k_{\varphi} } = \sqrt{ k_{r} } J_{k_{\varphi} }(k_{r}r) \Phi _{k_\varphi}(\varphi)$$

And for transforming only $r\in[0,a]$, the basis is given by:

$$\Psi _{k_{\varphi},n} = R_{k_{\varphi},n}(r) \Phi _{k_\varphi}(\varphi)$$

- the $n$-th basis corresponding to the $n$-th solution of $J_{k_{\varphi} }(k_{r}a)=0$

The orthogonality follows from orthogonality of Bessel and complex exponentials respectively.

The basis also satisfies our original Helmholtz Equation:

$$\nabla ^2 \Psi_{k_{\varphi},n}(r,\varphi)+k_{r}^2 \Psi_{k_{\varphi},n}(r,\varphi)=0$$

For the basis over finite domain as $\Psi_{k_{\varphi},n}$, the 2 components are:

- $k_{\varphi}=\dots-2,-1,0,1,2\dots$  is "the number of periods in the angular direction. In the polar fourier transform we consider only integer periods.
- $n=1,2,\dots:n-1$  corresponds to the number of zero crossings in the radial direction (As for the $n$-th root, the corresponding $J_{k_{\varphi} }$ has $n-1$ zeros).

Thus, these $n$ and $k_{\varphi}$ values forms a orthnormal basis for all fucntions over the disk, as:

![[Screenshot from 2024-06-01 18-50-34.png]]

### 4.2: The Transformation

Like traditional Fourier transform, the component of a basis in the function is given by the inner product of a basis with the function. Now for a particular period $k_{\varphi}$ and radial basis $k_{r}$, the inner product is:

$$P_{k_{r},k_{\varphi} } = \int _{0}^{2\pi} \int _{0}^\infty f(\mathrm{r,\varphi}) \Psi^*_{k_{\varphi},k_{r} } (r,\varphi) r\, dr  \, d\varphi $$

- This is the **polar Fourier Coefficient** (which like before, is a complex number giving amplitude and phase information. The phase information will only be for angular component)
- $\Psi^*$ is the complex conjugate of the original basis. Its under the same convention as the negative sign for the Fourier transform.

And for the inverse transform, since we only consider angular basis of integer periods, we have:

$$f(r,\varphi)=\int _{0}^\infty \sum_{k_{\varphi}=-\infty}^\infty P_{k_{\varphi},k_{r} } \Phi_{k_{r},k_{\varphi} }(r,\varphi) k_{r}\, dk_{r}$$

And for the finite domain case, we have:

$$P_{k_{\varphi},n} = \int _{0}^{2\pi} \int _{0}^a f(\mathrm{r,\varphi}) \Psi_{k_{\varphi},n}^* (r,\varphi) r\, dr  \, d\varphi $$

Where the inverse transform is, of course:

$$f(r,\varphi)=\sum_{n=1}^\infty \sum_{k_{\varphi}=-\infty}^\infty P_{k_{\varphi},n} \Psi _{k_{\varphi},n}(r,\varphi)$$

And here my friends, is the Fourier Transformation under polar coordinates.

### 4.3: Discrete Fourier Bessel Transform

Recall that our discrete Fourier transfrom looks like this:

$$F= \frac{1}{\sqrt{ 2\pi } }\begin{bmatrix}
\omega_{1}(x_{1}) & \dots & \omega_{1}(x_{k}) \\
&\dots \\
\omega_{n}(x_{1}) & \dots & \omega_{n}(x_{k})
\end{bmatrix} \begin{bmatrix}
f(x_{1}) \\
\dots \\
f(x_{k})
\end{bmatrix}=\frac{1}{\sqrt{ 2\pi } } W_{k} f_{k}$$

Which means, for the 2D Fourier transform, the discretization looks like:

$$F_{k_{x},k_{y} } = \frac{1}{2\pi} (W_{M\times M} f_{M\times N}) W_{N\times N}$$

Where $W_{M\times M}$ applys the transform on the rows, and $W_{N\times N}$ applys the transform on the columns. Then the inverse transform is defined as:

$$f_{M\times N} = \frac{1}{2\pi} (W_{M\times M}^H F_{M\times N}W^H _{N\times N})$$

where $H$ is the conjugate transpose (Hermitian)

In practice, for a signal of size $M\times N$, we use the "normalized" frequencies, where:

$$F=\frac{1}{ \sqrt{ MN } }(W_{M\times M} f_{M\times N}) W_{N\times N}$$

with the terms being:

$$W_{M\times M}: W_{mk} = \omega_{m}(x_{k})= e^{-i 2\pi (x_{k})/M}$$

Similarly for $W_{N\times N}$

In the radial and angular senario, our samples are of the form:

$$(r_{1}-r_{n}, \varphi_{1}-\varphi_{n})$$

Which means its a **radial grid** with more samples near the origin and less samples going outward.

The continuous transform looks like this:

$$P_{k_{\varphi},n} = \int _{0}^{2\pi} \int _{0}^a f(\mathrm{r,\varphi}) \Psi_{k_{\varphi},n}^* (r,\varphi) r\, dr  \, d\varphi $$

With the basis:

$$\begin{align}
&\Psi _{k_{\varphi},n} = R_{k_{\varphi},n}(r) \Phi _{k_\varphi}(\varphi)  \\
&\Phi_{k_{\varphi} }(\varphi) = \frac{1}{\sqrt{ 2\pi } }e^{ik_{\varphi} \varphi} \\

&R_{k_{\varphi},n} = \frac{1}{\sqrt{ \lVert J_{k_{\varphi} } \rVert _{n}^2 } }J_{k_{\varphi} }(k_{r,n,k_{\varphi} }r_{p})
\end{align}$$

We can represent this as the sum:

$$P_{k_{\varphi},n} = \sum_{m=0}^M  \sum_{p=0}^N f(r_{p},\varphi_{m}) \Psi_{k_{\varphi},n}^* (r_{p},\varphi_{m}) \,$$

Now, since the inner sum depends on $k_{\varphi}$, the order of the sums matters.

The inverse transform is also:

$$f(r,\varphi) = \sum_{m=0}^M  \sum_{p=0}^N P_{k_{\varphi},n}(r_{p},\varphi_{m}) \Psi_{k_{\varphi},n} (r_{p},\varphi_{m})$$

This is our discrete transform fornow.
Note: I will soon complete the matrix version of this guy, I am still working on deriving it.
## 6: All references

1: (Qing Wang, Olaf Ronneberger, Hans Burkhardt) *Fourier Analysis in Polar and Spherical Coordinates*: [link](https://lmb.informatik.uni-freiburg.de/papers/download/wa_report01_08.pdf)
- Note: This is the main reference.

2: (Joella Rae Deal) *Basics of Bessel Functions*: [link](https://pdxscholar.library.pdx.edu/cgi/viewcontent.cgi?article=1745&context=honorstheses)

3: (Russell Herman) *Introduction to Partial Differential Equations (Herman)* [link](https://math.libretexts.org/Bookshelves/Differential_Equations/Introduction_to_Partial_Differential_Equations_(Herman))
- Note: chapter 4 and chaper 5 are the main focus

4: (G. N. Watson) *A Treatise On The Theory of Bessel Functions*

