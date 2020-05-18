The code in this repository explores the use of finite volume methods for risk-bounded motion planning for mobile robots. This repository contains code to solve the advection equation, which, given a dynamic model, describes the evolution of the probability density of the robot's state as a function of time. The advection equation for a system with states $x\inX\subset\mathbb{R}^d$ and dynamics $\dot{x}=f(x)$ can be written as:

\begin{align}
\frac{\partial}{\partial t}p(x,t)+\sum_{j=1}^n\frac{\partial}{\partial x_j}(f_j(x)p(x,t))=0
\end{align}

where $p$ is the probability density. Note that at present the dynamics are deterministic and autonomous. Stochastic uncertainty and time-varying dynamics can be considered in the future.

We will assume that we are given:
1. The initial probability density of the ego robot, $p(x,0)$
2. Probability densities associated with other obstacles or entities $x_o\in X_o\subset \mathbb{R}^_{d_o}$ at any time
3. Trajectory dynamics $f$, and a risk function $g:X\times X_o \to \mathbb{R}$ that we want to evaluate at discrete time steps

The approach of the finite volume method is to divide the state space, $X$ into cells and keep track of the average volume in each cell. When we solve the advection equation we will ensure that, aside from leakage at the boundaries of $X$, the volume of the probability density will be conserved. Additionally this numerical method used is monotonic, meaning it will not introduce oscillations as a result of the state space/time discretization. These two properties make the finite volume method attractive compared to other numerical methods used to solve partial differential equations.

We will use the following notation 

$x\in X\subset \mathbb{R}^d$ state space of robot, dimension d.

$f:X\to \mathbb{R}^d$ dynamics of trajectory we want to verify

$Q$ vector containing the average value of the density in each cell

$n$ time step number

$\Delta t$ time step length

$T$ time horizon

$\Delta x_j$ width of cell in dimension $j\in\{1,...,d\}$

$\bar{u}_{j,i-1/2}$ the average velocity at the left border of a cell $i$ in dimension $j$

$\phi:\mathbb{R}\to\mathbb{R}$ a flux limiter function (we use van leer presently, will be updated to include more)


Let $Q^n$ be the cell averages at time step $n$. To compute the cell averages at time step $n+1$ we use the following algorithm

for $j = 1:d$
	
\begin{align}
\Delta Q_{i-1/2}^n&=Q_i^n-Q_{i-1}^n\\
I&= \begin{cases} i-1 & \text{ if } \bar{u}_{j,i-1/2} \geq 0\\ i+1 & \text{ if } \bar{u}_{j,i-1/2} <0 \end{cases}\\
\theta_{i-1/2} &= \frac{\Delta Q_{I-1/2}}{\Delta Q_{i-1/2}}\\
\delta_{i-1/2}^n&=\phi(\theta_{i-1/2})\Delta Q_{i-1/2}\\
\begin{split}
F_{i-1/2}^n &=\min(\bar{u}_{j,i-1/2},0)Q_{i}^n+\max(\bar{u}_{j,i+1/2},0)Q_{i-1}^n+\\
&+\frac{1}{2}|\bar{u}_{j,i-1/2}|\left(1-\left |\frac{\bar{u}_{j,i-1/2}\Delta t}{\Delta x_j}\right|\right)\delta_{i-1/2}^n
\end{split}\\
Q_i^{n+1} &= Q_i^n-\frac{\Delta t}{\Delta x_j}(F_{i+1/2}^n-F_{i-1/2}^n)
\end{align}
	
	
end

Then, at each timestep, we can evaluate the risk function on each cell and sum across all cells
$\sum_i Q_i^n \int g(x_i,x_o)p(x_o|t^n)dx dx_o$

Note that some of the examples in this repository use gaussian representations of obstacles and require the matlab package ftnorm to evaluate the moments, which can be found at:
https://www.tandfonline.com/doi/suppl/10.1080/10618600.2017.1322092
