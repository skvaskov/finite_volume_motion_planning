The code in this repository explores the use of finite volume methods for risk-bounded motion planning for mobile robots. This repository contains code to solve the advection equation, which, given a dynamic model, describes the evolution of the probability density of the robot's state as a function of time. The advection equation for a system with states <img src="/tex/ebbad7f1f9841ef914dc318c0c9fad1d.svg?invert_in_darkmode&sanitize=true" align=middle width=50.02787624999999pt height=27.91243950000002pt/> and dynamics <img src="/tex/5721f667885c7fe08087ae05e581a82e.svg?invert_in_darkmode&sanitize=true" align=middle width=63.31043069999999pt height=24.65753399999998pt/> can be written as:

<p align="center"><img src="/tex/92df0a93d533ed7bb156485418d92643.svg?invert_in_darkmode&sanitize=true" align=middle width=479.89774965pt height=47.1348339pt/></p>

where <img src="/tex/2ec6e630f199f589a2402fdf3e0289d5.svg?invert_in_darkmode&sanitize=true" align=middle width=8.270567249999992pt height=14.15524440000002pt/> is the probability density. Note that at present the dynamics are deterministic and autonomous. Stochastic uncertainty and time-varying dynamics can be considered in the future.

We will assume that we are given:
1. The initial probability density of the ego robot, <img src="/tex/8f8a04a5b29f50bc1487384702613871.svg?invert_in_darkmode&sanitize=true" align=middle width=45.97608014999999pt height=24.65753399999998pt/>
2. Probability densities associated with other obstacles or entities <img src="/tex/6936901ed1a2c5709435e83043e866f6.svg?invert_in_darkmode&sanitize=true" align=middle width=103.37176244999999pt height=22.648391699999998pt/> at any time
3. Trajectory dynamics <img src="/tex/190083ef7a1625fbc75f243cffb9c96d.svg?invert_in_darkmode&sanitize=true" align=middle width=9.81741584999999pt height=22.831056599999986pt/>, and a risk function <img src="/tex/232059c0bc77776e7f7420da0c57e44d.svg?invert_in_darkmode&sanitize=true" align=middle width=115.50066659999997pt height=22.648391699999998pt/> that we want to evaluate at discrete time steps

The approach of the finite volume method is to divide the state space, <img src="/tex/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode&sanitize=true" align=middle width=14.908688849999992pt height=22.465723500000017pt/> into cells and keep track of the average volume in each cell. When we solve the advection equation we will ensure that, aside from leakage at the boundaries of <img src="/tex/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode&sanitize=true" align=middle width=14.908688849999992pt height=22.465723500000017pt/>, the volume of the probability density will be conserved. Additionally this numerical method used is monotonic, meaning it will not introduce oscillations as a result of the state space/time discretization. These two properties make the finite volume method attractive compared to other numerical methods used to solve partial differential equations.

We will use the following notation 

<img src="/tex/7c5bf5474c8b8782731a2b814baca441.svg?invert_in_darkmode&sanitize=true" align=middle width=85.02767954999999pt height=27.91243950000002pt/> state space of robot, dimension d.

<img src="/tex/f3d25c125f1b8440a6c51bec1756b687.svg?invert_in_darkmode&sanitize=true" align=middle width=82.71035849999998pt height=27.91243950000002pt/> dynamics of trajectory we want to verify

<img src="/tex/1afcdb0f704394b16fe85fb40c45ca7a.svg?invert_in_darkmode&sanitize=true" align=middle width=12.99542474999999pt height=22.465723500000017pt/> vector containing the average value of the density in each cell

<img src="/tex/55a049b8f161ae7cfeb0197d75aff967.svg?invert_in_darkmode&sanitize=true" align=middle width=9.86687624999999pt height=14.15524440000002pt/> time step number

<img src="/tex/5a63739e01952f6a63389340c037ae29.svg?invert_in_darkmode&sanitize=true" align=middle width=19.634768999999988pt height=22.465723500000017pt/> time step length

<img src="/tex/2f118ee06d05f3c2d98361d9c30e38ce.svg?invert_in_darkmode&sanitize=true" align=middle width=11.889314249999991pt height=22.465723500000017pt/> time horizon

<img src="/tex/16fcdc2e4358a7f49fe2a486169a4238.svg?invert_in_darkmode&sanitize=true" align=middle width=29.198168999999993pt height=22.465723500000017pt/> width of cell in dimension <img src="/tex/941fed4982cce9d21aff5f034342c257.svg?invert_in_darkmode&sanitize=true" align=middle width=89.32558634999998pt height=24.65753399999998pt/>

<img src="/tex/02c8ac95474dd12996d5f96a5fb083fd.svg?invert_in_darkmode&sanitize=true" align=middle width=53.47633664999999pt height=18.666631500000015pt/> the average velocity at the left border of a cell <img src="/tex/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode&sanitize=true" align=middle width=5.663225699999989pt height=21.68300969999999pt/> in dimension <img src="/tex/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode&sanitize=true" align=middle width=7.710416999999989pt height=21.68300969999999pt/>

<img src="/tex/6cf3c7864d3d909d328c1285f1f5a052.svg?invert_in_darkmode&sanitize=true" align=middle width=72.8079231pt height=22.831056599999986pt/> a flux limiter function (we use van leer presently, will be updated to include more)


Let <img src="/tex/87c836923d5a3399b76cae2a1dac8b49.svg?invert_in_darkmode&sanitize=true" align=middle width=21.121448699999988pt height=22.465723500000017pt/> be the cell averages at time step <img src="/tex/55a049b8f161ae7cfeb0197d75aff967.svg?invert_in_darkmode&sanitize=true" align=middle width=9.86687624999999pt height=14.15524440000002pt/>. To compute the cell averages at time step <img src="/tex/3f18d8f60c110e865571bba5ba67dcc6.svg?invert_in_darkmode&sanitize=true" align=middle width=38.17727759999999pt height=21.18721440000001pt/> we use the following algorithm

for <img src="/tex/695e0ecd9088c76cddd0f52dc4a1151a.svg?invert_in_darkmode&sanitize=true" align=middle width=60.101642699999985pt height=22.831056599999986pt/>
	
<p align="center"><img src="/tex/95b73caad661823868db97c2605573db.svg?invert_in_darkmode&sanitize=true" align=middle width=548.80663365pt height=269.51802734999995pt/></p>
	
	
end

Then, at each timestep, we can evaluate the risk function on each cell and sum across all cells
<img src="/tex/215eefb9fd2bbd9546d3a472f3d852f5.svg?invert_in_darkmode&sanitize=true" align=middle width=223.6428975pt height=26.48417309999999pt/>

