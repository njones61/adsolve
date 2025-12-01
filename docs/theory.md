# 1D Advection-Diffusion Theory

## Problem Description

We consider a one-dimensional advection-diffusion problem in a vertical column of soil beneath a lake. The system consists of:

- **Domain**: A vertical column of porous media (soil/sediment) of length $L$ extending from the lake bottom ($z=0$) downward to the groundwater interface ($z=L$)
- **Flow direction**: Upward hydraulic gradient (groundwater flows upward toward the lake)
- **Boundary conditions**: 
  - Top boundary ($z=0$): Lake water with concentration $C_{lake}$
  - Bottom boundary ($z=L$): Groundwater with concentration $C_{gw}$ (typically $C_{lake} > C_{gw}$)
- **Transport mechanisms**: 
  - Advection: Upward flow carries solute from groundwater toward the lake
  - Diffusion/Dispersion: Molecular diffusion and mechanical dispersion transport solute from high to low concentration regions

The competition between upward advection (carrying low-concentration groundwater) and downward diffusion (from high-concentration lake water) creates a concentration profile that we wish to simulate.

## Governing Equation

The one-dimensional advection-diffusion equation for solute transport in porous media is:

>>$\dfrac{\partial C}{\partial t} = D_{eff} \dfrac{\partial^2 C}{\partial z^2} + v \dfrac{\partial C}{\partial z}$

where:

>>$C(z,t)$ i= the solute concentration [M/L³]<br>
>>$t$ = time [T]<br>
>>$z$ = the vertical coordinate (positive downward) [L]<br>
>>$D_{eff}$ = the effective dispersion coefficient [L²/T]<br>
>>$v$ = the pore water velocity magnitude in the upward direction (positive for upward flow) [L/T]

**Note**: Since $z$ increases downward and flow is upward, the advection term has a positive sign. The velocity $v$ is defined as positive for upward flow.

### Effective Dispersion Coefficient

The effective dispersion coefficient accounts for both molecular diffusion and mechanical dispersion:

>>$D_{eff} = D_m + D_{mech}$

where:

>>$D_m$ = the molecular diffusion coefficient [L²/T]<br>
>>$D_{mech} = \alpha_L v$ = the mechanical dispersion coefficient [L²/T]<br>
>>$\alpha_L$ = the longitudinal dispersivity [L]

Alternatively, we can express the effective dispersion as:

>>$D_{eff} = D_m + \alpha_L v$

### Pore Water Velocity

**Note**: The solver accepts the seepage velocity (pore water velocity) $v$ directly as an input parameter. However, if you have hydraulic parameters instead, you can calculate $v$ as described below.

The pore water velocity is related to the Darcy velocity by:

>>$v = \dfrac{q}{\theta}$

where:

>>$q$ = the Darcy velocity (specific discharge) [L/T]<br>
>>$\theta$ = the porosity (dimensionless)

The Darcy velocity can be calculated from the hydraulic gradient:

>>$q = -K \dfrac{dh}{dz}$

where:

>>$K$ = the hydraulic conductivity [L/T]<br>
>>$h$ = the hydraulic head [L]

For a constant upward gradient, we have:

>>$q = K \dfrac{\Delta h}{L}$

where 

>>$\Delta h$ = the head difference (positive for upward flow).

Combining these relationships, the seepage velocity can be calculated from hydraulic parameters as:

>>$v = \dfrac{K \Delta h}{\theta L}$

This formula allows you to compute $v$ if you have measurements of hydraulic conductivity, head difference, porosity, and column length.

## Boundary Conditions

### Top Boundary ($z=0$)

At the lake-sediment interface, we specify a Dirichlet boundary condition:

>>$C(0,t) = C_{lake}$

This assumes the lake concentration is constant and well-mixed.

### Bottom Boundary ($z=L$)

At the groundwater interface, we specify a Dirichlet boundary condition:

>>$C(L,t) = C_{gw}$

This assumes the groundwater concentration is constant.

## Initial Condition

The initial concentration profile can be specified as:

>>$C(z,0) = C_0(z)$

Common choices include:

- Constant initial concentration: $C_0(z) = C_{init}$
- Linear interpolation: $C_0(z) = C_{lake} + (C_{gw} - C_{lake}) \dfrac{z}{L}$
- Step function or other profiles

## Finite Difference Discretization

### Spatial Discretization

We divide the domain into $N$ equally spaced grid points:

>>$z_i = i \Delta z, \quad i = 0, 1, 2, \ldots, N$

where $\Delta z = L/N$ is the grid spacing.

The concentration at each grid point is denoted as:

>>$C_i(t) = C(z_i, t)$

### Temporal Discretization

Time is discretized as:

>>$t^n = n \Delta t, \quad n = 0, 1, 2, \ldots$

where $\Delta t$ is the time step.

The concentration at grid point $i$ and time step $n$ is denoted as:

>>$C_i^n = C(z_i, t^n)$

### Finite Difference Approximations

#### Second-Order Spatial Derivative (Diffusion Term)

Using central differences:

>>$\dfrac{\partial^2 C}{\partial z^2}\bigg|_{z_i} \approx \dfrac{C_{i+1} - 2C_i + C_{i-1}}{(\Delta z)^2}$

#### First-Order Spatial Derivative (Advection Term)

For the advection term, we have several options:

**Central Difference** (second-order accurate, but can be unstable for high Péclet numbers):

>>$\dfrac{\partial C}{\partial z}\bigg|_{z_i} \approx \dfrac{C_{i+1} - C_{i-1}}{2\Delta z}$

**Upwind Difference** (first-order accurate, stable for advection-dominated problems):

>>$\dfrac{\partial C}{\partial z}\bigg|_{z_i} \approx \dfrac{C_{i+1} - C_i}{\Delta z} \quad \text{(for } v > 0 \text{, upward flow)}$

Since we have upward flow ($v > 0$) and $z$ increases downward, the upwind scheme uses forward differences (information from below, i.e., larger $z$ values).

**Hybrid Scheme** (combines central and upwind based on Péclet number):

>>$\dfrac{\partial C}{\partial z}\bigg|_{z_i} \approx \begin{cases}
\dfrac{C_{i+1} - C_{i-1}}{2\Delta z} & \text{if } Pe < 2 \\
\dfrac{C_{i+1} - C_i}{\Delta z} & \text{if } Pe \geq 2
\end{cases}$

where the Péclet number is: $Pe = \dfrac{v \Delta z}{D_{eff}}$

### Time Discretization

#### Explicit (Forward Euler) Method

>>$\dfrac{C_i^{n+1} - C_i^n}{\Delta t} = D_{eff} \dfrac{C_{i+1}^n - 2C_i^n + C_{i-1}^n}{(\Delta z)^2} + v \dfrac{C_{i+1}^n - C_i^n}{\Delta z}$

Solving for $C_i^{n+1}$:

>>$C_i^{n+1} = C_i^n + \Delta t \left[ D_{eff} \dfrac{C_{i+1}^n - 2C_i^n + C_{i-1}^n}{(\Delta z)^2} + v \dfrac{C_{i+1}^n - C_i^n}{\Delta z} \right]$

**Stability condition**: The explicit method requires:

>>$\Delta t \leq \dfrac{(\Delta z)^2}{2D_{eff} + v \Delta z}$

#### Implicit (Backward Euler) Method

>>$\dfrac{C_i^{n+1} - C_i^n}{\Delta t} = D_{eff} \dfrac{C_{i+1}^{n+1} - 2C_i^{n+1} + C_{i-1}^{n+1}}{(\Delta z)^2} + v \dfrac{C_{i+1}^{n+1} - C_i^{n+1}}{\Delta z}$

This leads to a tridiagonal system of equations:

>>$-a_i C_{i-1}^{n+1} + b_i C_i^{n+1} - c_i C_{i+1}^{n+1} = C_i^n$

where:

>>$a_i = \dfrac{\Delta t}{(\Delta z)^2} D_{eff}$<br>
>>$b_i = 1 + \dfrac{2\Delta t}{(\Delta z)^2} D_{eff} + \dfrac{\Delta t}{\Delta z} v$<br>
>>$c_i = \dfrac{\Delta t}{(\Delta z)^2} D_{eff} + \dfrac{\Delta t}{\Delta z} v$

The implicit method is unconditionally stable but requires solving a linear system at each time step.

#### Crank-Nicolson Method (Recommended)

The Crank-Nicolson method averages the explicit and implicit schemes for second-order accuracy in time:

>>$\dfrac{C_i^{n+1} - C_i^n}{\Delta t} = \dfrac{1}{2} \left[ D_{eff} \dfrac{C_{i+1}^{n+1} - 2C_i^{n+1} + C_{i-1}^{n+1}}{(\Delta z)^2} + v \dfrac{C_{i+1}^{n+1} - C_i^{n+1}}{\Delta z} \right] + \dfrac{1}{2} \left[ D_{eff} \dfrac{C_{i+1}^n - 2C_i^n + C_{i-1}^n}{(\Delta z)^2} + v \dfrac{C_{i+1}^n - C_i^n}{\Delta z} \right]$

This also results in a tridiagonal system and is unconditionally stable with second-order accuracy.

## Solution Algorithm

### For Explicit Method

1. **Initialize**: Set $C_i^0 = C_0(z_i)$ for all $i$
2. **Apply boundary conditions**: $C_0^n = C_{lake}$, $C_N^n = C_{gw}$ for all $n$
3. **Time stepping**: 

>For each time step $n$:<br>

>>For each interior point $i = 1, 2, \ldots, N-1$:<br>

>>>>$C_i^{n+1} = C_i^n + \Delta t \left[ D_{eff} \dfrac{C_{i+1}^n - 2C_i^n + C_{i-1}^n}{(\Delta z)^2} + v \dfrac{C_{i+1}^n - C_i^n}{\Delta z} \right]$

>>Check stability condition

4. **Repeat** until desired time is reached

### For Implicit/Crank-Nicolson Method

1. **Initialize**: Set $C_i^0 = C_0(z_i)$ for all $i$
2. **Time stepping**: 

>For each time step $n$:<br>
>>Set up tridiagonal system with boundary conditions<br>
>>Solve tridiagonal system (using Thomas algorithm or similar)<br>
>>Update $C_i^{n+1}$ for all $i$

3. **Repeat** until desired time is reached

## Key Parameters

The following parameters can be varied to test different scenarios:

1. **Seepage velocity** ($v$): Direct input for pore water velocity (controls advection)
2. **Concentration difference** ($C_{lake} - C_{gw}$): Drives the diffusion process
3. **Dispersion coefficient** ($D_{eff}$): Combines molecular diffusion and mechanical dispersion
4. **Column length** ($L$): Domain size
5. **Longitudinal dispersivity** ($\alpha_L$): Controls mechanical dispersion
6. **Molecular diffusion coefficient** ($D_m$): Base diffusion rate

**Note**: The seepage velocity $v$ can be calculated from hydraulic parameters if needed:
- From Darcy velocity and porosity: $v = q/\theta$
- From hydraulic conductivity, head difference, and porosity: $v = (K \Delta h / L) / \theta$

## Steady-State Solution

At steady state ($\partial C/\partial t = 0$), the governing equation becomes:

>>$D_{eff} \dfrac{d^2 C}{dz^2} + v \dfrac{dC}{dz} = 0$

### Derivation

Let $G = \dfrac{dC}{dz}$. Then the steady ODE is:

>>$D_{eff} \dfrac{dG}{dz} + v G = 0 \quad \Rightarrow \quad \dfrac{dG}{dz} = -\dfrac{v}{D_{eff}} G$

Integrating:

>>$G(z) = G(0)\, e^{-\frac{v}{D_{eff}} z}$

Integrating once more:

>>$C(z) = C(0) + \int_0^z G(\xi)\, d\xi = C(0) + \dfrac{D_{eff}}{v} G(0)\left(1 - e^{-\frac{v}{D_{eff}} z}\right)$

Apply the Dirichlet boundary conditions $C(0) = C_{lake}$, $C(L) = C_{gw}$ to solve for $G(0)$:

>>$C_{gw} = C_{lake} + \dfrac{D_{eff}}{v} G(0)\left(1 - e^{-\frac{v}{D_{eff}} L}\right)
\;\Rightarrow\; G(0) = \dfrac{v}{D_{eff}} \dfrac{C_{gw} - C_{lake}}{1 - e^{-\frac{v}{D_{eff}} L}}$

Substitute back to obtain the closed form:

>>$C(z) = C_{lake} + (C_{gw} - C_{lake}) \dfrac{1 - e^{-Pe \dfrac{z}{L}}}{1 - e^{-Pe}}$

where $Pe = \dfrac{vL}{D_{eff}}$ is the Péclet number for the entire domain.

### Limiting cases

- $Pe \to 0$ (diffusion-dominated): using a first-order expansion, the solution approaches a linear profile

>>$C(z) \approx C_{lake} + (C_{gw} - C_{lake}) \dfrac{z}{L}$

- $Pe \to \infty$ (advection-dominated upward flow): the profile transitions sharply near the top boundary with most of the domain close to $C_{gw}$.

### Use in this package

- The function `solve_ss(params)` computes the steady-state profile directly from the analytical solution above using the same parameter dictionary as the transient solver.
- Time-related inputs such as `delta_t`, `t_max`, `save_history`, and `output_interval` are ignored.
- The spatial grid is the same uniform grid used in the transient solver: $z_i = i\,\Delta z$ with $\Delta z = L/N$. The function returns:
  - `C`: steady-state concentration on the grid
  - `z`: grid coordinates
  - `params`: the parameter dictionary augmented with derived quantities (e.g., `v`, `D_eff`, `Pe_global`)

This analytical steady solution is also useful as a benchmark for verifying the transient numerical methods.

This analytical solution can be used to verify the numerical solution.

## Summary of Required Variables

The following table summarizes all variables required to set up and run the 1D advection-diffusion solver:

### Physical Domain Parameters

| Variable | Symbol | Description | Units | Required |
|----------|--------|-------------|-------|----------|
| Column length | $L$ | Length of the vertical soil column | [L] | Yes |

### Hydraulic Parameters

| Variable | Symbol | Description | Units | Required |
|----------|--------|-------------|-------|----------|
| Seepage velocity | $v$ | Pore water velocity (average linear velocity) | [L/T] | Yes |

**Note**: The seepage velocity $v$ can be calculated from other hydraulic parameters if you have them:
- From Darcy velocity and porosity: $v = \dfrac{q}{\theta}$
- From hydraulic conductivity, head difference, and porosity: $v = \dfrac{K \Delta h / L}{\theta} = \dfrac{K \Delta h}{\theta L}$
  
  where:
  - $K$ = hydraulic conductivity [L/T]
  - $\Delta h$ = head difference (positive for upward flow) [L]
  - $\theta$ = porosity (dimensionless)
  - $q = K \dfrac{\Delta h}{L}$ = Darcy velocity [L/T]

### Transport Parameters

| Variable | Symbol | Description | Units | Required |
|----------|--------|-------------|-------|----------|
| Molecular diffusion coefficient | $D_m$ | Base molecular diffusion rate | [L²/T] | Yes |
| Longitudinal dispersivity | $\alpha_L$ | Characteristic length scale for mechanical dispersion | [L] | Yes |
| Mechanical dispersion coefficient | $D_{mech}$ | Calculated as $D_{mech} = \alpha_L v$ | [L²/T] | Calculated |
| Effective dispersion coefficient | $D_{eff}$ | Total dispersion, calculated as $D_{eff} = D_m + \alpha_L v$ | [L²/T] | Calculated |

### Boundary Conditions

| Variable | Symbol | Description | Units | Required |
|----------|--------|-------------|-------|----------|
| Lake concentration | $C_{lake}$ | Solute concentration at top boundary ($z=0$) | [M/L³] | Yes |
| Groundwater concentration | $C_{gw}$ | Solute concentration at bottom boundary ($z=L$) | [M/L³] | Yes |

### Initial Conditions

| Variable | Symbol | Description | Units | Required |
|----------|--------|-------------|-------|----------|
| Initial concentration profile | $C_0(z)$ | Initial concentration distribution | [M/L³] | Yes |
| Initial concentration (constant) | $C_{init}$ | Constant initial concentration (if used) | [M/L³] | Optional |

### Numerical Discretization Parameters

| Variable | Symbol | Description | Units | Required |
|----------|--------|-------------|-------|----------|
| Number of grid points | $N$ | Number of spatial discretization points | Dimensionless | Yes |
| Grid spacing | $\Delta z$ | Spatial step size, calculated as $\Delta z = L/N$ | [L] | Calculated |
| Time step | $\Delta t$ | Temporal step size | [T] | Yes |
| Total simulation time | $t_{max}$ | Maximum time for simulation | [T] | Yes |

### Derived/Calculated Quantities

| Variable | Symbol | Description | Units | Notes |
|----------|--------|-------------|-------|-------|
| Péclet number (local) | $Pe_{local}$ | Local Péclet number, $Pe = \dfrac{v \Delta z}{D_{eff}}$ | Dimensionless | Used for stability analysis |
| Péclet number (global) | $Pe$ | Global Péclet number, $Pe = \dfrac{vL}{D_{eff}}$ | Dimensionless | Characterizes transport regime |
| Concentration at grid point | $C_i^n$ | Concentration at spatial point $i$ and time step $n$ | [M/L³] | Solution variable |

### Notes

- **Required variables** must be specified by the user to run the simulation.
- **Calculated variables** are automatically computed from required variables.
- **Optional variables** may be used depending on the initial condition specification.
- Units: [L] = length, [T] = time, [M] = mass
- The solver computes the concentration profile $C(z,t)$ throughout the domain as a function of space and time.

