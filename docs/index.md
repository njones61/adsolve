# Advection-Diffusion Solver

Welcome to the documentation for the 1D advection-diffusion solver.

## Overview

This project implements solvers for one-dimensional advection-diffusion problems in porous media. It supports:
- Transient numerical simulation (finite differences) of concentration profile evolution
- Analytical steady-state solution for the same boundary value problem

The solver is designed to simulate solute transport in a vertical column of soil beneath a lake, where upward groundwater flow interacts with downward diffusion from the lake.

This code was developed by Norm Jones [njones61@gmail.com](mailto:njones61@gmail.com) and is hosted on the following github repo: [https://github.com/njones61/adsolve](https://github.com/njones61/adsolve).

## Features

- **Transient solver**: Numerical solution using finite difference methods (explicit, implicit, Crank–Nicolson)
- **Steady-state solver**: Closed-form analytical solution with Dirichlet boundaries
- **Parameter calibration**: Automated optimization to fit model parameters to observed data
- **Sensitivity analysis**: Quantify parameter importance and model sensitivity
- **Flexible Parameters**: Test various scenarios including:
  - Seepage velocity
  - Concentration differences
  - Dispersion coefficients (molecular diffusion + mechanical dispersion)
  - Column length
  - Longitudinal dispersivity

### Solution modes
- Use `solve_trans(params, method=...)` for time-dependent simulations and snapshots
- Use `solve_ss(params)` for the analytical steady-state profile (ignores time parameters)

## Documentation

- [Theory](theory.md): Detailed description of the governing equations and finite difference solution methods
- [Example](example.md): Complete working example with code and visualization
- [Calibration](calibrate.md): Parameter estimation and sensitivity analysis guide
- [API](api.md): Auto-generated API reference (mkdocstrings)

## Quick Start

### Installation

1. **Install dependencies**:

```bash
pip install -r requirements.txt
```

   Required packages:
   
   - `numpy` - Numerical computations
   - `scipy` - Tridiagonal system solver
   - `matplotlib` - Plotting (optional, for visualization)

### Basic Usage

The solver uses a dictionary-based parameter system. Here's a minimal example:

```python
from solve import solve_trans, solve_ss

# Define input parameters
params = {
    # Physical domain
    'L': 5.0,              # Column length [m]
    
    # Hydraulic parameters
    'v': 0.000667,         # Seepage velocity (pore water velocity) [m/d]
    
    # Transport parameters
    'D_m': 0.000175,       # Molecular diffusion coefficient [m²/d]
    'alpha_L': 0.5,        # Longitudinal dispersivity [m]
    
    # Boundary conditions
    'C_lake': 285.0,       # Lake concentration at top [M/L³]
    'C_gw': 20.0,          # Groundwater concentration at bottom [M/L³]
    
    # Numerical parameters
    'N': 100,              # Number of grid points
    'delta_t': 1.0,        # Time step [d]
    't_max': 730.0,        # Maximum time [d]
    
    # Options
    'output_interval': 60.0,  # Save snapshot every N days
}

# Transient solver
result = solve_trans(params, method='crank_nicolson', verbose=True)

# Access results
C = result['C']           # Final concentration profile
z = result['z']           # Spatial grid
t = result['t']           # Time array
params = result['params'] # Updated parameters (includes calculated values)

# Access snapshots if output_interval was set
if 'C_snapshots' in result:
    C_snapshots = result['C_snapshots']  # List of concentration profiles
    t_snapshots = result['t_snapshots']   # List of snapshot times
```

For the steady-state analytical solution, use:

```python
from solve import solve_ss

params = {
    'L': 5.0,
    'v': 0.000667,         # Seepage velocity (pore water velocity) [m/d]
    'D_m': 0.000175,
    'alpha_L': 0.5,
    'C_lake': 285.0,
    'C_gw': 20.0,
    'N': 100,
}

result_ss = solve_ss(params, verbose=True)
C_ss = result_ss['C']
z = result_ss['z']
```

### Running the Example

Run the included test case:

```bash
python main_trans.py
```

Run the steady-state example:

```bash
python main_ss.py
```

This will:

1. Solve the advection-diffusion problem with the default test parameters
2. Print a summary of results and derived parameters
3. Generate and display a plot showing the concentration profile evolution
4. Save the plot as `concentration_profile.png`

For a complete walkthrough with detailed code explanation and result interpretation, see the [Example](example.md) page.

### Time Stepping Methods

The solver supports three time stepping methods:

- **`'crank_nicolson'`** (recommended): Second-order accurate, unconditionally stable
- **`'implicit'`**: First-order accurate, unconditionally stable
- **`'explicit'`**: First-order accurate, conditionally stable (requires small time step)

### Key Parameters

See the [Theory](theory.md#summary-of-required-variables) documentation for a complete list of all required and optional parameters.

### Output

The solver returns a dictionary containing:

- `'C'`: Final concentration profile (numpy array)
- `'z'`: Spatial grid points (numpy array)
- `'t'`: Time array (numpy array)
- `'params'`: Updated parameter dictionary with calculated values
- `'C_snapshots'`: List of concentration profiles at output intervals (if `output_interval` is set)
- `'t_snapshots'`: List of snapshot times (if `output_interval` is set)

