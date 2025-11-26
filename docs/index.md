# Advection-Diffusion Solver

Welcome to the documentation for the 1D advection-diffusion solver. This code was developed by Norm Jones [njones61@gmail.com](mailto:njones61@gmail.com).

## Overview

This project implements a numerical solver for one-dimensional advection-diffusion problems in porous media. The solver is designed to simulate solute transport in a vertical column of soil beneath a lake, where upward groundwater flow interacts with downward diffusion from the lake.

## Features

- **1D Advection-Diffusion Solver**: Numerical solution using finite difference methods
- **Flexible Parameters**: Test various scenarios including:
  - Head differences
  - Concentration differences
  - Dispersion coefficients (molecular diffusion + mechanical dispersion)
  - Porosity
  - Column length
  - Hydraulic conductivity
  - Longitudinal dispersivity

## Documentation

- [Theory](theory.md): Detailed description of the governing equations and finite difference solution methods
- [Example](example.md): Complete working example with code and visualization

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
from solve import solve

# Define input parameters
params = {
    # Physical domain
    'L': 5.0,              # Column length [m]
    'porosity': 0.6,       # Porosity [dimensionless]
    
    # Hydraulic parameters
    'K': 0.01,             # Hydraulic conductivity [m/d]
    'delta_h': 0.2,        # Head difference [m] (positive for upward flow)
    
    # Transport parameters
    'D_m': 0.000175,       # Molecular diffusion coefficient [m²/d]
    'alpha_L': 0.5,        # Longitudinal dispersivity [m]
    
    # Boundary conditions
    'C_lake': 285.0,       # Lake concentration at top [M/L³]
    'C_gw': 20.0,         # Groundwater concentration at bottom [M/L³]
    
    # Initial condition
    'C_init': 20.0,        # Initial concentration [M/L³]
    
    # Numerical parameters
    'N': 100,              # Number of grid points
    'delta_t': 1.0,        # Time step [d]
    't_max': 730.0,        # Maximum time [d]
    
    # Options
    'output_interval': 60.0,  # Save snapshot every N days
}

# Run solver
result = solve(params, method='crank_nicolson', verbose=True)

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

### Running the Example

Run the included test case:

```bash
python main.py
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

