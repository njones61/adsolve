# 1D Advection-Diffusion Solver

A Python solver for one-dimensional advection-diffusion problems in porous media. This solver simulates solute transport in a vertical column of soil beneath a lake, where upward groundwater flow interacts with downward diffusion from the lake.

## Features

- **Transient solver**: Numerical solution using finite difference methods (explicit, implicit, Crank–Nicolson)
- **Steady-state solver**: Closed-form analytical solution with Dirichlet boundaries
- **Multiple Time Stepping Methods**: Explicit, Implicit, and Crank-Nicolson
- **Flexible Parameters**: Test various scenarios including head differences, concentration differences, dispersion coefficients, porosity, column length, etc.
- **Visualization**: Automatic plotting of concentration profile evolution
- **Documentation**: Comprehensive theory documentation with MkDocs

## Installation

```bash
pip install -r requirements.txt
```

## Quick Start

```python
from solve import solve_trans, solve_ss

params = {
    'L': 5.0,              # Column length [m]
    'porosity': 0.6,       # Porosity
    'K': 0.01,             # Hydraulic conductivity [m/d]
    'delta_h': 0.2,        # Head difference [m]
    'D_m': 0.000175,       # Molecular diffusion [m²/d]
    'alpha_L': 0.5,        # Longitudinal dispersivity [m]
    'C_lake': 285.0,       # Lake concentration
    'C_gw': 20.0,          # Groundwater concentration
    'N': 100,              # Number of grid points
    'delta_t': 1.0,        # Time step [d]
    't_max': 730.0,        # Maximum time [d]
    'output_interval': 60.0,  # Snapshot interval [d]
}

# Transient solution
result = solve_trans(params, method='crank_nicolson', verbose=True)

# Steady-state solution (analytical)
result_ss = solve_ss({
    'L': 5.0,
    'porosity': 0.6,
    'K': 0.01,
    'delta_h': 0.2,
    'D_m': 0.000175,
    'alpha_L': 0.5,
    'C_lake': 285.0,
    'C_gw': 20.0,
    'N': 100,
}, verbose=True)
```

Run the example:
```bash
python main_trans.py
python main_ss.py
```

## Documentation

Build and view the documentation:
```bash
mkdocs serve
```

Then open http://127.0.0.1:8000 in your browser.

## License

[Add your license here]

## Author

[Add your name/contact here]

