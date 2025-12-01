# Parameter Calibration

## Overview

The calibration module provides automated parameter estimation for the 1D advection-diffusion model by finding optimal parameter values that minimize the Root Mean Square Error (RMSE) between simulated and observed concentration profiles. This is particularly useful when field measurements are available and you need to determine unknown or uncertain transport parameters.

## Problem Formulation

Given observed concentration measurements $C_{obs}(z_i)$ at depths $z_i$, we seek to find parameter values $\mathbf{\theta}$ (a vector of parameters) that minimize the difference between the simulated steady-state profile $C_{sim}(z_i; \mathbf{\theta})$ and observations.

The objective function is the Root Mean Square Error (RMSE):

>>$\text{RMSE}(\mathbf{\theta}) = \sqrt{\dfrac{1}{n} \sum_{i=1}^{n} \left[ C_{sim}(z_i; \mathbf{\theta}) - C_{obs}(z_i) \right]^2}$

where:

>>$n$ = number of valid observation points<br>
>>$C_{sim}(z_i; \mathbf{\theta})$ = simulated concentration at depth $z_i$ with parameters $\mathbf{\theta}$<br>
>>$C_{obs}(z_i)$ = observed concentration at depth $z_i$

### Parameters That Can Be Calibrated

The following parameters can be optimized (individually or in combination):

- **$L$**: Column length [m]
- **$v$**: Seepage velocity (pore water velocity) [m/d]
- **$D_m$**: Molecular diffusion coefficient [m²/d]
- **$\alpha_L$**: Longitudinal dispersivity [m]

Each parameter can be constrained within specified bounds during optimization.

## Optimization Algorithm

The calibration uses the **L-BFGS-B** (Limited-memory Broyden–Fletcher–Goldfarb–Shanno with Bounds) algorithm from SciPy's `scipy.optimize.minimize`. This is a quasi-Newton optimization method that:

- Handles bound constraints efficiently
- Uses gradient information to converge quickly
- Is well-suited for smooth, differentiable objective functions

### Algorithm Steps

1. **Initialization**: Start with user-provided initial parameter values (or midpoints of bounds if not specified)

2. **Iteration**: At each iteration $k$:
   - Evaluate RMSE at current parameter values $\mathbf{\theta}_k$
   - Compute gradient approximation (via finite differences)
   - Update parameter values using L-BFGS-B update formula
   - Check convergence criteria

3. **Convergence**: Stop when:
   - Gradient norm is below tolerance, or
   - Parameter changes are below tolerance, or
   - Maximum iterations reached

4. **Output**: Return optimal parameters $\mathbf{\theta}^*$ and final RMSE

### Bound Constraints

Each parameter $p$ is constrained to lie within specified bounds:

>>$p_{min} \leq p \leq p_{max}$

The optimizer ensures all parameter values remain within these bounds throughout the search.

## Parameter Sensitivity Analysis

After optimization, the algorithm calculates parameter sensitivities to understand which parameters most strongly influence the RMSE. Three types of sensitivities are computed:

### 1. Absolute Sensitivities

The absolute sensitivity measures the change in RMSE per unit change in parameter:

>>$S_{abs}(p) = \dfrac{\partial \text{RMSE}}{\partial p} \approx \dfrac{\text{RMSE}(p + \Delta p) - \text{RMSE}(p)}{\Delta p}$

where $\Delta p$ is a small perturbation (typically 1% of the parameter value).

**Units**: [RMSE/parameter unit]

**Interpretation**: If $S_{abs}(L) = 0.5$, increasing $L$ by 1 m increases RMSE by 0.5.

### 2. Relative Sensitivities

The relative sensitivity measures the percent change in RMSE per 1% change in parameter:

>>$S_{rel}(p) = \dfrac{\partial \text{RMSE}/\text{RMSE}}{\partial p/p} = \dfrac{p}{\text{RMSE}} \dfrac{\partial \text{RMSE}}{\partial p}$

**Units**: Dimensionless

**Interpretation**: If $S_{rel}(v) = 2.0$, a 1% increase in $v$ causes a 2% increase in RMSE.

### 3. Normalized Sensitivities

Normalized sensitivities are scaled so they sum to 1.0:

>>$S_{norm}(p) = \dfrac{|S_{abs}(p)|}{\sum_{i} |S_{abs}(p_i)|}$

**Interpretation**: Higher values indicate more sensitive parameters. Useful for comparing relative importance when parameters have different units.

### Sensitivity Calculation Method

Sensitivities are calculated using forward finite differences at the optimal parameter values:

1. Perturb each parameter by 1% (or minimum absolute perturbation for very small values)
2. Recalculate RMSE with perturbed parameter
3. Compute sensitivity using finite difference approximation

## Usage

### Basic Example

```python
from solve import solve_ss, optimize_parameters
import pandas as pd

# Define input parameters
params = {
    # Solute specie (for plotting, data access)
    'specie': 'Sulfate',
    
    # Physical domain
    'L': 20.0,              # Column length [m]
    'L_solve': True,        # Optimize this parameter
    'L_min': 5,             # Minimum bound
    'L_max': 50,            # Maximum bound
    
    # Hydraulic parameters
    'v': 0.000667,          # Seepage velocity [m/d]
    'v_solve': True,        # Optimize this parameter
    'v_min': 0.00001,       # Minimum bound
    'v_max': 0.01,          # Maximum bound
    
    # Transport parameters
    'D_m': 0.0000175,       # Molecular diffusion [m²/d]
    'D_m_solve': True,      # Optimize this parameter
    'D_m_min': 0.0000001,   # Minimum bound
    'D_m_max': 0.001,       # Maximum bound
    
    'alpha_L': 0.2,         # Longitudinal dispersivity [m]
    'alpha_L_solve': True,  # Optimize this parameter
    'alpha_L_min': 0.01,    # Minimum bound
    'alpha_L_max': 1.0,     # Maximum bound
    
    # Boundary conditions
    'C_lake': 350.0,        # Lake concentration [M/L³]
    'C_gw': 20.0,           # Groundwater concentration [M/L³]
    
    # Numerical parameters
    'N': 200,               # Number of grid points
}

# Load observed data from Excel
observed_pd = pd.read_excel(
    'conc_data.xlsx',
    sheet_name='DW25',
    usecols='A:AL',
    nrows=27
)
# Create Depth [m] from the second column (assumed in cm), make it the index
depth_cm_series = observed_pd.iloc[:, 1]
observed_pd['Depth [m]'] = depth_cm_series.astype(float) / 100.0
observed_pd = observed_pd.set_index('Depth [m]')

# Perform optimization
opt_result = optimize_parameters(params, observed_pd, verbose=True)

# Extract results
params_optimal = opt_result['params_optimal']
rmse_optimal = opt_result['rmse_optimal']
optimal_values = opt_result['optimal_values']
sensitivities = opt_result['sensitivities']
```

### Optimization Parameters

For each parameter you wish to optimize, include the following in the `params` dictionary:

- **`{param}_solve`**: Boolean flag (`True` to optimize, `False` to keep fixed)
- **`{param}_min`**: Minimum bound for the parameter
- **`{param}_max`**: Maximum bound for the parameter

**Important**: At least one parameter must have `*_solve=True`, otherwise an error will be raised.

### Running the Calibration Script

The `main_calibrate.py` script provides a complete example:

```bash
python main_calibrate.py
```

This will:
1. Load observed data from `conc_data.xlsx`
2. Perform parameter optimization
3. Display iteration history table
4. Calculate and display parameter sensitivities
5. Generate a plot with calibrated profile and optimal parameter values
6. Save the plot as `concentration_profile_calibrated.png`

## Output Interpretation

### Optimization Results

The `optimize_parameters()` function returns a dictionary with:

- **`params_optimal`**: Dictionary with optimal parameter values
- **`rmse_optimal`**: Final RMSE value at optimum
- **`optimization_result`**: SciPy optimization result object
- **`param_names`**: List of optimized parameter names
- **`optimal_values`**: Dictionary of optimal parameter values
- **`iteration_history`**: List of dictionaries with parameter values and RMSE at each iteration
- **`sensitivities`**: Dictionary containing sensitivity analysis results

### Iteration History

The iteration history table shows:
- Iteration number
- RMSE at each iteration
- Parameter values at each iteration

This helps you understand:
- How quickly the algorithm converges
- Whether parameters are approaching bounds
- The path taken through parameter space

### Sensitivity Analysis

Use sensitivity results to:

1. **Identify critical parameters**: Parameters with high normalized sensitivity are most important for model accuracy
2. **Guide data collection**: Focus measurement efforts on sensitive parameters
3. **Assess calibration quality**: If all parameters are highly sensitive, the calibration may be ill-posed
4. **Understand model behavior**: Relative sensitivities show which physical processes dominate

### Example Sensitivity Interpretation

```
Normalized sensitivities:
  L          : 0.45
  v          : 0.30
  D_m        : 0.15
  alpha_L    : 0.10
```

This indicates:
- Column length ($L$) is the most sensitive parameter (45% of total sensitivity)
- Seepage velocity ($v$) is moderately sensitive (30%)
- Diffusion and dispersivity are less critical for this calibration

## Best Practices

### 1. Parameter Bounds

- Set bounds based on physical constraints and prior knowledge
- Avoid bounds that are too narrow (may prevent finding optimum)
- Avoid bounds that are too wide (may allow unrealistic values)

### 2. Initial Values

- Provide reasonable initial guesses when possible
- If uncertain, let the algorithm use midpoint of bounds
- Initial values are automatically clipped to bounds

### 3. Which Parameters to Optimize

- Optimize parameters that are uncertain or unknown
- Keep well-measured parameters fixed
- Consider parameter correlations (e.g., $D_m$ and $\alpha_L$ both affect dispersion)

### 4. Convergence

- Check optimization success: `result['optimization_result'].success`
- Review iteration history for convergence behavior
- If convergence fails, try:
  - Adjusting bounds
  - Changing initial values
  - Using a different optimization method

### 5. Validation

- Verify optimal parameters are physically reasonable
- Check that RMSE is acceptably low
- Compare calibrated profile with observations visually
- Consider cross-validation if multiple observation sets are available

## Mathematical Details

### Steady-State Solution

The calibration uses the analytical steady-state solution:

>>$C(z) = C_{lake} + (C_{gw} - C_{lake}) \dfrac{1 - e^{-Pe \cdot z/L}}{1 - e^{-Pe}}$

where the global Péclet number is:

>>$Pe = \dfrac{v L}{D_{eff}} = \dfrac{v L}{D_m + \alpha_L v}$

This solution is evaluated at observation depths and compared to measurements.

### Interpolation

When observation depths don't exactly match grid points, the simulated profile is interpolated to observation depths using:
- **PCHIP** (Piecewise Cubic Hermite Interpolating Polynomial) when available
- **Linear interpolation** as fallback

PCHIP preserves monotonicity and is more accurate for smooth profiles.

### RMSE Calculation

Only valid observation points are used in RMSE calculation:
- Observations with finite (non-NaN) values
- Depths within the model domain
- Depths where interpolation is valid (no extrapolation)

## API Reference

See the [API documentation](api.md) for detailed function signatures:

- `optimize_parameters()`: Main optimization function
- `calculate_rmse()`: RMSE calculation function
- `calculate_parameter_sensitivities()`: Sensitivity analysis function

## References

- SciPy Optimization: [scipy.optimize.minimize](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html)
- L-BFGS-B Algorithm: Byrd, R. H., Lu, P., Nocedal, J., & Zhu, C. (1995). A limited memory algorithm for bound constrained optimization. *SIAM Journal on Scientific Computing*, 16(5), 1190-1208.

