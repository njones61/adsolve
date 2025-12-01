"""
Parameter sensitivity analysis for the 1D advection-diffusion solver.

Tests the impact of changing a selected parameter by running multiple simulations
and plotting the concentration profiles at t_max for comparison.
"""

import numpy as np
import matplotlib.pyplot as plt
from solve import solve_trans


def run_sensitivity_analysis(param_name, param_min, param_max, n_runs, 
                             base_params=None, method='crank_nicolson', verbose=False):
    """
    Run sensitivity analysis by varying a single parameter.
    
    Parameters
    ----------
    param_name : str
        Name of the parameter to vary (e.g., 'v', 'D_m')
    param_min : float
        Minimum value for the parameter
    param_max : float
        Maximum value for the parameter
    n_runs : int
        Number of runs (parameter values will be evenly spaced)
    base_params : dict, optional
        Base parameter dictionary. If None, uses default values from main.py
    method : str, optional
        Time stepping method (default: 'crank_nicolson')
    verbose : bool, optional
        Print progress for each run (default: False)
        
    Returns
    -------
    dict
        Dictionary containing:
        - 'param_values': array of tested parameter values
        - 'C_profiles': list of concentration profiles at t_max
        - 'z': spatial grid
        - 'param_name': name of the varied parameter
    """
    # Use default parameters if not provided
    if base_params is None:
        base_params = {
            # Physical domain
            'L':5.0,              # Column length [m]
            
            # Hydraulic parameters
            'v': 0.000667,          # Seepage velocity (pore water velocity) [m/d]
            
            # Transport parameters
            'D_m': 0.000175,         # Molecular diffusion coefficient [m²/d]
            'alpha_L': 0.5,          # Longitudinal dispersivity [m]
            
            # Boundary conditions
            'C_lake': 285.0,         # Lake concentration at top [M/L³]
            'C_gw': 20.0,            # Groundwater concentration at bottom [M/L³]
            
            # Initial condition
            'C_init': 20.0,          # Initial concentration (constant, same as C_gw) [M/L³]
            
            # Solute specie
            'specie': 'chloride',
            
            # Numerical parameters
            'N': 100,                # Number of grid points
            'delta_t': 1.0,          # Time step [d]
            't_max': 730.0,          # Maximum time [d] (2 years)
            
            # Options
            'save_history': False,   # Save full time history
            'output_interval': None, # Don't save snapshots for sensitivity analysis
        }
    
    # Generate parameter values
    param_values = np.linspace(param_min, param_max, n_runs)
    
    # Storage for results
    C_profiles = []
    z = None
    
    print(f"\n{'='*60}")
    print(f"PARAMETER SENSITIVITY ANALYSIS")
    print(f"{'='*60}")
    print(f"Parameter: {param_name}")
    print(f"Range: [{param_min:.4f}, {param_max:.4f}]")
    print(f"Number of runs: {n_runs}")
    print(f"{'='*60}\n")
    
    # Run simulations
    for i, param_val in enumerate(param_values):
        if verbose:
            print(f"Run {i+1}/{n_runs}: {param_name} = {param_val:.4f}")
        
        # Create parameter dictionary for this run
        params = base_params.copy()
        params[param_name] = param_val
        
        # Run solver
        result = solve_trans(params, method=method, verbose=False)
        
        # Store final concentration profile
        C_profiles.append(result['C'].copy())
        
        # Store spatial grid (same for all runs)
        if z is None:
            z = result['z']
        
        if verbose:
            print(f"  Complete. Final concentration at top: {result['C'][0]:.2f}")
    
    print(f"\nAll {n_runs} runs completed!")
    
    return {
        'param_values': param_values,
        'C_profiles': C_profiles,
        'z': z,
        'param_name': param_name,
        't_max': base_params['t_max'],
        'specie': base_params.get('specie', 'chloride')
    }


def plot_sensitivity_results(results, save=False):
    """
    Plot concentration profiles for all parameter values.
    
    Parameters
    ----------
    results : dict
        Results dictionary from run_sensitivity_analysis
    save_path : str, optional
        Path to save the plot. If None, displays the plot.
    """
    param_values = results['param_values']
    C_profiles = results['C_profiles']
    z = results['z']
    param_name = results['param_name']
    t_max = results['t_max']
    specie = str(results.get('specie', 'chloride')).title()
    
    # Create figure
    plt.figure(figsize=(6, 10))
    
    # Use a colormap for different parameter values
    colors = plt.cm.viridis(np.linspace(0, 1, len(param_values)))
    
    # Plot each profile
    for i, (C, param_val) in enumerate(zip(C_profiles, param_values)):
        plt.plot(C, z, color=colors[i], linewidth=1.5, 
                label=f'{param_name} = {param_val:.4f}')
    
    plt.xlabel('Concentration [M/L³]', fontsize=12)
    plt.ylabel('Depth z [m]', fontsize=12)
    plt.title(f'Parameter Sensitivity ({specie}): {param_name}\nConcentration Profile at t = {t_max:.0f} days', 
              fontsize=14, fontweight='bold')
    plt.legend(loc='best', fontsize=9, ncol=1, framealpha=0.9)
    plt.grid(True, alpha=0.3)
    plt.gca().invert_yaxis()  # Invert so z=0 is at top
    plt.tight_layout()
    
    save_path = f'param_sensitivity_{param_name}.png'
    if save:
        plt.savefig(save_path, dpi=150)
        print(f"\nPlot saved to '{save_path}'")
    
    plt.show()


def main():
    """Run sensitivity analysis for a selected parameter."""
    
    # Test case: vary D_m parameter
    results = run_sensitivity_analysis(
        param_name='D_m',
        param_min=0.00001,
        param_max=0.001,
        n_runs=10,
        method='crank_nicolson',
        verbose=True
    )
    
    # Plot results
    plot_sensitivity_results(results, save=True)
    
    return results


if __name__ == '__main__':
    results = main()

