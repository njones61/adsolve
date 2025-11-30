"""
Main script to compute the 1D advection-diffusion steady-state solution.

Test case: Vertical soil column under a lake with upward groundwater flow.
"""

import numpy as np
from solve import solve_ss


def main():
    """Run the steady-state computation."""
    
    # Define input parameters (time-related parameters are ignored)
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
        'C_gw': 20.0,          # Groundwater concentration at bottom [M/L³]
        
        # Numerical parameters
        'N': 100,              # Number of grid points
        
        # Time parameters (ignored)
        'delta_t': 1.0,
        't_max': 730.0,
    }
    
    # Compute steady-state solution
    result = solve_ss(params, verbose=True)
    
    # Extract results
    C = result['C']
    z = result['z']
    params = result['params']
    
    # Print summary
    print("\n" + "="*60)
    print("STEADY-STATE SOLUTION SUMMARY")
    print("="*60)
    print(f"Concentration at top (z=0):     {C[0]:.2f}")
    print(f"Concentration at bottom (z=L): {C[-1]:.2f}")
    print(f"Concentration at midpoint:     {C[len(C)//2]:.2f}")
    print(f"\nDerived parameters:")
    print(f"  Darcy velocity (q):           {params['q']:.6f} m/d")
    print(f"  Pore water velocity (v):      {params['v']:.6f} m/d")
    print(f"  Effective dispersion (D_eff): {params['D_eff']:.6f} m²/d")
    print(f"  Péclet number (global):      {params['Pe_global']:.4f}")
    print("="*60)
    
    # Optional: Plot results (if matplotlib is available)
    try:
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(6, 10))
        
        # Plot steady-state profile
        plt.plot(C, z, 'b-', linewidth=2.5, label='Steady-state')
        
        # Plot boundary concentrations
        plt.axvline(x=params['C_lake'], color='r', linestyle=':', 
                    linewidth=1.5, alpha=0.7, label='Lake boundary')
        plt.axvline(x=params['C_gw'], color='g', linestyle=':', 
                    linewidth=1.5, alpha=0.7, label='Groundwater boundary')
        
        plt.xlabel('Concentration [M/L³]')
        plt.ylabel('Depth z [m]')
        plt.title('1D Advection-Diffusion: Steady-State Concentration Profile')
        plt.legend(loc='best', fontsize=9)
        plt.grid(True, alpha=0.3)
        plt.gca().invert_yaxis()  # Invert so z=0 is at top
        plt.tight_layout()
        plt.savefig('concentration_profile_ss.png', dpi=150)
        plt.show()
        print("\nPlot saved to 'concentration_profile_ss.png'")
    except ImportError:
        print("\nMatplotlib not available. Skipping plot.")
    
    return result


if __name__ == '__main__':
    result = main()


