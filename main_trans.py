"""
Main script to run the 1D advection-diffusion solver (transient).

Test case: Vertical soil column under a lake with upward groundwater flow.
"""

import numpy as np
from solve import solve_trans, initialize_concentration


def main():
    """Run the transient test case."""
    
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
        'C_gw': 20.0,          # Groundwater concentration at bottom [M/L³]
        
        # Initial condition (ignored by initializer; default is robust)
        'C_init': 20.0,
        
        # Solute specie
        'specie': 'Cl-',
        
        # Numerical parameters
        'N': 100,              # Number of grid points
        'delta_t': 1.0,        # Time step [d]
        't_max': 730.0,        # Maximum time [d] (2 years)
        
        # Options
        'save_history': False,   # Save full time history (set True to save all time steps)
        'output_interval': 60.0, # Save snapshot every N days (None to disable)
    }
    
    # Run transient solver with Crank-Nicolson method (recommended)
    result = solve_trans(params, method='crank_nicolson', verbose=True)
    
    # Extract results
    C = result['C']
    z = result['z']
    t = result['t']
    params = result['params']
    
    # Get snapshots if available
    C_snapshots = result.get('C_snapshots', None)
    t_snapshots = result.get('t_snapshots', None)
    
    # Print summary
    print("\n" + "="*60)
    print("TRANSIENT SOLUTION SUMMARY")
    print("="*60)
    print(f"Final time: {t[-1]:.2f} days")
    print(f"Concentration at top (z=0):     {C[0]:.2f}")
    print(f"Concentration at bottom (z=L): {C[-1]:.2f}")
    print(f"Concentration at midpoint:     {C[len(C)//2]:.2f}")
    print(f"\nDerived parameters:")
    print(f"  Darcy velocity (q):           {params['q']:.6f} m/d")
    print(f"  Pore water velocity (v):      {params['v']:.6f} m/d")
    print(f"  Effective dispersion (D_eff): {params['D_eff']:.6f} m²/d")
    print(f"  Péclet number (global):      {params['Pe_global']:.4f}")
    if C_snapshots is not None:
        print(f"\nSnapshots saved at {len(t_snapshots)} time points")
    print("="*60)
    
    # Optional: Plot results (if matplotlib is available)
    try:
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(6, 10))
        
        # Determine plotting depth
        depth_limit = params.get('plot_depth', params.get('L', None))
        if depth_limit is None:
            depth_limit = z[-1]
        try:
            depth_limit = float(depth_limit)
        except Exception:
            depth_limit = z[-1]
        L_val = params.get('L', z[-1])
        try:
            L_val = float(L_val)
        except Exception:
            L_val = z[-1]
        if depth_limit <= 0 or depth_limit > L_val:
            depth_limit = L_val
        mask = z <= depth_limit
        z_plot = z[mask]
        
        # Plot initial condition (from initializer)
        C_init = initialize_concentration(params)
        plt.plot(C_init[mask], z_plot, 'k--', linewidth=1.5, alpha=0.7, label='Initial (t=0)')
        
        # Plot snapshots if available
        if C_snapshots is not None and len(C_snapshots) > 0:
            colors = plt.cm.viridis(np.linspace(0, 1, len(C_snapshots)))
            for i, (C_snap, t_snap) in enumerate(zip(C_snapshots, t_snapshots)):
                if i == len(C_snapshots) - 1:
                    plt.plot(C_snap[mask], z_plot, color=colors[i], linewidth=2.5, 
                             label=f't = {t_snap:.0f} d (final)')
                else:
                    plt.plot(C_snap[mask], z_plot, color=colors[i], linewidth=1.5, 
                             alpha=0.8, label=f't = {t_snap:.0f} d')
        else:
            plt.plot(C[mask], z_plot, 'b-', linewidth=2, label=f'Final (t={t[-1]:.0f} d)')
        
        # Plot boundary concentrations
        plt.axvline(x=params['C_lake'], color='r', linestyle=':', 
                    linewidth=1.5, alpha=0.7, label='Lake boundary')
        plt.axvline(x=params['C_gw'], color='g', linestyle=':', 
                    linewidth=1.5, alpha=0.7, label='Groundwater boundary')
        
        plt.xlabel('Concentration [M/L³]')
        plt.ylabel('Depth z [m]')
        specie_name = str(params.get('specie', 'chloride')).title()
        plt.title(f'1D Advection-Diffusion: Transient Concentration Profile ({specie_name})')
        plt.legend(loc='best', fontsize=9)
        plt.grid(True, alpha=0.3)
        plt.gca().invert_yaxis()  # Invert so z=0 is at top
        plt.ylim(depth_limit, 0)
        plt.tight_layout()
        plt.savefig('concentration_profile.png', dpi=150)
        plt.show()
        print("\nPlot saved to 'concentration_profile.png'")
    except ImportError:
        print("\nMatplotlib not available. Skipping plot.")
    
    return result


if __name__ == '__main__':
    result = main()


