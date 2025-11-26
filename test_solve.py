"""
Test script to compare all three time stepping methods.

Runs the same problem with explicit, implicit, and Crank-Nicolson methods
and plots the results for comparison.
"""

import numpy as np
import matplotlib.pyplot as plt
from solve import solve, initialize_concentration


def main():
    """Run test case with all three methods."""
    
    # Define input parameters (same as main.py)
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
        
        # Initial condition
        'C_init': 20.0,        # Initial concentration (constant, same as C_gw) [M/L³]
        
        # Numerical parameters
        'N': 100,              # Number of grid points
        'delta_t': 1.0,        # Time step [d]
        't_max': 730.0,        # Maximum time [d] (2 years)
        
        # Options
        'save_history': False,
        'output_interval': 60.0,  # Save snapshot every 60 days
    }
    
    methods = ['explicit', 'implicit', 'crank_nicolson']
    results = {}
    
    print("="*70)
    print("COMPARING TIME STEPPING METHODS")
    print("="*70)
    print(f"Problem: {params['L']} m column, {params['t_max']} days simulation")
    print(f"Time step: {params['delta_t']} days, Grid points: {params['N']+1}")
    print("="*70)
    
    # Run each method
    for method in methods:
        print(f"\nRunning {method.upper()} method...")
        print("-" * 70)
        
        try:
            result = solve(params.copy(), method=method, verbose=True)
            results[method] = result
            
            # Print summary
            C = result['C']
            params_calc = result['params']
            print(f"\n{method.upper()} Results:")
            print(f"  Final concentration at top (z=0):     {C[0]:.2f}")
            print(f"  Final concentration at bottom (z=L): {C[-1]:.2f}")
            print(f"  Final concentration at midpoint:     {C[len(C)//2]:.2f}")
            print(f"  Péclet number (global):              {params_calc['Pe_global']:.4f}")
            
        except Exception as e:
            print(f"ERROR: {method} method failed: {e}")
            results[method] = None
    
    # Create comparison plots
    if len([r for r in results.values() if r is not None]) > 0:
        create_comparison_plots(results, params)
    else:
        print("\nNo successful results to plot.")
    
    return results


def create_comparison_plots(results, params):
    """Create comparison plots for all methods."""
    
    # Get z from results (all should have the same z)
    z = None
    params_calc = None
    for result in results.values():
        if result is not None:
            z = result['z']
            params_calc = result['params']  # Use calculated params
            break
    
    if z is None:
        print("No results available for plotting.")
        return
    
    # Get initial condition using calculated parameters
    C_init = initialize_concentration(params_calc)
    
    if z is None:
        return
    
    # Create figure with subplots
    fig = plt.figure(figsize=(16, 10))
    
    # Plot 1: Final concentration profiles comparison
    ax1 = plt.subplot(2, 2, 1)
    
    # Plot initial condition
    ax1.plot(C_init, z, 'k--', linewidth=1.5, alpha=0.7, label='Initial (t=0)')
    
    # Plot final profiles for each method
    colors = {'explicit': 'blue', 'implicit': 'red', 'crank_nicolson': 'green'}
    linestyles = {'explicit': '-', 'implicit': '--', 'crank_nicolson': '-'}
    
    for method in ['explicit', 'implicit', 'crank_nicolson']:
        if results[method] is not None:
            C = results[method]['C']
            ax1.plot(C, z, color=colors[method], linestyle=linestyles[method], 
                    linewidth=2, label=f'{method.replace("_", " ").title()}')
    
    # Plot boundary concentrations
    ax1.axvline(x=params['C_lake'], color='r', linestyle=':', 
               linewidth=1.5, alpha=0.5, label='Lake boundary')
    ax1.axvline(x=params['C_gw'], color='g', linestyle=':', 
               linewidth=1.5, alpha=0.5, label='Groundwater boundary')
    
    ax1.set_xlabel('Concentration [M/L³]')
    ax1.set_ylabel('Depth z [m]')
    ax1.set_title('Final Concentration Profiles: Method Comparison')
    ax1.legend(loc='best', fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.invert_yaxis()
    
    # Plot 2: Difference from Crank-Nicolson (if available)
    ax2 = plt.subplot(2, 2, 2)
    
    if results['crank_nicolson'] is not None:
        C_ref = results['crank_nicolson']['C']
        
        for method in ['explicit', 'implicit']:
            if results[method] is not None:
                C = results[method]['C']
                diff = C - C_ref
                ax2.plot(diff, z, color=colors[method], linestyle=linestyles[method],
                        linewidth=2, label=f'{method.replace("_", " ").title()} - Crank-Nicolson')
        
        ax2.axvline(x=0, color='k', linestyle='--', linewidth=1, alpha=0.5)
        ax2.set_xlabel('Concentration Difference [M/L³]')
        ax2.set_ylabel('Depth z [m]')
        ax2.set_title('Difference from Crank-Nicolson Method')
        ax2.legend(loc='best', fontsize=9)
        ax2.grid(True, alpha=0.3)
        ax2.invert_yaxis()
    else:
        ax2.text(0.5, 0.5, 'Crank-Nicolson\nnot available\nfor comparison', 
                ha='center', va='center', transform=ax2.transAxes)
        ax2.set_title('Difference from Crank-Nicolson Method')
    
    # Plot 3: Evolution comparison (if snapshots available)
    ax3 = plt.subplot(2, 2, 3)
    
    # Plot initial condition
    ax3.plot(C_init, z, 'k--', linewidth=1, alpha=0.5, label='Initial')
    
    # Plot snapshots for each method (show every other snapshot to avoid clutter)
    snapshot_indices = [0, 2, 4, 6, 8, -1]  # First, middle, and last snapshots
    
    for method in ['explicit', 'implicit', 'crank_nicolson']:
        if results[method] is not None:
            C_snapshots = results[method].get('C_snapshots', None)
            t_snapshots = results[method].get('t_snapshots', None)
            
            if C_snapshots is not None and len(C_snapshots) > 0:
                # Use lighter colors for evolution
                color_map = {'explicit': 'lightblue', 'implicit': 'lightcoral', 
                           'crank_nicolson': 'lightgreen'}
                
                for idx in snapshot_indices:
                    if 0 <= idx < len(C_snapshots):
                        C_snap = C_snapshots[idx]
                        t_snap = t_snapshots[idx]
                        alpha = 0.6 if idx != len(C_snapshots) - 1 else 0.9
                        linewidth = 1 if idx != len(C_snapshots) - 1 else 2
                        
                        if idx == len(C_snapshots) - 1:
                            label = f'{method.replace("_", " ").title()} (final)'
                        else:
                            label = None
                        
                        ax3.plot(C_snap, z, color=color_map[method], 
                               linewidth=linewidth, alpha=alpha, label=label)
    
    # Plot boundary concentrations
    ax3.axvline(x=params['C_lake'], color='r', linestyle=':', 
               linewidth=1, alpha=0.3)
    ax3.axvline(x=params['C_gw'], color='g', linestyle=':', 
               linewidth=1, alpha=0.3)
    
    ax3.set_xlabel('Concentration [M/L³]')
    ax3.set_ylabel('Depth z [m]')
    ax3.set_title('Concentration Evolution Comparison')
    ax3.legend(loc='best', fontsize=8, ncol=2)
    ax3.grid(True, alpha=0.3)
    ax3.invert_yaxis()
    
    # Plot 4: Concentration at midpoint over time (if history available)
    ax4 = plt.subplot(2, 2, 4)
    
    for method in ['explicit', 'implicit', 'crank_nicolson']:
        if results[method] is not None:
            C_snapshots = results[method].get('C_snapshots', None)
            t_snapshots = results[method].get('t_snapshots', None)
            
            if C_snapshots is not None and len(C_snapshots) > 0:
                midpoint_idx = len(z) // 2
                C_midpoint = [C[midpoint_idx] for C in C_snapshots]
                
                ax4.plot(t_snapshots, C_midpoint, 'o-', color=colors[method],
                        linewidth=2, markersize=6, label=f'{method.replace("_", " ").title()}')
    
    ax4.set_xlabel('Time [days]')
    ax4.set_ylabel('Concentration at Midpoint [M/L³]')
    ax4.set_title('Concentration at Midpoint vs Time')
    ax4.legend(loc='best', fontsize=9)
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('method_comparison.png', dpi=150)
    print("\nComparison plot saved to 'method_comparison.png'")
    plt.show()
    
    # Print summary statistics
    print("\n" + "="*70)
    print("COMPARISON SUMMARY")
    print("="*70)
    
    if results['crank_nicolson'] is not None:
        C_ref = results['crank_nicolson']['C']
        print(f"\nReference: Crank-Nicolson method")
        print(f"  Final concentration at midpoint: {C_ref[len(C_ref)//2]:.4f}")
        
        for method in ['explicit', 'implicit']:
            if results[method] is not None:
                C = results[method]['C']
                diff = C - C_ref
                max_diff = np.max(np.abs(diff))
                rms_diff = np.sqrt(np.mean(diff**2))
                print(f"\n{method.replace('_', ' ').title()} method:")
                print(f"  Final concentration at midpoint: {C[len(C)//2]:.4f}")
                print(f"  Maximum difference from Crank-Nicolson: {max_diff:.6f}")
                print(f"  RMS difference from Crank-Nicolson: {rms_diff:.6f}")
    
    print("="*70)


if __name__ == '__main__':
    results = main()

