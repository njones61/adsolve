"""
Main script to compute the 1D advection-diffusion steady-state solution.

Test case: Vertical soil column under a lake with upward groundwater flow.
"""

import numpy as np
from solve import solve_ss
import pandas as pd


def main():
    """Run the steady-state computation."""
    
    # Define input parameters (time-related parameters are ignored)
    params = {

        # Solute specie (for plotting, data access)
        'specie': 'Cl-',

        # Physical domain
        'L': 5.0,              # Column length [m]
        'porosity': 0.6,       # Porosity [dimensionless]

        'plot_depth': 1.0,     # Depth to plot concentration profile. If omitted, plot to L [m].
        
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
    
    # Attempt to read observed concentration data from Excel first
    observed_pd = None
    try:
        # Read the specified range from the given sheet
        # Use header at first row and only first 27 rows and columns A:AL
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
        print("Observed data loaded from 'conc_data.xlsx' (DW25, A1:AL27).")
    except FileNotFoundError:
        print("Observed data file 'conc_data.xlsx' not found. Proceeding without overlay.")
    except Exception as exc:
        print(f"Unable to read observed data (DW25 A1:AL27): {exc}. Proceeding without overlay.")
    
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
        
        # Determine plotting depth
        depth_limit = params.get('plot_depth', params.get('L', None))
        if depth_limit is None:
            depth_limit = z[-1]
        try:
            depth_limit = float(depth_limit)
        except Exception:
            depth_limit = z[-1]
        # Clamp to [0, L]
        L_val = params.get('L', z[-1])
        try:
            L_val = float(L_val)
        except Exception:
            L_val = z[-1]
        if depth_limit <= 0 or depth_limit > L_val:
            depth_limit = L_val
        
        # Slice model arrays to requested depth
        mask = z <= depth_limit
        z_plot = z[mask]
        C_plot = C[mask]
        
        # Plot steady-state profile
        plt.plot(C_plot, z_plot, 'b-', linewidth=2.5, label='Steady-state')
        
        # Overlay observed points if available and specie column exists
        specie_val = params.get('specie')
        specie_label = str(specie_val).strip() if isinstance(specie_val, str) else None
        specie_title = specie_label.title() if specie_label else None
        rmse_value = None
        if observed_pd is not None and specie_label:
            # Match column name case-sensitively first; if not found, try case-insensitive match
            obs_col = None
            if specie_label in observed_pd.columns:
                obs_col = specie_label
            else:
                # Case-insensitive search
                lowered = {c.lower(): c for c in observed_pd.columns}
                if specie_label.lower() in lowered:
                    obs_col = lowered[specie_label.lower()]
            if obs_col is not None:
                try:
                    # Restrict observed points to depth_limit
                    observed_slice = observed_pd[observed_pd.index <= depth_limit]
                    # Drop NaNs in the observation column
                    observed_slice = observed_slice[[obs_col]].dropna()
                    obs_depths = observed_slice.index.values.astype(float)
                    obs_values = observed_slice[obs_col].values.astype(float)
                    
                    # Interpolate simulated profile to observation depths (prefer PCHIP)
                    try:
                        from scipy.interpolate import PchipInterpolator
                        interp_fn = PchipInterpolator(z_plot, C_plot, extrapolate=False)
                        sim_at_obs = interp_fn(obs_depths)
                    except Exception:
                        # Fallback to linear interpolation
                        sim_at_obs = np.interp(obs_depths, z_plot, C_plot)
                    
                    # Compute RMSE over valid pairs
                    valid_mask = np.isfinite(sim_at_obs) & np.isfinite(obs_values)
                    if valid_mask.any():
                        diffs = sim_at_obs[valid_mask] - obs_values[valid_mask]
                        rmse_value = float(np.sqrt(np.mean(diffs**2)))
                    
                    plt.plot(
                        observed_slice[obs_col].values,
                        observed_slice.index.values,
                        linestyle='None',
                        marker='o',
                        markersize=5,
                        markerfacecolor='none',
                        markeredgecolor='k',
                        label=f'Observed ({specie_title})'
                    )
                except Exception as exc:
                    print(f"Could not overlay observed points for '{specie_label}': {exc}")
            else:
                print(f"Selected specie '{specie_label}' not found in observed data columns.")
        elif observed_pd is not None and not specie_label:
            print("No specie provided; skipping observed overlay.")
        
        # Plot boundary concentrations
        plt.axvline(x=params['C_lake'], color='r', linestyle=':', 
                    linewidth=1.5, alpha=0.7, label='Lake boundary')
        plt.axvline(x=params['C_gw'], color='g', linestyle=':', 
                    linewidth=1.5, alpha=0.7, label='Groundwater boundary')
        
        plt.xlabel('Concentration [M/L³]')
        plt.ylabel('Depth z [m]')
        plot_title = '1D Advection-Diffusion: Steady-State Concentration Profile'
        if specie_title:
            plot_title += f' ({specie_title})'
        plt.title(plot_title)
        # Add RMSE label if available
        if rmse_value is not None:
            plt.text(
                0.98, 0.02, f'RMSE = {rmse_value:.3g}',
                transform=plt.gca().transAxes,
                ha='right', va='bottom',
                fontsize=10,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.6, edgecolor='gray')
            )
        plt.legend(loc='best', fontsize=9)
        plt.grid(True, alpha=0.3)
        plt.gca().invert_yaxis()  # Invert so z=0 is at top
        # Apply depth limit on y-axis
        plt.ylim(depth_limit, 0)
        plt.tight_layout()
        plt.savefig('concentration_profile_ss.png', dpi=150)
        plt.show()
        print("\nPlot saved to 'concentration_profile_ss.png'")
    except ImportError:
        print("\nMatplotlib not available. Skipping plot.")
    
    return result


if __name__ == '__main__':
    result = main()


