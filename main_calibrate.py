"""
Calibration script to find optimal parameters that minimize RMSE.

Test case: Vertical soil column under a lake with upward groundwater flow.
"""

import numpy as np
from solve import solve_ss, optimize_parameters
import pandas as pd


def main():
    """Run the calibration/optimization."""
    
    # Define input parameters (time-related parameters are ignored)
    params = {

        # Solute specie (for plotting, data access)
        'specie': 'Na+',

        'plot_depth': 2.0,     # Depth to plot concentration profile. If omitted, plot to L [m].

        # Physical domain
        'L': 5.0,              # Column length [m]
        'L_solve': True,
        'L_min': 2,
        'L_max': 50,

        # Hydraulic parameters
        'v': 0.000667,         # Seepage velocity (pore water velocity) [m/d]
        'v_solve': True,
        'v_min': 0.00001,
        'v_max': 0.01,

        # Transport parameters
        'D_m': 1.33e-9,       # Molecular diffusion coefficient [m²/d]
        'D_m_solve': True,
        'D_m_min': 1e-12,
        'D_m_max': 0.001,

        'alpha_L': 0.2,        # Longitudinal dispersivity [m]
        'alpha_L_solve': True,
        'alpha_L_min': 0.01,
        'alpha_L_max': 1.0,

        # Boundary conditions
        'C_lake': 220,       # Lake concentration at top [M/L³]
        'C_gw': 136,          # Groundwater concentration at bottom [M/L³]
        
        # Plotting parameters
        'N': 200,              # Number of points for plotting (evaluated over plot_depth)

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
        print("Observed data file 'conc_data.xlsx' not found. Cannot perform calibration.")
        return None
    except Exception as exc:
        print(f"Unable to read observed data (DW25 A1:AL27): {exc}. Cannot perform calibration.")
        return None
    
    # Perform optimization
    print("\n" + "="*60)
    print("PARAMETER CALIBRATION")
    print("="*60)
    opt_result = optimize_parameters(params, observed_pd, verbose=True)
    
    # Extract optimized parameters
    params_optimal = opt_result['params_optimal']
    rmse_optimal = opt_result['rmse_optimal']
    optimal_values = opt_result['optimal_values']
    iteration_history = opt_result.get('iteration_history', [])
    sensitivities = opt_result.get('sensitivities', {})
    
    # Display iteration history table
    if iteration_history:
        print("\n" + "="*80)
        print("OPTIMIZATION ITERATION HISTORY")
        print("="*80)
        df_history = pd.DataFrame(iteration_history)
        # Format column names for better display
        df_display = df_history.copy()
        df_display.columns = [col.replace('_', ' ').title() if col != 'iteration' else 'Iter' 
                              for col in df_display.columns]
        # Format the display
        with pd.option_context('display.max_columns', None, 
                               'display.width', None,
                               'display.max_colwidth', None,
                               'display.float_format', lambda x: f'{x:.6g}'):
            print(df_display.to_string(index=False))
        print("="*80)
    
    # Compute steady-state solution with optimized parameters
    result = solve_ss(params_optimal, verbose=True)
    
    # Extract results
    C = result['C']
    z = result['z']
    params_final = result['params']
    
    # Print summary
    print("\n" + "="*60)
    print("CALIBRATED STEADY-STATE SOLUTION SUMMARY")
    print("="*60)
    print(f"Optimal RMSE: {rmse_optimal:.6f}")
    print(f"\nOptimal parameter values:")
    for name, value in optimal_values.items():
        print(f"  {name}: {value:.6g}")
    
    # Display parameter sensitivities
    if sensitivities:
        print("\n" + "="*60)
        print("PARAMETER SENSITIVITY ANALYSIS")
        print("="*60)
        print("Sensitivity: Change in RMSE per unit change in parameter")
        print("\nAbsolute sensitivities [RMSE/parameter unit]:")
        for name in optimal_values.keys():
            if name in sensitivities.get('absolute_sensitivities', {}):
                sens = sensitivities['absolute_sensitivities'][name]
                print(f"  {name:12s}: {sens:12.6g}")
        
        print("\nRelative sensitivities [% change in RMSE per 1% change in parameter]:")
        for name in optimal_values.keys():
            if name in sensitivities.get('relative_sensitivities', {}):
                sens = sensitivities['relative_sensitivities'][name]
                print(f"  {name:12s}: {sens:12.6g}")
        
        print("\nNormalized sensitivities (sum = 1.0, higher = more sensitive):")
        for name in optimal_values.keys():
            if name in sensitivities.get('normalized_sensitivities', {}):
                sens = sensitivities['normalized_sensitivities'][name]
                print(f"  {name:12s}: {sens:12.6f}")
        print("="*60)
    
    print(f"\nConcentration at top (z=0):     {C[0]:.2f}")
    print(f"Concentration at bottom (z=L): {C[-1]:.2f}")
    print(f"Concentration at midpoint:     {C[len(C)//2]:.2f}")
    print(f"\nDerived parameters:")
    print(f"  Effective dispersion (D_eff): {params_final['D_eff']:.6f} m²/d")
    print(f"  Péclet number (global):      {params_final['Pe_global']:.4f}")
    print("="*60)
    
    # Optional: Plot results (if matplotlib is available)
    try:
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(6, 10))
        
        # Determine plotting depth
        depth_limit = params_final.get('plot_depth', params_final.get('L', None))
        if depth_limit is None:
            depth_limit = z[-1]
        try:
            depth_limit = float(depth_limit)
        except Exception:
            depth_limit = z[-1]
        # Clamp to [0, L]
        L_val = params_final.get('L', z[-1])
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
        plt.plot(C_plot, z_plot, 'b-', linewidth=2.5, label='Calibrated steady-state')
        
        # Overlay observed points if available and specie column exists
        specie_val = params_final.get('specie')
        specie_label = str(specie_val).strip() if isinstance(specie_val, str) else None
        specie_title = specie_label.title() if specie_label else None
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
        plt.axvline(x=params_final['C_lake'], color='r', linestyle=':', 
                    linewidth=1.5, alpha=0.7, label='Lake boundary')
        plt.axvline(x=params_final['C_gw'], color='g', linestyle=':', 
                    linewidth=1.5, alpha=0.7, label='Groundwater boundary')
        
        plt.xlabel('Concentration [M/L³]')
        plt.ylabel('Depth z [m]')
        plot_title = 'Calibrated Steady-State Concentration Profile'
        if specie_title:
            plot_title += f' ({specie_title})'
        plt.title(plot_title)
        
        # Create box label with RMSE and optimal parameter values
        label_lines = [f'RMSE = {rmse_optimal:.3g}']
        label_lines.append('')
        label_lines.append('Optimal parameters:')
        for name, value in optimal_values.items():
            # Format parameter name nicely
            param_label = name.replace('_', ' ').title()
            label_lines.append(f'{param_label} = {value:.4g}')
        
        label_text = '\n'.join(label_lines)
        plt.text(
            0.98, 0.02, label_text,
            transform=plt.gca().transAxes,
            ha='right', va='bottom',
            fontsize=9,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray', pad=0.5)
        )
        
        plt.legend(loc='best', fontsize=9)
        plt.grid(True, alpha=0.3)
        plt.gca().invert_yaxis()  # Invert so z=0 is at top
        # Apply depth limit on y-axis
        plt.ylim(depth_limit, 0)
        plt.tight_layout()
        plt.savefig('concentration_profile_calibrated.png', dpi=150)
        plt.show()
        print("\nPlot saved to 'concentration_profile_calibrated.png'")
    except ImportError:
        print("\nMatplotlib not available. Skipping plot.")
    
    return {
        'result': result,
        'optimization': opt_result
    }


if __name__ == '__main__':
    result = main()

