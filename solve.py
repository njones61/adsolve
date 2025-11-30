"""
1D Advection-Diffusion Solver

Solves the 1D advection-diffusion equation in a vertical soil column:
    dC/dt = D_eff * d²C/dz² + v * dC/dz

with boundary conditions:
    C(0, t) = C_lake
    C(L, t) = C_gw
"""

import numpy as np
from scipy.linalg import solve_banded


def calculate_derived_parameters(params):
    """
    Calculate derived parameters from input parameters.
    
    Parameters
    ----------
    params : dict
        Dictionary containing input parameters
        
    Returns
    -------
    dict
        Dictionary with calculated parameters added
    """
    # Extract input parameters
    L = params['L']
    theta = params['porosity']
    K = params['K']
    delta_h = params['delta_h']
    D_m = params['D_m']
    alpha_L = params['alpha_L']
    N = params['N']
    
    # Calculate grid spacing
    delta_z = L / N
    
    # Calculate Darcy velocity
    q = K * delta_h / L
    
    # Calculate pore water velocity
    v = q / theta
    
    # Calculate mechanical dispersion
    D_mech = alpha_L * v
    
    # Calculate effective dispersion
    D_eff = D_m + D_mech
    
    # Calculate local Péclet number
    Pe_local = v * delta_z / D_eff
    
    # Calculate global Péclet number
    Pe_global = v * L / D_eff
    
    # Add calculated parameters to dictionary
    params['delta_z'] = delta_z
    params['q'] = q
    params['v'] = v
    params['D_mech'] = D_mech
    params['D_eff'] = D_eff
    params['Pe_local'] = Pe_local
    params['Pe_global'] = Pe_global
    
    return params


def initialize_concentration(params):
    """
    Initialize the concentration profile (default, simple choice).
    
    Parameters
    ----------
    params : dict
        Dictionary containing parameters including boundary concentrations
        
    Returns
    -------
    numpy.ndarray
        Initial concentration array of size N+1
    """
    N = params['N']
    C_gw = params['C_gw']
    C_lake = params['C_lake']
    L = params['L']
    delta_z = params['delta_z']
    
    # Create spatial grid
    z = np.linspace(0, L, N + 1)
    
    # Simple, robust default: start uniform at groundwater concentration
    C = np.full(N + 1, C_gw, dtype=float)
    
    # Apply boundary conditions at t=0
    C[0] = C_lake
    C[-1] = C_gw
    
    return C


def explicit_step(C, params):
    """
    Perform one explicit (Forward Euler) time step.
    
    Parameters
    ----------
    C : numpy.ndarray
        Current concentration profile
    params : dict
        Dictionary containing parameters
        
    Returns
    -------
    numpy.ndarray
        Updated concentration profile
    """
    N = params['N']
    delta_t = params['delta_t']
    delta_z = params['delta_z']
    D_eff = params['D_eff']
    v = params['v']
    C_lake = params['C_lake']
    C_gw = params['C_gw']
    
    # Create new concentration array
    C_new = C.copy()
    
    # Update interior points
    for i in range(1, N):
        # Diffusion term (central difference)
        diffusion = D_eff * (C[i+1] - 2*C[i] + C[i-1]) / (delta_z**2)
        
        # Advection term (upwind difference for v > 0)
        advection = v * (C[i+1] - C[i]) / delta_z
        
        # Update concentration
        C_new[i] = C[i] + delta_t * (diffusion + advection)
    
    # Apply boundary conditions
    C_new[0] = C_lake
    C_new[-1] = C_gw
    
    return C_new


def implicit_step(C, params):
    """
    Perform one implicit (Backward Euler) time step using tridiagonal solver.
    
    Parameters
    ----------
    C : numpy.ndarray
        Current concentration profile
    params : dict
        Dictionary containing parameters
        
    Returns
    -------
    numpy.ndarray
        Updated concentration profile
    """
    N = params['N']
    delta_t = params['delta_t']
    delta_z = params['delta_z']
    D_eff = params['D_eff']
    v = params['v']
    C_lake = params['C_lake']
    C_gw = params['C_gw']
    
    # Set up tridiagonal system: A * C_new = b
    # For interior points i = 1, 2, ..., N-1:
    # -a_i * C_{i-1} + b_i * C_i - c_i * C_{i+1} = C_i^n
    
    # Coefficients
    a = np.zeros(N - 1)  # Lower diagonal (not used for first point)
    b = np.zeros(N - 1)  # Main diagonal
    c = np.zeros(N - 1)  # Upper diagonal (not used for last point)
    rhs = np.zeros(N - 1)  # Right-hand side
    
    # Fill coefficients for interior points
    for i in range(1, N):
        idx = i - 1  # Index in the system (0 to N-2)
        
        # Lower diagonal (coefficient of C_{i-1})
        if i > 1:
            a[idx] = -delta_t * D_eff / (delta_z**2)
        
        # Main diagonal (coefficient of C_i)
        b[idx] = 1 + 2*delta_t*D_eff/(delta_z**2) + delta_t*v/delta_z
        
        # Upper diagonal (coefficient of C_{i+1})
        if i < N - 1:
            c[idx] = -delta_t*D_eff/(delta_z**2) - delta_t*v/delta_z
        
        # Right-hand side
        rhs[idx] = C[i]
        
        # Account for boundary conditions in RHS
        if i == 1:
            # C[0] = C_lake affects first equation
            rhs[idx] += delta_t*D_eff/(delta_z**2) * C_lake
        if i == N - 1:
            # C[N] = C_gw affects last equation
            rhs[idx] += (delta_t*D_eff/(delta_z**2) + delta_t*v/delta_z) * C_gw
    
    # Solve tridiagonal system
    # Format for solve_banded: [upper, main, lower] diagonals
    ab = np.zeros((3, N - 1))
    ab[0, 1:] = c[:-1]  # Upper diagonal (shifted)
    ab[1, :] = b        # Main diagonal
    ab[2, :-1] = a[1:]  # Lower diagonal (shifted)
    
    # Solve
    C_interior = solve_banded((1, 1), ab, rhs)
    
    # Construct full solution
    C_new = np.zeros(N + 1)
    C_new[0] = C_lake
    C_new[1:N] = C_interior
    C_new[N] = C_gw
    
    return C_new


def crank_nicolson_step(C, params):
    """
    Perform one Crank-Nicolson time step.
    
    Parameters
    ----------
    C : numpy.ndarray
        Current concentration profile
    params : dict
        Dictionary containing parameters
        
    Returns
    -------
    numpy.ndarray
        Updated concentration profile
    """
    N = params['N']
    delta_t = params['delta_t']
    delta_z = params['delta_z']
    D_eff = params['D_eff']
    v = params['v']
    C_lake = params['C_lake']
    C_gw = params['C_gw']
    
    # Calculate explicit part (at time n)
    explicit_part = np.zeros(N - 1)
    for i in range(1, N):
        idx = i - 1
        # Diffusion term
        diffusion = D_eff * (C[i+1] - 2*C[i] + C[i-1]) / (delta_z**2)
        # Advection term
        advection = v * (C[i+1] - C[i]) / delta_z
        explicit_part[idx] = C[i] + 0.5 * delta_t * (diffusion + advection)
    
    # Set up implicit part (at time n+1)
    # Tridiagonal system similar to implicit method
    a = np.zeros(N - 1)
    b = np.zeros(N - 1)
    c = np.zeros(N - 1)
    rhs = np.zeros(N - 1)
    
    for i in range(1, N):
        idx = i - 1
        
        # Lower diagonal
        if i > 1:
            a[idx] = -0.5 * delta_t * D_eff / (delta_z**2)
        
        # Main diagonal
        b[idx] = 1 + delta_t*D_eff/(delta_z**2) + 0.5*delta_t*v/delta_z
        
        # Upper diagonal
        if i < N - 1:
            c[idx] = -0.5*delta_t*D_eff/(delta_z**2) - 0.5*delta_t*v/delta_z
        
        # Right-hand side (explicit part)
        rhs[idx] = explicit_part[idx]
        
        # Account for boundary conditions
        if i == 1:
            rhs[idx] += 0.5*delta_t*D_eff/(delta_z**2) * C_lake
        if i == N - 1:
            rhs[idx] += (0.5*delta_t*D_eff/(delta_z**2) + 0.5*delta_t*v/delta_z) * C_gw
    
    # Solve tridiagonal system
    ab = np.zeros((3, N - 1))
    ab[0, 1:] = c[:-1]
    ab[1, :] = b
    ab[2, :-1] = a[1:]
    
    C_interior = solve_banded((1, 1), ab, rhs)
    
    # Construct full solution
    C_new = np.zeros(N + 1)
    C_new[0] = C_lake
    C_new[1:N] = C_interior
    C_new[N] = C_gw
    
    return C_new


def check_stability(params):
    """
    Check stability condition for explicit method.
    
    Parameters
    ----------
    params : dict
        Dictionary containing parameters
        
    Returns
    -------
    bool
        True if stable, False otherwise
    float
        Maximum stable time step
    """
    delta_z = params['delta_z']
    D_eff = params['D_eff']
    v = params['v']
    delta_t = params['delta_t']
    
    # Stability condition: delta_t <= delta_z^2 / (2*D_eff + v*delta_z)
    max_delta_t = (delta_z**2) / (2*D_eff + v*delta_z)
    
    is_stable = delta_t <= max_delta_t
    
    return is_stable, max_delta_t


def solve_trans(params, method='crank_nicolson', verbose=True):
    """
    Solve the 1D advection-diffusion problem (transient).
    
    Parameters
    ----------
    params : dict
        Dictionary containing all input parameters
    method : str, optional
        Time stepping method: 'explicit', 'implicit', or 'crank_nicolson' (default)
    verbose : bool, optional
        Print progress information (default: True)
        
    Returns
    -------
    dict
        Dictionary containing:
        - 'C': concentration array at final time [N+1]
        - 'C_history': concentration history [n_steps+1, N+1] (if save_history=True)
        - 't': time array [n_steps+1]
        - 'z': spatial grid [N+1]
        - 'params': updated parameters dictionary
    """
    # Calculate derived parameters
    params = calculate_derived_parameters(params)
    
    # Extract parameters
    N = params['N']
    delta_t = params['delta_t']
    t_max = params['t_max']
    L = params['L']
    
    # Create time array
    n_steps = int(t_max / delta_t)
    t = np.linspace(0, t_max, n_steps + 1)
    
    # Create spatial grid
    z = np.linspace(0, L, N + 1)
    
    # Initialize concentration
    C = initialize_concentration(params)
    
    # Check stability for explicit method
    if method == 'explicit':
        is_stable, max_delta_t = check_stability(params)
        if not is_stable:
            if verbose:
                print(f"WARNING: Time step {delta_t:.6f} may be unstable.")
                print(f"Maximum stable time step: {max_delta_t:.6f}")
    
    # Storage for history (optional)
    save_history = params.get('save_history', False)
    if save_history:
        C_history = np.zeros((n_steps + 1, N + 1))
        C_history[0, :] = C.copy()
    
    # Storage for output snapshots at intervals
    output_interval = params.get('output_interval', None)
    C_snapshots = []
    t_snapshots = []
    next_output_time = 0.0
    
    # Save initial condition
    if output_interval is not None:
        C_snapshots.append(C.copy())
        t_snapshots.append(0.0)
        next_output_time = output_interval
    
    # Time stepping
    if verbose:
        print(f"Solving with {method} method...")
        print(f"Time steps: {n_steps}, Grid points: {N+1}")
        print(f"Péclet number (local): {params['Pe_local']:.4f}")
        print(f"Péclet number (global): {params['Pe_global']:.4f}")
        if output_interval is not None:
            print(f"Output interval: {output_interval} days")
    
    # Select time stepping function
    if method == 'explicit':
        step_func = explicit_step
    elif method == 'implicit':
        step_func = implicit_step
    elif method == 'crank_nicolson':
        step_func = crank_nicolson_step
    else:
        raise ValueError(f"Unknown method: {method}. Use 'explicit', 'implicit', or 'crank_nicolson'")
    
    # Time loop
    for n in range(n_steps):
        C = step_func(C, params)
        current_time = t[n + 1]
        
        if save_history:
            C_history[n + 1, :] = C.copy()
        
        # Save snapshot at output intervals
        if output_interval is not None:
            # Check if we should save at this time step
            # Save if we've reached or passed the next output time
            # or if it's the final time step
            save_snapshot = False
            
            if n == n_steps - 1:  # Always save final state
                # Only save if not already saved at this time
                if len(t_snapshots) == 0 or abs(t_snapshots[-1] - current_time) > delta_t / 2:
                    save_snapshot = True
            elif current_time >= next_output_time - delta_t / 2:
                # Save snapshot when we reach the output interval
                save_snapshot = True
                # Update next output time
                next_output_time += output_interval
            
            if save_snapshot:
                C_snapshots.append(C.copy())
                t_snapshots.append(current_time)
        
        if verbose and (n + 1) % (n_steps // 10) == 0:
            print(f"  Step {n+1}/{n_steps} (t = {current_time:.2f})")
    
    # Prepare output
    result = {
        'C': C,
        't': t,
        'z': z,
        'params': params
    }
    
    if save_history:
        result['C_history'] = C_history
    
    if output_interval is not None:
        result['C_snapshots'] = C_snapshots
        result['t_snapshots'] = t_snapshots
    
    if verbose:
        print("Solution complete!")
    
    return result


def solve_ss(params, verbose=True):
    """
    Compute the steady-state solution to the 1D advection-diffusion equation
    with Dirichlet boundary conditions using the analytical expression.
    
    Parameters
    ----------
    params : dict
        Dictionary containing all input parameters. Time-related parameters
        (e.g., 'delta_t', 't_max') are ignored.
    verbose : bool, optional
        Print summary information (default: True)
    
    Returns
    -------
    dict
        Dictionary containing:
        - 'C': steady concentration array [N+1]
        - 'z': spatial grid [N+1]
        - 'params': updated parameters dictionary (with derived values)
    """
    # Calculate derived parameters
    params = calculate_derived_parameters(params)
    
    # Extract needed parameters
    N = params['N']
    L = params['L']
    D_eff = params['D_eff']
    v = params['v']
    C_lake = params['C_lake']
    C_gw = params['C_gw']
    
    # Spatial grid
    z = np.linspace(0, L, N + 1)
    
    # Global Péclet number
    Pe = v * L / D_eff if D_eff != 0 else np.inf
    params['Pe_global'] = Pe
    
    # Handle small Pe with linear limit to avoid 0/0
    if np.isfinite(Pe) and abs(Pe) < 1e-10:
        # Linear interpolation between boundaries
        C = C_lake + (C_gw - C_lake) * (z / L)
    else:
        # General analytical solution
        # Avoid numerical issues if denominator is tiny
        denom = 1.0 - np.exp(-Pe)
        # If denom is too small, fall back to linear
        if not np.isfinite(denom) or abs(denom) < 1e-14:
            C = C_lake + (C_gw - C_lake) * (z / L)
        else:
            C = C_lake + (C_gw - C_lake) * (1.0 - np.exp(-Pe * z / L)) / denom
    
    # Enforce boundary values exactly
    C[0] = C_lake
    C[-1] = C_gw
    
    if verbose:
        print("Steady-state solution computed.")
        print(f"Grid points: {N+1}")
        print(f"Péclet number (global): {params['Pe_global']:.6f}")
    
    return {
        'C': C,
        'z': z,
        'params': params
    }

