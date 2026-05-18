"""Numerical methods for the filleted infinite potential well.

Implements:
  1. Shooting method using rifusaki-catkit ODE integrators + root-finders.
  2. Finite Difference Method (FDM) using scipy.linalg.eigh_tridiagonal.

All quantities are dimensionless: energies in units of E₀ = ℏ²/(2m a₀²),
lengths in units of a₀.
"""

import numpy as np
from numpy.typing import NDArray
from typing import Tuple, Optional
import scipy.linalg as la
import rifusaki_catkit.integrators as integ
import rifusaki_catkit.roots as rts

from .potential import v_tilde


# ============================================================
# Shooting Method
# ============================================================

def schrodinger_ode(
    x_tilde: float,
    y: NDArray[np.floating],
    E_tilde: float,
    alpha: float,
    beta: float,
    gamma: float,
) -> NDArray[np.floating]:
    """RHS of the 1D Schrödinger equation as a first-order ODE system.

    y = [ψ, ψ']
    y' = [ψ', (Ṽ - Ẽ) ψ]

    In dimensionless form:  ψ'' = (Ṽ(x̃) - Ẽ) ψ
    """
    psi, psi_prime = y[0], y[1]
    V = v_tilde(x_tilde, alpha, beta, gamma)
    if np.isinf(V):
        # At or beyond the wall, force ψ to decay
        return np.array([psi_prime, 1e10 * psi])
    return np.array([psi_prime, (V - E_tilde) * psi])


def shoot(
    E_tilde: float,
    alpha: float,
    beta: float,
    gamma: float,
    parity: str = "even",
    h: float = 0.001,
) -> Tuple[NDArray[np.floating], NDArray[np.floating]]:
    """Integrate the Schrödinger equation from x̃=0 to x̃=L.

    Parameters
    ----------
    E_tilde : float
        Trial dimensionless energy.
    alpha, beta, gamma : float
        Well geometry parameters.
    parity : str
        "even" or "odd".
    h : float
        Step size for the ODE integrator.

    Returns
    -------
    t : ndarray
        x̃ positions.
    y : ndarray
        [ψ(x̃), ψ'(x̃)] at each position.
    """
    L = alpha + gamma

    # Initial conditions at x̃ = 0
    if parity == "even":
        y0 = np.array([1.0, 0.0])   # ψ(0) = 1, ψ'(0) = 0
    else:
        y0 = np.array([0.0, 1.0])   # ψ(0) = 0, ψ'(0) = 1

    def f(t: float, y: NDArray[np.floating]) -> NDArray[np.floating]:
        return schrodinger_ode(t, y, E_tilde, alpha, beta, gamma)

    t, y = integ.integrate(f, (0.0, L), y0, h=h, step_fn=integ.rk4_step)
    return t, y


def boundary_value(
    E_tilde: float,
    alpha: float,
    beta: float,
    gamma: float,
    parity: str = "even",
    h: float = 0.001,
) -> float:
    """Return ψ(L; Ẽ) for a given trial energy (the miss-distance).

    Used as the root-finding target: find Ẽ such that boundary_value(Ẽ) = 0.
    """
    _, y = shoot(E_tilde, alpha, beta, gamma, parity, h)
    return y[-1, 0]  # ψ at x̃ = L


def find_eigenvalue_shooting(
    alpha: float,
    beta: float,
    gamma: float,
    parity: str = "even",
    E_low: Optional[float] = None,
    E_high: Optional[float] = None,
    h: float = 0.001,
    tol: float = 1e-8,
    max_iter: int = 100,
) -> float:
    """Find a bound-state eigenvalue Ẽ using the shooting method.

    Uses catkit's bisection root-finding. The bracket [E_low, E_high] must
    contain exactly one root (sign change in ψ(L; Ẽ)).

    Parameters
    ----------
    alpha, beta, gamma : float
        Well geometry.
    parity : str
        "even" or "odd".
    E_low, E_high : float or None
        Energy bracket. If None, a bracket is searched automatically.
    h : float
        ODE step size.
    tol : float
        Convergence tolerance.
    max_iter : int
        Maximum iterations.

    Returns
    -------
    float
        Dimensionless energy eigenvalue.
    """
    L = alpha + gamma

    def f(E: float) -> float:
        return boundary_value(E, alpha, beta, gamma, parity, h)

    if E_low is not None and E_high is not None:
        # Use provided bracket
        return rts.find_root(
            f, df=None,
            state0=(E_low, E_high),
            step_fn=rts.bisection_step,
            tol=tol, max_iter=max_iter,
        )

    # Automatic bracket search
    # Start from just above the well bottom
    E_a = -beta + 0.01
    E_b = -beta + 0.5
    # Expand upward until sign change or hit E=0
    fa = f(E_a)
    for _ in range(50):
        fb = f(E_b)
        if fa * fb <= 0:
            break
        fa = fb
        E_a = E_b
        E_b += 0.5
        if E_b > -0.01:  # don't go above E=0 (bound states only)
            E_b = -0.01
            break
    else:
        raise RuntimeError(f"Could not bracket eigenvalue for {parity} parity")

    return rts.find_root(
        f, df=None,
        state0=(E_a, E_b),
        step_fn=rts.bisection_step,
        tol=tol, max_iter=max_iter,
    )


def find_all_eigenvalues_shooting(
    alpha: float,
    beta: float,
    gamma: float,
    n_states: int = 4,
    h: float = 0.001,
) -> list[dict]:
    """Find the lowest n_states eigenvalues using the shooting method.

    Returns a list of dicts with keys: n, parity, E.
    """
    results = []
    parities = []
    for n in range(1, n_states + 1):
        parities.append("even" if n % 2 == 1 else "odd")

    # Bracket for each state
    E_prev = -beta  # start from well bottom
    for n, parity in enumerate(parities, 1):
        # Search for the n-th root
        # For ground state, search from -beta to something higher
        # For excited states, search above previous eigenvalue

        def f(E: float) -> float:
            return boundary_value(E, alpha, beta, gamma, parity, h)

        # Scan upward to find a bracket
        E_a = E_prev + 0.001
        E_b = E_a + 0.2
        fa = f(E_a)
        found = False
        for _ in range(100):
            fb = f(E_b)
            if fa * fb <= 0:
                found = True
                break
            fa = fb
            E_a = E_b
            E_b += 0.2
            if E_b > -0.001:
                break
        if not found:
            # Try finer scan
            E_b = E_a + 0.02
            for _ in range(200):
                fb = f(E_b)
                if fa * fb <= 0:
                    found = True
                    break
                fa = fb
                E_a = E_b
                E_b += 0.02
                if E_b > -0.001:
                    break

        if not found:
            raise RuntimeError(f"Could not bracket eigenvalue n={n} ({parity})")

        E_val = rts.find_root(
            f, df=None,
            state0=(E_a, E_b),
            step_fn=rts.bisection_step,
            tol=1e-8, max_iter=100,
        )
        results.append({"n": n, "parity": parity, "E": E_val})
        E_prev = E_val  # next state is above this one

    return results


# ============================================================
# Finite Difference Method (FDM)
# ============================================================

def build_hamiltonian_fdm(
    alpha: float,
    beta: float,
    gamma: float,
    N: int = 2000,
    parity: str = "even",
) -> Tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]:
    """Build the tridiagonal FDM Hamiltonian matrix for eigh_tridiagonal.

    The domain is x̃ ∈ [0, L] with L = α + γ.
    Boundary conditions:
      - Even parity: ψ'(0) = 0, ψ(L) = 0
      - Odd parity:  ψ(0) = 0, ψ(L) = 0

    Uses central differences: ψ'' ≈ (ψ_{i+1} - 2ψ_i + ψ_{i-1}) / h²

    For even parity, the ghost-point Neumann BC at x=0 makes the raw
    tridiagonal matrix non-symmetric (H_{0,1} ≠ H_{1,0}).  We apply a
    similarity transform D H D⁻¹ with D = diag(1, √2, √2, …) before
    passing the matrix to eigh_tridiagonal.  This restores symmetry
    without changing eigenvalues.  The eigenvectors must be unscaled
    back (ψ_i = v_i / √2 for i ≥ 1) after the eigensolve.

    Parameters
    ----------
    alpha, beta, gamma : float
        Well geometry.
    N : int
        Number of interior grid points on [0, L].
    parity : str
        "even" or "odd".

    Returns
    -------
    diagonal : ndarray
        Main diagonal of the Hamiltonian (after similarity transform).
    off_diagonal : ndarray
        Off-diagonal (sub-diagonal / super-diagonal, after transform).
    x_grid : ndarray
        Spatial grid points.
    """
    L = alpha + gamma
    h = L / N  # grid spacing

    if parity == "even":
        # Interior points: i = 0, 1, ..., N-1  (ψ_N = 0 at boundary)
        M = N
        i_start = 0
    else:
        # Interior points: i = 1, 2, ..., N-1  (ψ_0 = 0, ψ_N = 0)
        M = N - 1
        i_start = 1

    x_grid = np.array([i * h for i in range(i_start, i_start + M)])
    V_grid = v_tilde(x_grid, alpha, beta, gamma)

    # Standard interior: diag = 2/h² + V_i, off = -1/h²
    diagonal = np.full(M, 2.0 / h**2) + V_grid
    off_diagonal = np.full(M - 1, -1.0 / h**2)

    if parity == "even":
        # Neumann BC at x=0: ψ'(0) = 0 → ψ_{-1} = ψ₁
        # ψ''₀ = (ψ₁ - 2ψ₀ + ψ₁)/h² = 2(ψ₁ - ψ₀)/h²
        # Raw stencil gives H_{0,1} = -2/h², which is non-symmetric with
        # the standard interior off-diagonal (-1/h²), breaking eigh_tridiagonal.
        # Apply similarity transform D H D⁻¹ with D = diag(1, √2, √2, …):
        #   H'_{0,1} = 1 · (-2/h²) / √2 = -√2/h²  ← symmetric with H'_{1,0}
        # Diagonal elements are unchanged by the transform.
        # Eigenvectors are unscaled back after eigh_tridiagonal.
        diagonal[0] = 2.0 / h**2 + V_grid[0]
        off_diagonal[0] = -np.sqrt(2.0) / h**2

    return diagonal, off_diagonal, x_grid


def solve_fdm(
    alpha: float,
    beta: float,
    gamma: float,
    N: int = 2000,
    n_states: int = 6,
) -> list[dict]:
    """Solve the eigenvalue problem using the FDM.

    Returns a list of dicts with keys: n, parity, E, psi (eigenvector), x_grid.
    """
    results = []

    for parity in ["even", "odd"]:
        diag, off_diag, x_grid = build_hamiltonian_fdm(
            alpha, beta, gamma, N, parity,
        )
        # Solve tridiagonal eigenvalue problem
        eigenvalues, eigenvectors = la.eigh_tridiagonal(diag, off_diag)

        # Unscale eigenvectors from the similarity transform (even parity only)
        if parity == "even":
            for j in range(eigenvectors.shape[1]):
                eigenvectors[1:, j] /= np.sqrt(2.0)

        # Take the lowest eigenvalues
        n_eigs = min(n_states // 2 + 1, len(eigenvalues))
        for j in range(n_eigs):
            E_val = eigenvalues[j]
            psi = eigenvectors[:, j]

            # Only consider bound states (E < 0)
            if E_val >= 0.0:
                continue

            # Normalize the eigenvector on the full domain [-L, L]
            # For even parity, ψ extends symmetrically; for odd, antisymmetrically.
            psi_normalized = _normalize_eigenvector(psi, x_grid, parity)

            # Determine quantum number
            if parity == "even":
                n = 2 * j + 1  # n=1, 3, 5, ...
            else:
                n = 2 * j + 2  # n=2, 4, 6, ...

            results.append({
                "n": n,
                "parity": parity,
                "E": float(E_val),
                "psi": psi_normalized,
                "x_grid": x_grid,
            })

    # Sort by energy
    results.sort(key=lambda r: r["E"])
    # Re-number sequentially
    for i, r in enumerate(results):
        r["n"] = i + 1

    return results[:n_states]


def _normalize_eigenvector(
    psi_half: NDArray[np.floating],
    x_grid: NDArray[np.floating],
    parity: str,
) -> NDArray[np.floating]:
    """Normalise the half-domain eigenvector on the full domain [-L, L].

    Uses Simpson's rule for integration.
    """
    h = x_grid[1] - x_grid[0]
    M = len(psi_half)

    # Build full-domain psi (include x=0 and x=L for integration)
    if parity == "even":
        # ψ(-x) = ψ(x), ψ(0) is included
        # Full grid: x_{-M+1}, ..., x_{-1}, x_0, x_1, ..., x_{M-1}
        # But we have values at x_0 (i=0), ..., x_{M-1}
        # Need to integrate from -L to L
        # Symmetry: ∫₋Lᴸ ψ² dx = 2 ∫₀ᴸ ψ² dx
        # Simpson on [0, L] with M points (including x=0)
        y_sq = psi_half**2

        # Simpson integration on the half-domain
        if M >= 3:
            norm_sq = 2.0 * (h / 3.0) * (
                y_sq[0] + y_sq[-1] +
                4.0 * np.sum(y_sq[1:-1:2]) +
                2.0 * np.sum(y_sq[2:-2:2])
            )
        else:
            norm_sq = 2.0 * h * np.sum(y_sq)  # fallback
    else:
        # ψ(0) = 0, so first point is x = h
        # Full grid integration with ψ(-x) = -ψ(x) → ψ²(-x) = ψ²(x)
        y_sq = np.concatenate([[0.0], psi_half**2])
        x_full = np.concatenate([[0.0], x_grid])
        # Simpson on [0, L]
        M_full = len(y_sq)
        if M_full >= 3:
            norm_sq = 2.0 * (h / 3.0) * (
                y_sq[0] + y_sq[-1] +
                4.0 * np.sum(y_sq[1:-1:2]) +
                2.0 * np.sum(y_sq[2:-2:2])
            )
        else:
            norm_sq = 2.0 * h * np.sum(y_sq)

    norm_factor = 1.0 / max(np.sqrt(norm_sq), 1e-30)
    return psi_half * norm_factor
