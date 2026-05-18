"""Analytical approximate methods for the filleted infinite potential well.

Implements:
  1. First-order perturbation theory using the shifted infinite square well
     as the unperturbed Hamiltonian.
  2. Variational method with parameterized trial wavefunctions for the
     ground state (even parity) and first excited state (odd parity).

All quantities are dimensionless: energies in units of E₀ = ℏ²/(2m a₀²),
lengths in units of a₀.
"""

import numpy as np
from numpy.typing import NDArray
from typing import Callable, Tuple
import rifusaki_catkit.roots as rts

from .potential import v_tilde


# ============================================================
# Numerical quadrature (Simpson's rule)
# ============================================================

def simpson(
    f: Callable[[float], float],
    a: float,
    b: float,
    n: int = 10001,
) -> float:
    """Composite Simpson's rule integration over [a, b].

    Uses an odd number of points (n must be odd; if even, adds 1).
    """
    if n % 2 == 0:
        n += 1
    x = np.linspace(a, b, n)
    h = (b - a) / (n - 1)
    y = np.array([f(xi) for xi in x])
    return (h / 3.0) * (y[0] + y[-1] + 4.0 * np.sum(y[1:-1:2]) + 2.0 * np.sum(y[2:-2:2]))


# ============================================================
# Perturbation Theory
# ============================================================

def infinite_well_wavefunction(
    x_tilde: float,
    n: int,
    L: float,
) -> float:
    """Unperturbed wavefunction for infinite square well [-L, L].

    ψ_n^{(0)}(x̃) = √(1/L) sin(nπ(x̃+L) / (2L))

    Normalised: ∫₋Lᴸ |ψ|² dx̃ = 1.
    Parity: n odd → even, n even → odd.
    """
    return np.sqrt(1.0 / L) * np.sin(n * np.pi * (x_tilde + L) / (2.0 * L))


def infinite_well_energy(n: int, L: float, beta: float = 0.0) -> float:
    """Dimensionless unperturbed energy for shifted infinite well.

    E_n^{(0)} = n²π² / (8 L²) - β

    The well has constant potential -β inside [-L, L] and ∞ outside.
    """
    return (n * np.pi) ** 2 / (8.0 * L ** 2) - beta


def first_order_correction(
    n: int,
    alpha: float,
    beta: float,
    gamma: float,
    n_quad: int = 10001,
) -> float:
    """First-order perturbation energy correction E_n^{(1)}.

    Uses the shifted infinite well as H₀ (V = -β inside [-L, L], ∞ outside).
    The perturbation is ΔṼ = Ṽ_filleted(x̃) - (-β) for |x̃| < L.

    Parameters
    ----------
    n : int
        Quantum number (1, 2, 3, ...).
    alpha, beta, gamma : float
        Well geometry parameters.
    n_quad : int
        Number of quadrature points.

    Returns
    -------
    float
        E_n^{(1)} = ⟨ψ_n^{(0)}|ΔṼ|ψ_n^{(0)}⟩
    """
    L = alpha + gamma

    def integrand(x_tilde: float) -> float:
        psi = infinite_well_wavefunction(x_tilde, n, L)
        v = v_tilde(x_tilde, alpha, beta, gamma)
        # Perturbation: ΔṼ = Ṽ_filleted - (-β) = Ṽ_filleted + β
        delta_v = v + beta
        return psi**2 * delta_v

    # Integrate over [0, L] and double (symmetry)
    half = simpson(integrand, 0.0, L, n_quad)
    return 2.0 * half


def perturbation_energies(
    alpha: float,
    beta: float,
    gamma: float,
    n_states: int = 6,
    n_quad: int = 10001,
) -> list[dict]:
    """Compute perturbation theory energies for the lowest n_states.

    Returns a list of dicts with keys: n, parity, E0, E1, E_total.
    """
    L = alpha + gamma
    results = []
    for n in range(1, n_states + 1):
        E0 = infinite_well_energy(n, L, beta)
        E1 = first_order_correction(n, alpha, beta, gamma, n_quad)
        parity = "even" if n % 2 == 1 else "odd"
        results.append({
            "n": n,
            "parity": parity,
            "E0": E0,
            "E1": E1,
            "E_total": E0 + E1,
        })
    return results


# ============================================================
# Variational Method
# ============================================================

def trial_ground(
    x_tilde: float,
    lam: float,
    L: float,
) -> float:
    """Trial wavefunction for the ground state (even parity).

    ψ(x̃; λ) = (1 - x̃²/L²) · exp(-λ x̃² / (2L²))

    This satisfies ψ(L) = 0 and ψ'(0) = 0 automatically.
    """
    x = abs(x_tilde)
    if x >= L:
        return 0.0
    return (1.0 - x**2 / L**2) * np.exp(-lam * x**2 / (2.0 * L**2))


def trial_ground_derivative(
    x_tilde: float,
    lam: float,
    L: float,
) -> float:
    """Derivative dψ/dx̃ of the ground-state trial function."""
    x = abs(x_tilde)
    if x >= L:
        return 0.0
    pref = (1.0 - x**2 / L**2)
    exp_part = np.exp(-lam * x**2 / (2.0 * L**2))
    dpref = -2.0 * x / L**2
    dexp = -lam * x / L**2 * exp_part
    deriv = dpref * exp_part + pref * dexp
    return deriv * np.sign(x_tilde)  # restore sign for x < 0


def trial_excited(
    x_tilde: float,
    lam: float,
    L: float,
) -> float:
    """Trial wavefunction for the first excited state (odd parity).

    ψ(x̃; λ) = x̃ · (1 - x̃²/L²) · exp(-λ x̃² / (2L²))

    This satisfies ψ(0) = 0 and ψ(L) = 0 automatically.
    """
    x = abs(x_tilde)
    if x >= L:
        return 0.0
    return x_tilde * (1.0 - x**2 / L**2) * np.exp(-lam * x**2 / (2.0 * L**2))


def trial_excited_derivative(
    x_tilde: float,
    lam: float,
    L: float,
) -> float:
    """Derivative dψ/dx̃ of the excited-state trial function."""
    x = abs(x_tilde)
    if x >= L:
        return 0.0
    sign = np.sign(x_tilde)
    pref = (1.0 - x**2 / L**2)
    exp_part = np.exp(-lam * x**2 / (2.0 * L**2))
    # ψ = x̃ * pref * exp
    # ψ' = sign*pref*exp + x̃*(-2x/L²)*exp + x̃*pref*(-lam*x/L²)*exp
    term1 = sign * pref * exp_part
    term2 = x_tilde * (-2.0 * x / L**2) * exp_part
    term3 = x_tilde * pref * (-lam * x / L**2) * exp_part
    return term1 + term2 + term3


def normalization_integral(
    fn: Callable[[float], float],
    lam: float,
    L: float,
    n_quad: int = 10001,
) -> float:
    """Compute ∫₋Lᴸ |ψ|² dx̃."""
    def integrand(x_tilde: float) -> float:
        return fn(x_tilde) ** 2
    return 2.0 * simpson(integrand, 0.0, L, n_quad)


def kinetic_integral(
    fn: Callable[[float], float],
    fn_deriv: Callable[[float], float],
    lam: float,
    L: float,
    n_quad: int = 10001,
) -> float:
    """Compute ∫₋Lᴸ |ψ'|² dx̃ (kinetic energy numerator)."""
    def integrand(x_tilde: float) -> float:
        return fn_deriv(x_tilde) ** 2
    return 2.0 * simpson(integrand, 0.0, L, n_quad)


def potential_integral(
    fn: Callable[[float], float],
    lam: float,
    L: float,
    alpha: float,
    beta: float,
    gamma: float,
    n_quad: int = 10001,
) -> float:
    """Compute ∫₋Lᴸ |ψ|² Ṽ(x̃) dx̃."""
    def integrand(x_tilde: float) -> float:
        return fn(x_tilde) ** 2 * v_tilde(x_tilde, alpha, beta, gamma)
    return 2.0 * simpson(integrand, 0.0, L, n_quad)


def variational_energy(
    lam: float,
    L: float,
    alpha: float,
    beta: float,
    gamma: float,
    trial_fn: Callable[[float, float, float], float],
    trial_deriv: Callable[[float, float, float], float],
    n_quad: int = 10001,
) -> float:
    """Compute the energy functional Ẽ[λ] for a given trial function.

    Ẽ[λ] = (⟨T⟩ + ⟨V⟩) / ⟨ψ|ψ⟩
    where ⟨T⟩ = ∫ |ψ'|² dx̃, ⟨V⟩ = ∫ |ψ|² Ṽ dx̃.
    """
    # Pre-bind lambda and L
    psi = lambda x: trial_fn(x, lam, L)
    psi_p = lambda x: trial_deriv(x, lam, L)

    norm = normalization_integral(psi, lam, L, n_quad)
    T = kinetic_integral(psi, psi_p, lam, L, n_quad)
    V = potential_integral(psi, lam, L, alpha, beta, gamma, n_quad)

    return (T + V) / norm


def optimize_variational(
    L: float,
    alpha: float,
    beta: float,
    gamma: float,
    trial_fn: Callable[[float, float, float], float],
    trial_deriv: Callable[[float, float, float], float],
    lam_initial: float = 0.5,
    n_quad: int = 10001,
    tol: float = 1e-6,
    max_iter: int = 50,
) -> dict:
    """Minimize Ẽ[λ] with respect to λ using catkit's secant root-finding.

    Finds dẼ/dλ = 0 via finite differences.

    Returns dict with keys: lam_opt, E_opt, iterations.
    """
    def energy(lam: float) -> float:
        return variational_energy(
            lam, L, alpha, beta, gamma,
            trial_fn, trial_deriv, n_quad,
        )

    def denergy(lam: float) -> float:
        """Central difference derivative dE/dλ."""
        h = 1e-6 * max(abs(lam), 1.0)
        return (energy(lam + h) - energy(lam - h)) / (2.0 * h)

    # Try secant method first
    try:
        lam0 = lam_initial
        lam1 = lam_initial * 1.1
        lam_opt = rts.find_root(
            denergy,
            df=None,
            state0=(lam0, lam1),
            step_fn=rts.secant_step,
            tol=tol,
            max_iter=max_iter,
        )
        E_opt = energy(lam_opt)
        return {"lam_opt": lam_opt, "E_opt": E_opt, "method": "secant"}
    except (RuntimeError, ZeroDivisionError):
        # Fall back to bisection with a bracket search
        pass

    # Bracket search: expand interval until sign change in denergy
    lam_a = lam_initial * 0.5
    lam_b = lam_initial * 2.0
    for _ in range(20):
        fa = denergy(lam_a)
        fb = denergy(lam_b)
        if fa * fb <= 0:
            break
        if abs(fa) < abs(fb):
            lam_a -= (lam_b - lam_a)
        else:
            lam_b += (lam_b - lam_a)
        if lam_a <= 0.01:
            lam_a = 0.01
    else:
        raise RuntimeError("Could not bracket the minimum of E[λ]")

    lam_opt = rts.find_root(
        denergy,
        df=None,
        state0=(lam_a, lam_b),
        step_fn=rts.bisection_step,
        tol=tol,
        max_iter=max_iter,
    )
    E_opt = energy(lam_opt)
    return {"lam_opt": lam_opt, "E_opt": E_opt, "method": "bisection"}


def variational_solution(
    alpha: float,
    beta: float,
    gamma: float,
    n_quad: int = 10001,
) -> dict:
    """Full variational solution for ground and first excited states.

    Returns dict with keys: ground, excited.
    """
    L = alpha + gamma

    # Ground state (even)
    ground = optimize_variational(
        L, alpha, beta, gamma,
        trial_fn=trial_ground,
        trial_deriv=trial_ground_derivative,
        lam_initial=0.5,
        n_quad=n_quad,
    )

    # First excited state (odd)
    excited = optimize_variational(
        L, alpha, beta, gamma,
        trial_fn=trial_excited,
        trial_deriv=trial_excited_derivative,
        lam_initial=0.5,
        n_quad=n_quad,
    )

    return {
        "ground": {
            "parity": "even",
            "lam_opt": ground["lam_opt"],
            "E_opt": ground["E_opt"],
            "method": ground["method"],
        },
        "excited": {
            "parity": "odd",
            "lam_opt": excited["lam_opt"],
            "E_opt": excited["E_opt"],
            "method": excited["method"],
        },
    }
