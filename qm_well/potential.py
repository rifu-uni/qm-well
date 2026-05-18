"""The filleted infinite potential well potential V(x).

Defined piecewise with circular-arc fillets at the four corners of a
rectangular well.  All functions use dimensionless coordinates
x̃ = x / a₀ and dimensionless potential Ṽ = (2m a₀² / ℏ²) V.
"""

import numpy as np
from numpy.typing import NDArray


def v_tilde(
    x_tilde: float | NDArray[np.floating],
    alpha: float,
    beta: float,
    gamma: float,
) -> float | NDArray[np.floating]:
    """Dimensionless potential Ṽ(x̃) for the filleted infinite well.

    The well is symmetric about x̃ = 0.  The right half (x̃ ≥ 0) is:

        Ṽ(x̃) = -β                             0 ≤ x̃ ≤ α-γ
        Ṽ(x̃) = -(β-γ) - √[γ² - (x̃-α+γ)²]     α-γ < x̃ ≤ α
        Ṽ(x̃) = -γ + √[γ² - (x̃-α-γ)²]          α < x̃ ≤ α+γ
        Ṽ(x̃) = ∞                              x̃ > α+γ

    The left half follows by symmetry Ṽ(-x̃) = Ṽ(x̃).

    Parameters
    ----------
    x_tilde : float or ndarray
        Dimensionless position(s).
    alpha : float
        Half-width of the main rectangular section.
    beta : float
        Depth of the well.
    gamma : float
        Corner rounding (fillet) radius.  Must satisfy γ < α, 2γ < β.

    Returns
    -------
    float or ndarray
        Dimensionless potential Ṽ(x̃).  Returns np.inf where x̃ > α+γ.
    """
    x = np.abs(np.asarray(x_tilde, dtype=float))
    scalar_input = np.ndim(x) == 0

    result = np.full_like(x, np.inf, dtype=float)

    # Region I: flat bottom  0 ≤ x ≤ α-γ
    mask1 = x <= alpha - gamma
    result[mask1] = -beta

    # Region II: lower fillet  α-γ < x ≤ α
    mask2 = (x > alpha - gamma) & (x <= alpha)
    dx2 = x[mask2] - (alpha - gamma)
    result[mask2] = -(beta - gamma) - np.sqrt(gamma**2 - dx2**2)

    # Region III: upper fillet  α < x ≤ α+γ
    mask3 = (x > alpha) & (x <= alpha + gamma)
    dx3 = x[mask3] - (alpha + gamma)
    result[mask3] = -gamma + np.sqrt(gamma**2 - dx3**2)

    if scalar_input:
        return float(result.item())
    return result


def v_tilde_derivative(
    x_tilde: float,
    alpha: float,
    beta: float,
    gamma: float,
) -> float:
    """Derivative dṼ/dx̃ (for x̃ ≥ 0) used in ODE integration.
    
    Returns the analytical derivative for stability in the shooting method.
    The derivative is zero in the flat region and infinite at x = α+γ.
    """
    x = abs(x_tilde)
    # Region I
    if x <= alpha - gamma:
        return 0.0
    # Region II
    elif x <= alpha:
        dx = x - (alpha - gamma)
        denom = np.sqrt(max(gamma**2 - dx**2, 1e-30))
        return dx / denom  # chain rule: derivative of -√(γ²-dx²) = dx/√(γ²-dx²)
    # Region III
    elif x < alpha + gamma:
        dx = x - (alpha + gamma)
        denom = np.sqrt(max(gamma**2 - dx**2, 1e-30))
        return -dx / denom  # derivative of +√(γ²-dx²) = -dx/√(γ²-dx²)
    else:
        return 0.0  # infinite wall, derivative not meaningful


def rectangular_v_tilde(
    x_tilde: float | NDArray[np.floating],
    alpha: float,
    beta: float,
    gamma: float,
) -> float | NDArray[np.floating]:
    """Reference rectangular well (no fillets) for perturbation theory.

    V_rect = -β for |x̃| < α,  0 for α < |x̃| < α+γ,  ∞ for |x̃| > α+γ.
    """
    x = np.abs(np.asarray(x_tilde, dtype=float))
    scalar_input = np.ndim(x) == 0

    result = np.full_like(x, np.inf, dtype=float)
    mask_deep = x <= alpha
    mask_shallow = (x > alpha) & (x <= alpha + gamma)

    result[mask_deep] = -beta
    result[mask_shallow] = 0.0

    if scalar_input:
        return float(result.item())
    return result


def perturbation_delta_v(
    x_tilde: float | NDArray[np.floating],
    alpha: float,
    beta: float,
    gamma: float,
) -> float | NDArray[np.floating]:
    """Perturbation ΔṼ = Ṽ_filleted - Ṽ_rectangular.

    This is the difference that represents the effect of the filleted corners.
    """
    vf = v_tilde(x_tilde, alpha, beta, gamma)
    vr = rectangular_v_tilde(x_tilde, alpha, beta, gamma)
    return vf - vr
