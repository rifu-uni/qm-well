"""Physical constants and dimensionless unit conversion for the filleted well.

All quantities are expressed in terms of the Bohr radius a₀ as the natural
length scale and E₀ = ℏ²/(2m a₀²) as the natural energy scale (one Rydberg
for the electron mass).
"""

import numpy as np

# === Physical Constants (SI) ===

HBAR = 1.054571817e-34          # Reduced Planck constant (J·s)
M_E = 9.1093837015e-31          # Electron rest mass (kg)
A0 = 5.29177210903e-11          # Bohr radius (m)
E_RY = HBAR**2 / (2.0 * M_E * A0**2)  # Rydberg energy ≈ 2.17987×10⁻¹⁸ J = 13.606 eV

# === Dimensionless Conversion ===

def energy_scale(mass: float = M_E) -> float:
    """Return the energy scale E₀ = ℏ²/(2m a₀²) in joules.
    
    Parameters
    ----------
    mass : float
        Particle mass in kg. Defaults to electron mass.
    
    Returns
    -------
    float
        Energy scale in joules.
    """
    return HBAR**2 / (2.0 * mass * A0**2)


def to_dimensionless_energy(E: float, mass: float = M_E) -> float:
    """Convert physical energy (joules) to dimensionless units.
    
    Parameters
    ----------
    E : float
        Energy in joules.
    mass : float
        Particle mass in kg.
    
    Returns
    -------
    float
        Dimensionless energy Ẽ = (2m a₀²/ℏ²) E.
    """
    return E / energy_scale(mass)


def from_dimensionless_energy(E_tilde: float, mass: float = M_E) -> float:
    """Convert dimensionless energy back to joules.
    
    Parameters
    ----------
    E_tilde : float
        Dimensionless energy.
    mass : float
        Particle mass in kg.
    
    Returns
    -------
    float
        Energy in joules.
    """
    return E_tilde * energy_scale(mass)


def to_dimensionless_position(x: float) -> float:
    """Convert physical position (meters) to dimensionless units x̃ = x / a₀."""
    return x / A0


def from_dimensionless_position(x_tilde: float) -> float:
    """Convert dimensionless position back to meters."""
    return x_tilde * A0


# === Default well parameters (dimensionless) ===
# These are sensible defaults for a non-trivial demonstration.

DEFAULT_ALPHA = 3.0     # half-width of the main rectangular section
DEFAULT_BETA = 5.0      # depth of the well
DEFAULT_GAMMA = 1.0     # corner rounding (fillet) radius
