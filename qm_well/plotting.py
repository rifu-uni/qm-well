"""Visualization for the filleted infinite potential well solution.

Generates publication-quality plots of:
  1. The potential profile V(x)
  2. Energy levels from all methods overlaid on the potential
  3. Wavefunctions superimposed on the energy levels
  4. Comparison of methods (perturbation, variational, FDM, shooting)
"""

import numpy as np
from numpy.typing import NDArray
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from typing import Optional

from .potential import v_tilde, rectangular_v_tilde
from .analytical import (
    trial_ground, trial_ground_derivative,
    trial_excited, trial_excited_derivative,
    variational_energy,
)
from .numerical import solve_fdm


# ============================================================
# Style configuration
# ============================================================

COLORS = {
    "potential": "#2c3e50",
    "potential_rect": "#7f8c8d",
    "energy_fdm": "#e74c3c",
    "energy_shoot": "#3498db",
    "energy_pert": "#f39c12",
    "energy_var": "#2ecc71",
    "wavefunction_even": "#8e44ad",
    "wavefunction_odd": "#e67e22",
    "well_fill": "#ecf0f1",
}

LINEWIDTH = 1.5
FIGSIZE = (10, 6)
DPI = 150


def set_style():
    """Configure matplotlib for clean, publication-style plots."""
    plt.rcParams.update({
        "figure.figsize": FIGSIZE,
        "figure.dpi": DPI,
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 9,
        "lines.linewidth": LINEWIDTH,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def plot_potential(
    alpha: float,
    beta: float,
    gamma: float,
    n_points: int = 2000,
    show_rectangular: bool = True,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Plot the filleted potential V(x) with optional rectangular reference.

    Parameters
    ----------
    alpha, beta, gamma : float
        Well geometry parameters.
    n_points : int
        Number of points for the smooth potential curve.
    show_rectangular : bool
        If True, also plot the reference rectangular well.
    ax : matplotlib Axes or None
        Axes to plot on. Creates new figure if None.

    Returns
    -------
    matplotlib Axes
    """
    if ax is None:
        _, ax = plt.subplots()

    L = alpha + gamma
    x_tilde = np.linspace(-L, L, n_points)
    V = v_tilde(x_tilde, alpha, beta, gamma)

    # Mask infinite values for plotting
    mask = np.isfinite(V)
    ax.plot(x_tilde[mask], V[mask], color=COLORS["potential"],
            label="Filleted potential $\\tilde{V}(\\tilde{x})$", linewidth=2)

    if show_rectangular:
        V_rect = rectangular_v_tilde(x_tilde, alpha, beta, gamma)
        mask_rect = np.isfinite(V_rect)
        ax.plot(x_tilde[mask_rect], V_rect[mask_rect],
                color=COLORS["potential_rect"], linestyle="--",
                label="Rectangular well (reference)", alpha=0.7)

    # Vertical walls
    for xw in [-L, L]:
        ax.axvline(x=xw, color="black", linestyle=":", alpha=0.5)
    ax.axhline(y=0, color="black", linestyle=":", alpha=0.3)

    ax.set_xlabel("$\\tilde{x} = x / a_0$")
    ax.set_ylabel("$\\tilde{V}(\\tilde{x})$")
    ax.set_title(f"Filleted Infinite Potential Well "
                 f"($\\alpha={alpha},\\ \\beta={beta},\\ \\gamma={gamma}$)")
    ax.legend(loc="lower center", ncol=2)
    ax.set_aspect("equal")

    return ax


def plot_energy_levels(
    fdm_results: list[dict],
    alpha: float,
    beta: float,
    gamma: float,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Plot FDM energy levels overlaid on the potential profile.

    Parameters
    ----------
    fdm_results : list[dict]
        Each dict has keys: n, parity, E, psi, x_grid.
    alpha, beta, gamma : float
        Well geometry.
    ax : matplotlib Axes or None

    Returns
    -------
    matplotlib Axes
    """
    if ax is None:
        _, ax = plt.subplots()

    # Plot potential
    L = alpha + gamma
    x_plot = np.linspace(-L, L, 2000)
    V = v_tilde(x_plot, alpha, beta, gamma)
    mask = np.isfinite(V)
    ax.plot(x_plot[mask], V[mask], color=COLORS["potential"],
            linewidth=2.5, zorder=1)

    # Fill the well
    ax.fill_between(x_plot[mask], V[mask], 0,
                     color=COLORS["well_fill"], alpha=0.3, zorder=0)

    # Draw energy level lines
    for result in fdm_results:
        E = result["E"]
        n = result["n"]
        parity = result["parity"]
        color = COLORS["wavefunction_even"] if parity == "even" else COLORS["wavefunction_odd"]

        ax.axhline(y=E, xmin=0.05, xmax=0.3, color=color,
                   linestyle="-", linewidth=1.5, alpha=0.8, zorder=2)
        ax.text(-L * 0.68, E, f"$n={n}$", fontsize=9, color=color,
                va="center", ha="right")

    # Vertical walls
    for xw in [-L, L]:
        ax.axvline(x=xw, color="black", linestyle=":", alpha=0.5)

    ax.set_xlabel("$\\tilde{x} = x / a_0$")
    ax.set_ylabel("$\\tilde{E},\\ \\tilde{V}$")
    ax.set_title("Energy Levels (FDM)")
    ax.set_ylim(beta * -1.3, max(0.5, max(r["E"] for r in fdm_results) * 1.3))
    ax.set_aspect("equal")

    return ax


def plot_wavefunctions(
    fdm_results: list[dict],
    alpha: float,
    beta: float,
    gamma: float,
    n_states: int = 4,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Plot wavefunctions superimposed on energy levels.

    Each wavefunction is scaled and shifted to its energy level.

    Parameters
    ----------
    fdm_results : list[dict]
        FDM results with psi and x_grid.
    alpha, beta, gamma : float
        Well geometry.
    n_states : int
        Number of states to plot.
    ax : matplotlib Axes or None

    Returns
    -------
    matplotlib Axes
    """
    if ax is None:
        _, ax = plt.subplots()

    L = alpha + gamma

    # Plot potential
    x_plot = np.linspace(-L, L, 2000)
    V = v_tilde(x_plot, alpha, beta, gamma)
    mask = np.isfinite(V)
    ax.plot(x_plot[mask], V[mask], color=COLORS["potential"],
            linewidth=2.5, zorder=1)
    ax.fill_between(x_plot[mask], V[mask], 0,
                     color=COLORS["well_fill"], alpha=0.3, zorder=0)

    # Amplitude for visual scaling of wavefunctions
    wave_scale = 0.8 * abs(beta) / (n_states + 2)

    for result in fdm_results[:n_states]:
        E = result["E"]
        parity = result["parity"]
        psi_half = result["psi"]
        x_half = result["x_grid"]

        # Mirror to full domain [-L, L]
        if parity == "even":
            x_full = np.concatenate([-x_half[::-1], x_half])
            psi_full = np.concatenate([psi_half[::-1], psi_half])
        else:
            x_full = np.concatenate([-x_half[::-1], x_half])
            psi_full = np.concatenate([-psi_half[::-1], psi_half])

        # Scale and shift
        psi_scaled = psi_full * wave_scale + E

        color = COLORS["wavefunction_even"] if parity == "even" else COLORS["wavefunction_odd"]
        ax.plot(x_full, psi_scaled, color=color, linewidth=LINEWIDTH,
                zorder=3, label=f"$n={result['n']}$ ({parity})" if result['n'] <= n_states else "")

        # Energy level line
        ax.axhline(y=E, xmin=0.05, xmax=0.25, color=color,
                   linestyle="--", linewidth=0.8, alpha=0.6)

        # Zero line for wavefunction
        ax.plot(x_full, np.full_like(x_full, E), color=color,
                linestyle=":", linewidth=0.5, alpha=0.3, zorder=2)

    # Vertical walls
    for xw in [-L, L]:
        ax.axvline(x=xw, color="black", linestyle=":", alpha=0.5)

    ax.set_xlabel("$\\tilde{x} = x / a_0$")
    ax.set_ylabel("$\\tilde{E},\\ \\tilde{V}$")
    ax.set_title("Wavefunctions and Energy Levels")
    ax.legend(loc="lower right", ncol=2, fontsize=8)
    ax.set_aspect("equal")

    return ax


def plot_method_comparison(
    pert_results: list[dict],
    var_results: dict,
    fdm_results: list[dict],
    shoot_results: list[dict],
    alpha: float,
    beta: float,
    gamma: float,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Compare energy eigenvalues from all four methods.

    Parameters
    ----------
    pert_results : list[dict]
        Perturbation theory results (keys: n, E_total).
    var_results : dict
        Variational results (keys: ground, excited; each has E_opt).
    fdm_results : list[dict]
        FDM results (keys: n, E).
    shoot_results : list[dict]
        Shooting results (keys: n, E).
    alpha, beta, gamma : float
        Well geometry.
    ax : matplotlib Axes or None

    Returns
    -------
    matplotlib Axes
    """
    if ax is None:
        _, ax = plt.subplots()

    n_states = len(fdm_results)
    ns = np.arange(1, n_states + 1)

    # Extract energies
    E_pert = np.array([r["E_total"] for r in pert_results[:n_states]], dtype=float)
    E_fdm = np.array([r["E"] for r in fdm_results], dtype=float)
    E_shoot = np.array([r["E"] for r in shoot_results[:n_states]], dtype=float)

    # Variational energies (only ground and first excited calculated)
    E_var_list = []
    for n in ns:
        if n == 1:
            E_var_list.append(var_results["ground"]["E_opt"])
        elif n == 2:
            E_var_list.append(var_results["excited"]["E_opt"])
        else:
            E_var_list.append(np.nan)
    E_var = np.array(E_var_list, dtype=float)

    # Plot with offsets for clarity
    width = 0.2
    x_offsets = np.array([-1.5, -0.5, 0.5, 1.5]) * width

    labels = ["Perturbation", "Variational", "FDM", "Shooting"]
    colors_list = [COLORS["energy_pert"], COLORS["energy_var"],
                   COLORS["energy_fdm"], COLORS["energy_shoot"]]
    energies_list = [E_pert, E_var, E_fdm, E_shoot]

    for i, (label, color, E_vals) in enumerate(zip(labels, colors_list, energies_list)):
        valid = ~np.isnan(E_vals)
        ax.bar(ns[valid] + x_offsets[i], E_vals[valid], width,
               color=color, label=label, alpha=0.8, edgecolor="white")

    ax.set_xlabel("Quantum number $n$")
    ax.set_ylabel("$\\tilde{E}_n$")
    ax.set_title(f"Method Comparison ($\\alpha={alpha},\\ \\beta={beta},\\ \\gamma={gamma}$)")
    ax.set_xticks(ns)
    ax.legend()
    ax.axhline(y=0, color="black", linestyle=":", alpha=0.3)

    return ax


def plot_variational_convergence(
    alpha: float,
    beta: float,
    gamma: float,
    n_lam: int = 100,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Plot the variational energy E[λ] as a function of λ.

    Shows the minimization landscape for ground and excited states.

    Parameters
    ----------
    alpha, beta, gamma : float
        Well geometry.
    n_lam : int
        Number of λ values to sample.
    ax : matplotlib Axes or None

    Returns
    -------
    matplotlib Axes
    """
    if ax is None:
        _, ax = plt.subplots()

    L = alpha + gamma
    lam_vals = np.linspace(0.1, 5.0, n_lam)

    # Ground state
    E_ground = np.array([
        variational_energy(lam, L, alpha, beta, gamma,
                           trial_ground, trial_ground_derivative,
                           n_quad=5001)
        for lam in lam_vals
    ])

    # Excited state
    E_excited = np.array([
        variational_energy(lam, L, alpha, beta, gamma,
                           trial_excited, trial_excited_derivative,
                           n_quad=5001)
        for lam in lam_vals
    ])

    ax.plot(lam_vals, E_ground, color=COLORS["wavefunction_even"],
            label="Ground state (even)", linewidth=2)
    ax.plot(lam_vals, E_excited, color=COLORS["wavefunction_odd"],
            label="First excited (odd)", linewidth=2)

    # Mark minima
    i_min_g = np.argmin(E_ground)
    i_min_e = np.argmin(E_excited)
    ax.scatter([lam_vals[i_min_g]], [E_ground[i_min_g]],
               color=COLORS["wavefunction_even"], s=80, zorder=5)
    ax.scatter([lam_vals[i_min_e]], [E_excited[i_min_e]],
               color=COLORS["wavefunction_odd"], s=80, zorder=5)

    ax.set_xlabel("$\\lambda$")
    ax.set_ylabel("$\\tilde{E}[\\lambda]$")
    ax.set_title("Variational Energy Functional")
    ax.legend()

    return ax


def plot_all(
    alpha: float,
    beta: float,
    gamma: float,
    pert_results: list[dict],
    var_results: dict,
    fdm_results: list[dict],
    shoot_results: list[dict],
    save_path: Optional[str] = None,
) -> Figure:
    """Generate a comprehensive 4-panel figure summarizing all results.

    Parameters
    ----------
    alpha, beta, gamma : float
        Well geometry.
    pert_results : list[dict]
        Perturbation theory results.
    var_results : dict
        Variational method results.
    fdm_results : list[dict]
        FDM results.
    shoot_results : list[dict]
        Shooting method results.
    save_path : str or None
        If provided, save figure to this path.

    Returns
    -------
    matplotlib Figure
    """
    set_style()
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        f"Filleted Infinite Potential Well  "
        f"($\\alpha={alpha},\\ \\beta={beta},\\ \\gamma={gamma}$)",
        fontsize=14, fontweight="bold",
    )

    # Panel 1: Potential + energy levels
    plot_energy_levels(fdm_results, alpha, beta, gamma, ax=axes[0, 0])
    axes[0, 0].set_title("Potential Profile & Energy Levels")

    # Panel 2: Wavefunctions
    plot_wavefunctions(fdm_results, alpha, beta, gamma,
                       n_states=4, ax=axes[0, 1])
    axes[0, 1].set_title("Wavefunctions")

    # Panel 3: Method comparison
    plot_method_comparison(pert_results, var_results, fdm_results, shoot_results,
                           alpha, beta, gamma, ax=axes[1, 0])
    axes[1, 0].set_title("Energy Comparison Across Methods")

    # Panel 4: Variational landscape
    plot_variational_convergence(alpha, beta, gamma, ax=axes[1, 1])
    axes[1, 1].set_title("Variational $\\tilde{E}[\\lambda]$")

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Figure saved to {save_path}")

    return fig
