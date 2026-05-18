#!/usr/bin/env python3
"""Main driver for the Filleted Infinite Potential Well solution.

Runs all methods (analytical and numerical), generates plots, and prints
a comprehensive summary of results.

Usage:
    python scripts/solve.py [--alpha A] [--beta B] [--gamma G] [--save]
"""

import sys
import os
import argparse
import time
import matplotlib
matplotlib.use("Agg")

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from qm_well.units import (
    A0, E_RY, HBAR, M_E,
    DEFAULT_ALPHA, DEFAULT_BETA, DEFAULT_GAMMA,
)
from qm_well.potential import v_tilde
from qm_well.analytical import (
    perturbation_energies,
    variational_solution,
)
from qm_well.numerical import (
    solve_fdm,
    find_all_eigenvalues_shooting,
)
from qm_well.plotting import plot_all


# ============================================================
# Pretty printing
# ============================================================

SEPARATOR = "=" * 65


def print_header(text: str) -> None:
    print(f"\n{SEPARATOR}")
    print(f"  {text}")
    print(SEPARATOR)


def print_energy_table(
    title: str,
    results: list[dict],
    key: str = "E",
    extra_keys: list[str] | None = None,
) -> None:
    """Print a formatted energy table.

    Parameters
    ----------
    title : str
        Table title.
    results : list[dict]
        Each dict must have keys: n, parity, and the energy key.
    key : str
        Key for the energy value.
    extra_keys : list[str] or None
        Additional keys to display.
    """
    print(f"\n  {title}")
    print(f"  {'n':>3s}  {'Parity':>6s}  {'Ẽ (dimensionless)':>20s}", end="")
    if extra_keys:
        for ek in extra_keys:
            print(f"  {ek:>20s}", end="")
    print()
    print(f"  {'-'*3}  {'-'*6}  {'-'*20}", end="")
    if extra_keys:
        for _ in extra_keys:
            print(f"  {'-'*20}", end="")
    print()

    for row in results:
        E_val = row.get(key, float("nan"))
        print(f"  {row['n']:3d}  {row['parity']:>6s}  {E_val:20.8f}", end="")
        if extra_keys:
            for ek in extra_keys:
                val = row.get(ek, float("nan"))
                print(f"  {val:20.8f}", end="")
        print()


def print_comparison(
    pert: list[dict],
    var: dict,
    fdm: list[dict],
    shoot: list[dict],
) -> None:
    """Print a comparison table of all methods."""
    print(f"\n  {'n':>3s}  {'Parity':>6s}  {'Perturb.':>14s}  {'Variational':>14s}  "
          f"{'FDM':>14s}  {'Shooting':>14s}")
    print(f"  {'-'*3}  {'-'*6}  {'-'*14}  {'-'*14}  {'-'*14}  {'-'*14}")

    n_states = max(len(pert), len(fdm), len(shoot))
    for i in range(n_states):
        n = i + 1
        parity = "even" if n % 2 == 1 else "odd"

        # Perturbation
        if i < len(pert):
            ep = f"{pert[i]['E_total']:14.8f}"
        else:
            ep = f"{'—':>14s}"

        # Variational (only ground n=1 and excited n=2)
        if n == 1:
            ev = f"{var['ground']['E_opt']:14.8f}"
        elif n == 2:
            ev = f"{var['excited']['E_opt']:14.8f}"
        else:
            ev = f"{'—':>14s}"

        # FDM
        if i < len(fdm):
            ef = f"{fdm[i]['E']:14.8f}"
        else:
            ef = f"{'—':>14s}"

        # Shooting
        if i < len(shoot):
            es = f"{shoot[i]['E']:14.8f}"
        else:
            es = f"{'—':>14s}"

        print(f"  {n:3d}  {parity:>6s}  {ep}  {ev}  {ef}  {es}")


# ============================================================
# Main
# ============================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filleted Infinite Potential Well — Quantum Mechanics Solution",
    )
    parser.add_argument("--alpha", type=float, default=DEFAULT_ALPHA,
                        help="Half-width parameter (default: 3.0)")
    parser.add_argument("--beta", type=float, default=DEFAULT_BETA,
                        help="Depth parameter (default: 5.0)")
    parser.add_argument("--gamma", type=float, default=DEFAULT_GAMMA,
                        help="Fillet radius parameter (default: 1.0)")
    parser.add_argument("--n-states", type=int, default=6,
                        help="Number of bound states to compute (default: 6)")
    parser.add_argument("--n-fdm", type=int, default=2000,
                        help="FDM grid points (default: 2000)")
    parser.add_argument("--shoot-step", type=float, default=0.002,
                        help="Shooting method step size (default: 0.002)")
    parser.add_argument("--save", action="store_true",
                        help="Save plots to output/ directory")
    parser.add_argument("--no-plots", action="store_true",
                        help="Skip plot generation")
    args = parser.parse_args()

    alpha, beta, gamma = args.alpha, args.beta, args.gamma
    n_states = args.n_states

    # Validate parameters
    if gamma >= alpha:
        print("ERROR: Must have γ < α (fillet radius less than half-width)")
        sys.exit(1)
    if 2 * gamma >= beta:
        print("ERROR: Must have 2γ < β (fillet corners must not overlap)")
        sys.exit(1)

    # ============================================================
    # Parameter summary
    # ============================================================
    print_header("Filleted Infinite Potential Well — Quantum Solution")
    print(f"\n  Parameters:")
    print(f"    α (half-width)     = {alpha:.2f} a₀")
    print(f"    β (depth)          = {beta:.2f} E₀")
    print(f"    γ (fillet radius)  = {gamma:.2f} a₀")
    print(f"    Well width (2L)    = {2*(alpha+gamma):.2f} a₀")
    print(f"    Well depth         = {-beta:.2f} E₀")
    print(f"\n  Physical scale:")
    print(f"    a₀ (Bohr radius)   = {A0:.4e} m")
    print(f"    E₀ = ℏ²/(2mₑ a₀²)  = {E_RY:.4e} J = {E_RY/1.602176634e-19:.3f} eV")
    print(f"    For electron mass: E₀ = 1 Rydberg ≈ 13.606 eV")

    # ============================================================
    # 1. Analytical: Perturbation Theory
    # ============================================================
    print_header("1. First-Order Perturbation Theory")
    t0 = time.perf_counter()
    pert_results = perturbation_energies(alpha, beta, gamma, n_states=n_states, n_quad=10001)
    t_pert = time.perf_counter() - t0
    print_energy_table("Perturbation Theory Energies", pert_results,
                       key="E_total", extra_keys=["E0", "E1"])
    print(f"\n  Computation time: {t_pert:.2f} s")

    # ============================================================
    # 2. Analytical: Variational Method
    # ============================================================
    print_header("2. Variational Method")
    t0 = time.perf_counter()
    var_results = variational_solution(alpha, beta, gamma, n_quad=10001)
    t_var = time.perf_counter() - t0

    for state_name, state_data in var_results.items():
        print(f"\n  {state_name.capitalize()} state ({state_data['parity']} parity):")
        print(f"    Optimal λ       = {state_data['lam_opt']:.6f}")
        print(f"    Ẽ[λ_opt]        = {state_data['E_opt']:.8f}")
        print(f"    Root method      = {state_data['method']}")
    print(f"\n  Computation time: {t_var:.2f} s")

    # ============================================================
    # 3. Numerical: Finite Difference Method
    # ============================================================
    print_header("3. Finite Difference Method (FDM)")
    t0 = time.perf_counter()
    fdm_results = solve_fdm(alpha, beta, gamma, N=args.n_fdm, n_states=n_states)
    t_fdm = time.perf_counter() - t0

    if not fdm_results:
        print("\n  No bound states found (E < 0) for these parameters.")
        print(f"  Try increasing β (depth) or α (width).")
        print(f"\n{SEPARATOR}")
        print("  Done.")
        print(f"{SEPARATOR}\n")
        return

    print_energy_table("FDM Energies", fdm_results)
    print(f"\n  Grid points: N = {args.n_fdm}")
    print(f"  Grid spacing: h = {(alpha+gamma)/args.n_fdm:.6f} a₀")
    print(f"  Computation time: {t_fdm:.2f} s")

    # ============================================================
    # 4. Numerical: Shooting Method (catkit)
    # ============================================================
    print_header("4. Shooting Method (catkit RK4 + Bisection)")
    t0 = time.perf_counter()
    n_shoot = min(n_states, 4)  # Shooting is slower, limit to 4 states
    try:
        shoot_results = find_all_eigenvalues_shooting(
            alpha, beta, gamma, n_states=n_shoot, h=args.shoot_step,
        )
        t_shoot = time.perf_counter() - t0
        print_energy_table("Shooting Method Energies", shoot_results)
        print(f"\n  Step size: h = {args.shoot_step:.4f}")
        print(f"  Computation time: {t_shoot:.2f} s")
    except RuntimeError as e:
        print(f"\n  WARNING: Shooting method failed: {e}")
        print(f"  (This can happen for higher excited states with fine-tuning needed)")
        shoot_results = []
        t_shoot = 0.0

    # ============================================================
    # 5. Method Comparison
    # ============================================================
    print_header("5. Method Comparison")
    if shoot_results:
        print_comparison(pert_results, var_results, fdm_results, shoot_results)
    else:
        print_comparison(pert_results, var_results, fdm_results, fdm_results)

    # ============================================================
    # 6. Physical Energies (for electron mass)
    # ============================================================
    print_header("6. Physical Energies (FDM, electron mass)")
    print(f"\n  Converting with E₀ = {E_RY/1.602176634e-19:.3f} eV:")
    print(f"\n  {'n':>3s}  {'Parity':>6s}  {'Ẽ (dimless)':>14s}  {'E (eV)':>14s}  {'E (J)':>18s}")
    print(f"  {'-'*3}  {'-'*6}  {'-'*14}  {'-'*14}  {'-'*18}")
    for row in fdm_results[:n_states]:
        E_tilde = row["E"]
        E_eV = E_tilde * (E_RY / 1.602176634e-19)
        E_J = E_tilde * E_RY
        print(f"  {row['n']:3d}  {row['parity']:>6s}  {E_tilde:14.8f}  {E_eV:14.6f}  {E_J:18.8e}")

    # ============================================================
    # 7. Plots
    # ============================================================
    if not args.no_plots:
        print_header("7. Generating Plots")

        # Create output directory if saving
        save_path = None
        if args.save:
            os.makedirs("output", exist_ok=True)
            save_path = f"output/well_alpha{alpha}_beta{beta}_gamma{gamma}.png"

        try:
            fig = plot_all(
                alpha, beta, gamma,
                pert_results, var_results, fdm_results, shoot_results or fdm_results,
                save_path=save_path,
            )
            if save_path:
                print(f"\n  Plot saved to: {save_path}")
            else:
                print(f"\n  Plots generated. Use --save to save to file.")
                print(f"  Use --no-plots to skip plot generation.")
        except Exception as e:
            print(f"\n  WARNING: Plot generation failed: {e}")

    # ============================================================
    # Summary
    # ============================================================
    print_header("Summary")
    print(f"\n  Ground state energy (FDM):  {fdm_results[0]['E']:.8f} E₀")
    if len(fdm_results) > 1:
        print(f"  First excited (FDM):        {fdm_results[1]['E']:.8f} E₀")
    print(f"  Number of bound states:      {len(fdm_results)} (E < 0)")

    print(f"\n  Method agreement (|ΔE| max for first {n_shoot} states):")
    for i in range(min(n_shoot, len(fdm_results))):
        if shoot_results and i < len(shoot_results):
            diff = abs(fdm_results[i]["E"] - shoot_results[i]["E"])
            print(f"    FDM vs Shooting  n={i+1}:  {diff:.2e}")

    print(f"\n  Total computation time: {t_pert + t_var + t_fdm + t_shoot:.2f} s")
    print(f"\n{SEPARATOR}")
    print("  Done.")
    print(f"{SEPARATOR}\n")


if __name__ == "__main__":
    main()
