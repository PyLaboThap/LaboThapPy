import numpy as np
import CoolProp.CoolProp as CP
from labothappy.correlations.void_fraction.void_fraction import compute_void_fraction

# ----------------------------------------------------------------------
# Initialize CoolProp fluid (R410A in this example)
# ----------------------------------------------------------------------
fluid_name = "R410A"
AS = CP.AbstractState("HEOS", fluid_name)

# Test conditions
P = 1.2e6  # pressure [Pa]

d_hyd = 0.008       # hydraulic diameter [m]
m_dot = 300 * (np.pi * d_hyd ** 2 / 4.0)


params = {
    "d_hyd": d_hyd,
}

models = [
    "Homogeneous", "Zivi", "Fauske", "Premoli",
    "Lockhart-Martinelli", "Hughmark",
    "Graham", "Armand-Treschev", "Bankoff", "Rouhani-Axelsson",
    "DIX", "Woldesemayat-Ghajar", "Cioncolini-Thome"
]

qualities = np.linspace(0.05, 0.95, 10)

# ----------------------------------------------------------------------
# Run compute_void_fraction for each model, across qualities
# ----------------------------------------------------------------------
results = {name: [] for name in models}

for x in qualities:
    # Set state with pressure and quality
    AS.update(CP.PQ_INPUTS, P, x)
    for name in models:
        alpha = compute_void_fraction(AS, params, m_dot, void_fraction_model=name)
        results[name].append(alpha)

# ----------------------------------------------------------------------
# Print comparison table
# ----------------------------------------------------------------------
header = f"{'x':>6}" + "".join(f"{name:>24}" for name in models)
print(header)
print("-" * len(header))

for i, x in enumerate(qualities):
    row = f"{x:6.2f}" + "".join(f"{results[name][i]:24.4f}" for name in models)
    print(row)

# ----------------------------------------------------------------------
# Optional: plot to compare visually
# ----------------------------------------------------------------------
try:
    import matplotlib.pyplot as plt

    plt.figure(figsize=(9, 6))
    for name, values in results.items():
        plt.plot(qualities, values, marker='o', label=name)

    plt.xlabel("Vapor quality x [-]")
    plt.ylabel("Void fraction alpha [-]")
    plt.title("Void fraction correlations comparison")
    plt.legend(fontsize=8, ncol=2)
    plt.grid(True)
    plt.tight_layout()
    plt.show()
except ImportError:
    print("\n(matplotlib not installed, skipping plot)")
