"""
example.py
Demonstration of preconditioning filter design.

Author: Yongjie Zhuang (Python refactor)
"""

import numpy as np
from pathlib import Path
import sys
import os

# --- Add src to Python path ---
this_dir = Path(__file__).resolve().parent
src_path = this_dir.parent / "src"
sys.path.append(str(src_path))

from precondition import compute_precondition_filters

# === Paths ===
save_path = this_dir / "figures"
save_path.mkdir(exist_ok=True)

# === Parameters (same as in the MATLAB paper example) ===
Nr, Ne, Ns = 2, 2, 2
N_Finv, N_mininv = 10, 10
plot_flag = True

fs = 1000
c0 = 340
M = np.array([[0.75, 0.3], [0.3, 1.0]])
l1, l2 = 2.0, 3.0
N1 = int(round(l1 * fs / c0))  # expected 6
N2 = int(round(l2 * fs / c0))  # expected 9

# === Build system matrices ===
Fxx = np.zeros((16, Nr, Nr))
Fxx[0, :, :] = M

N_Ge = max(16, max(N1, N2))
h_Ge = np.zeros((N_Ge, Ne, Ns))
h_Ge[N1, 0, 0] = 1
h_Ge[N2, 0, 1] = l1 / l2
h_Ge[N2, 1, 0] = l1 / l2
h_Ge[N1, 1, 1] = 1

# === Compute preconditioning filters ===
Fxx_inv, Ge_min_inv, Ge_all = compute_precondition_filters(
    Fxx, h_Ge, N_Finv, N_mininv, plot=plot_flag, save_path=save_path
)

# === Display results ===
print("Fxx_inv:\n", Fxx_inv)
print("Ge_min_inv:\n", Ge_min_inv)
print("Ge_all:\n", Ge_all)
print(f"Preconditioning filter computation complete. Figures saved to: {save_path}")
