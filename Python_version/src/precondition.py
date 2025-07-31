"""
precondition.py
Preconditioning filter design for multichannel ANC (Python version).

Author: Yongjie Zhuang
"""

import numpy as np
from numpy.linalg import inv
from scipy import signal
from pathlib import Path

from specfact import spectral_factorization
from plotting import plot_frequency_responses, plot_impulse_responses
from invfreqz_ls import invfreqz_ls


def compute_precondition_filters(Fxx: np.ndarray,
                                 h_Ge: np.ndarray,
                                 N_Finv: int,
                                 N_mininv: int,
                                 plot: bool = True,
                                 save_path: Path = None):
    """
    Compute preconditioning filters.

    Parameters
    ----------
    Fxx : ndarray
        Reference filter: shape (N_xx, Nr, Nr)
    h_Ge : ndarray
        Secondary path filter: shape (N_Ge, Ne, Ns)
    N_Finv : int
        Length of Fxx inverse filter
    N_mininv : int
        Length of minimum-phase inverse filter
    plot : bool, optional
        Whether to plot frequency/impulse responses
    save_path : Path or None
        Where to save figures (None disables saving)

    Returns
    -------
    Fxx_inv : ndarray
        Inverse whitening filter, shape (N_Finv, Nr, Nr)
    Ge_min_inv : ndarray
        Inverse minimum-phase filter, shape (N_mininv, Ns, Ns)
    Ge_all : ndarray
        Combined secondary path filter, shape (N_mininv, Ne, Ns)
    """
    if save_path is not None:
        save_path.mkdir(parents=True, exist_ok=True)

    N_xx, Nr, Nr2 = Fxx.shape
    N_Ge, Ne, Ns = h_Ge.shape
    if Nr != Nr2:
        raise ValueError("Fxx must be square in its last two dimensions")

    den_fac = 8
    w = np.linspace(0, np.pi, N_Finv * den_fac)

    # === Fxx frequency response ===
    Fxx_freq = np.zeros((Nr, Nr, len(w)), dtype=complex)
    for ii in range(Nr):
        for jj in range(Nr):
            _, H = signal.freqz(Fxx[:, ii, jj], worN=w)
            Fxx_freq[ii, jj, :] = H

    # Invert at each frequency point
    Fxx_freq_inv = np.zeros_like(Fxx_freq)
    for k in range(len(w)):
        Fxx_freq_inv[:, :, k] = inv(Fxx_freq[:, :, k])

    # Fit FIR inverse filter using custom least-squares invfreqz
    Fxx_inv = np.zeros((N_Finv, Nr, Nr))
    for ii in range(Nr):
        for jj in range(Nr):
            H_target = Fxx_freq_inv[ii, jj, :]
            hn, _ = invfreqz_ls(H_target, w, N_Finv, M=1, wt=np.ones_like(w))
            Fxx_inv[:, ii, jj] = np.real(hn)

    # === Ge_min via spectral factorization ===
    GeGe = np.zeros((2 * N_Ge - 1, Ns, Ns))
    for ii in range(Ns):
        for jj in range(ii, Ns):
            for kk in range(Ne):
                GeGe[:, ii, jj] += np.convolve(
                    np.flip(h_Ge[:, kk, ii]), h_Ge[:, kk, jj]
                )
    GeGe_p = np.zeros((N_Ge, Ns, Ns))
    for ii in range(Ns):
        for jj in range(ii, Ns):
            GeGe_p[:, ii, jj] = GeGe[N_Ge - 1:, ii, jj]
            GeGe_p[:, jj, ii] = GeGe[:N_Ge, ii, jj][::-1]

    Ge_min_H, _, _ = spectral_factorization(GeGe_p)
    Ge_min = np.zeros((N_Ge, Ns, Ns))
    for ii in range(Ns):
        for jj in range(Ns):
            Ge_min[:, ii, jj] = Ge_min_H[:, jj, ii]

    # === Inverse of Ge_min ===
    w_mininv = np.linspace(0, np.pi, N_mininv * den_fac)
    Ge_min_freq = np.zeros((Ns, Ns, len(w_mininv)), dtype=complex)
    for ii in range(Ns):
        for jj in range(Ns):
            _, H = signal.freqz(Ge_min[:, ii, jj], worN=w_mininv)
            Ge_min_freq[ii, jj, :] = H
    for k in range(len(w_mininv)):
        Ge_min_freq[:, :, k] = inv(Ge_min_freq[:, :, k])

    Ge_min_inv = np.zeros((N_mininv, Ns, Ns))
    for ii in range(Ns):
        for jj in range(Ns):
            H_target = 1.0 / Ge_min_freq[ii, jj, :]
            hn, _ = invfreqz_ls(H_target, w_mininv, N_mininv, M=1, wt=np.ones_like(w_mininv))
            Ge_min_inv[:, ii, jj] = np.real(hn)

    # === Ge_all ===
    Ge_freq = np.zeros((Ne, Ns, len(w_mininv)), dtype=complex)
    for ii in range(Ne):
        for jj in range(Ns):
            _, H = signal.freqz(h_Ge[:, ii, jj], worN=w_mininv)
            Ge_freq[ii, jj, :] = H

    Ge_all_freq = np.zeros((Ne, Ns, len(w_mininv)), dtype=complex)
    for k in range(len(w_mininv)):
        Ge_all_freq[:, :, k] = Ge_freq[:, :, k] @ Ge_min_freq[:, :, k]

    Ge_all = np.zeros((N_mininv, Ne, Ns))
    for ii in range(Ne):
        for jj in range(Ns):
            H_target = Ge_all_freq[ii, jj, :]
            hn, _ = invfreqz_ls(H_target, w_mininv, N_mininv, M=1, wt=np.ones_like(w_mininv))
            Ge_all[:, ii, jj] = np.real(hn)

    # === Plotting ===
    if plot:
        if save_path is None:
            save_path = Path.cwd()
        plot_frequency_responses(Fxx_freq, Fxx_inv, "Inverse of Fxx", save_path / "Fxx_freq.png")
        plot_impulse_responses(Fxx_inv, "Inverse of Fxx Impulse", save_path / "Fxx_impulse.png")
        plot_frequency_responses(Ge_min_freq, Ge_min_inv, "Inverse of G_min", save_path / "Gmin_freq.png")
        plot_impulse_responses(Ge_min_inv, "Inverse of G_min Impulse", save_path / "Gmin_impulse.png")
        plot_frequency_responses(Ge_all_freq, Ge_all, "Ge_all", save_path / "Geall_freq.png")
        plot_impulse_responses(Ge_all, "Ge_all Impulse", save_path / "Geall_impulse.png")

    return Fxx_inv, Ge_min_inv, Ge_all
