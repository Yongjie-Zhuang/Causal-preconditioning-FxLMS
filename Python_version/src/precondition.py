"""
precondition.py
Preconditioning filter design for multichannel ANC (Python version).

Author: Yongjie Zhuang
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv
from scipy import signal
from pathlib import Path
from scipy.linalg import solve_discrete_are, cholesky

def plot_frequency_responses(freq_resp, filters, title_prefix, save_path: Path):
    """
    Plot magnitude and phase of frequency responses.

    Parameters
    ----------
    freq_resp : ndarray
        Frequency responses from measurement (complex), shape (Nout, Nin, Nfreq)
    filters : ndarray
        Time-domain filters for comparison, shape (Ntap, Nout, Nin)
    title_prefix : str
        Prefix for figure title
    save_path : Path
        Path to save the figure (PNG)
    """
    Nout, Nin, Nfreq = freq_resp.shape
    Ntap = filters.shape[0]

    fig, axes = plt.subplots(Nout, Nin * 2, figsize=(12, 3 * Nout))
    if Nout == 1 and Nin == 1:
        axes = np.array([[axes[0], axes[1]]])
    elif Nout == 1:
        axes = np.array([axes])
    elif Nin == 1:
        axes = np.array([[a] for a in axes])

    for i in range(Nout):
        for j in range(Nin):
            # --- Measured frequency response ---
            H_target = freq_resp[i, j, :]

            # --- Fitted filter frequency response ---
            w_filter, H_fit = signal.freqz(filters[:, i, j], worN=Nfreq)

            # Magnitude subplot
            idx_mag = axes[i, j * 2]
            idx_mag.plot(np.linspace(0, np.pi, Nfreq),
                         20 * np.log10(np.abs(H_target) + 1e-12),
                         'r-', label='Measured')
            idx_mag.plot(w_filter,
                         20 * np.log10(np.abs(H_fit) + 1e-12),
                         'b--', label='Fitted')
            idx_mag.set_title(f'{title_prefix} Mag: ({i + 1},{j + 1})')
            idx_mag.set_xlabel('Frequency (rad/sample)')
            idx_mag.set_ylabel('Magnitude (dB)')
            idx_mag.grid(True)

            # Phase subplot
            idx_phase = axes[i, j * 2 + 1]
            idx_phase.plot(np.linspace(0, np.pi, Nfreq),
                           np.angle(H_target),
                           'r-', label='Measured')
            idx_phase.plot(w_filter,
                           np.angle(H_fit),
                           'b--', label='Fitted')
            idx_phase.set_title(f'{title_prefix} Phase: ({i + 1},{j + 1})')
            idx_phase.set_xlabel('Frequency (rad/sample)')
            idx_phase.set_ylabel('Phase (rad)')
            idx_phase.set_ylim([-np.pi, np.pi])
            idx_phase.grid(True)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right')
    fig.suptitle(f'{title_prefix} Frequency Response')
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(save_path, dpi=300)
    plt.close(fig)


def plot_impulse_responses(filters, title_prefix, save_path: Path):
    """
    Plot impulse responses.

    Parameters
    ----------
    filters : ndarray
        Time-domain filters, shape (Ntap, Nout, Nin)
    title_prefix : str
        Prefix for figure title
    save_path : Path
        Path to save the figure (PNG)
    """
    Ntap, Nout, Nin = filters.shape
    fig, axes = plt.subplots(Nout, Nin, figsize=(8, 3 * Nout))
    if Nout == 1 and Nin == 1:
        axes = np.array([[axes]])
    elif Nout == 1:
        axes = np.array([axes])
    elif Nin == 1:
        axes = np.array([[a] for a in axes])

    for i in range(Nout):
        for j in range(Nin):
            idx = axes[i, j]
            idx.plot(filters[:, i, j], linewidth=2)
            idx.set_title(f'{title_prefix} Impulse: ({i + 1},{j + 1})')
            idx.set_xlabel('Samples')
            idx.set_ylabel('Amplitude')
            idx.grid(True)

    fig.suptitle(f'{title_prefix} Impulse Response')
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(save_path, dpi=300)
    plt.close(fig)


def invfreqz_ls(H, wi, N, M=1, wt=None):
    """
    Least-squares approximation of invfreqz.
    
    Parameters
    ----------
    H : ndarray
        Target frequency response (complex), shape (K,)
    wi : ndarray
        Frequency grid (rad/sample), shape (K,)
    N : int
        Number of numerator coefficients
    M : int
        Number of denominator coefficients (default=1 -> FIR)
    wt : ndarray or None
        Optional weighting array (default = ones)
        
    Returns
    -------
    hn : ndarray
        Numerator coefficients
    hm : ndarray
        Denominator coefficients
    """
    if wt is None:
        wt = np.ones_like(H, dtype=float)

    # Apply weighting
    H_wt = H * wt

    # Numerator Fourier matrix
    FN = wt[:, None] * np.exp(-1j * np.outer(wi, np.arange(N)))

    if M == 1:
        # FIR case (denominator = 1)
        hm = np.array([1.0])
        FN_H = FN.conj().T
        F_H_F = np.real(FN_H @ FN)
        F_H_H = np.real(FN_H @ H_wt)
        try:
            hn = np.linalg.solve(F_H_F, F_H_H)
        except np.linalg.LinAlgError:
            # Use pseudoinverse if singular
            hn = np.linalg.pinv(F_H_F) @ F_H_H
    else:
        # IIR case (denominator coefficients also estimated)
        FM = -H_wt[:, None] * np.exp(-1j * np.outer(wi, np.arange(M)))
        B = np.concatenate((FM[:, 1:], FN), axis=1)
        b = FM[:, 0]
        B_H = B.conj().T
        B_H_B = np.real(B_H @ B)
        B_H_b = np.real(B_H @ b)
        temp = np.linalg.solve(B_H_B, -B_H_b)
        hm = np.zeros(M)
        hm[0] = 1
        hm[1:] = temp[0:M - 1]
        hn = temp[M - 1:]

    return hn, hm

def spectral_factorization(Sxx: np.ndarray):
    """
    Multichannel spectral factorization.
    Equivalent to MATLAB mul_specFact.

    Parameters
    ----------
    Sxx : ndarray
        Power spectral density polynomial matrix.
        Shape: (q+1, Nx, Nx), indexed from z^0 to z^-q.

    Returns
    -------
    F : ndarray
        Spectral factor filter, shape (q+1, Nx, Nx)
    L : ndarray
        Gain matrix, shape (Nx, Nx, q)
    Re : ndarray
        Residual matrix, shape (Nx, Nx)
    """
    Nq = Sxx.shape[0] - 1
    Nx = Sxx.shape[1]
    if Nx != Sxx.shape[2]:
        raise ValueError("Sxx must be square (Nx x Nx) for each coefficient")

    # Construct Riccati matrices
    F_ric = np.kron(np.diag(np.ones(Nq - 1), -1), np.eye(Nx))
    h_vec = np.zeros(Nq)
    h_vec[-1] = 1  # Last element is 1
    h = np.kron(h_vec, np.eye(Nx))
    N_bar = np.zeros((Nx * Nq, Nx))
    for ii in range(Nq):
        N_bar[ii*Nx:(ii+1)*Nx, :] = Sxx[Nq - ii, :, :]

    R0 = Sxx[0, :, :]

    # Solve Discrete Algebraic Riccati Equation (DARE)
    Sigma = solve_discrete_are(F_ric.T, h.T,
                               np.zeros((Nx * Nq, Nx * Nq)),
                               -R0, e = np.eye(Nx * Nq), s=-N_bar)
    Re = R0 - h @ Sigma @ h.T
    # Ensure Re is positive definite
    if np.any(np.linalg.eigvals(Re) <= 0):
        Re = Re + 1e-6 * np.eye(Nx)
    g = (N_bar - F_ric @ Sigma @ h.T) @ np.linalg.inv(Re)

    # Reshape to L
    L = np.zeros((Nx, Nx, Nq))
    for ii in range(Nq):
        L[:, :, Nq - 1 - ii] = g[ii*Nx:(ii+1)*Nx, :]

    # Build F filter
    chol_Re = cholesky(Re, lower=True)
    F = np.zeros((Nq + 1, Nx, Nx))
    F[0, :, :] = chol_Re
    for ii in range(1, Nq + 1):
        F[ii, :, :] = L[:, :, ii - 1] @ chol_Re

    return F, L, Re

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
            GeGe_p[:, jj, ii] = GeGe[N_Ge - 1::-1, ii, jj]

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
    # Invert frequency response
    Ge_min_freq_inv = np.zeros_like(Ge_min_freq)
    for k in range(len(w_mininv)):
            Ge_min_freq_inv[:, :, k] = inv(Ge_min_freq[:, :, k])

    Ge_min_inv = np.zeros((N_mininv, Ns, Ns))
    for ii in range(Ns):
        for jj in range(Ns):
            H_target = Ge_min_freq_inv[ii, jj, :]
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
        Ge_all_freq[:, :, k] = Ge_freq[:, :, k] @ Ge_min_freq_inv[:, :, k]

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
        plot_frequency_responses(Fxx_freq_inv, Fxx_inv, "Inverse of Fxx Frequency Response", save_path / "Inverse_Fxx_FreqResp.png")
        plot_impulse_responses(Fxx_inv, "Inverse of Fxx Impulse Response", save_path / "Inverse_Fxx_ImpulseResp.png")
        plot_frequency_responses(Ge_min_freq_inv, Ge_min_inv, "Inverse of G_min Frequency Response", save_path / "Inverse_Gmin_FreqResp.png")
        plot_impulse_responses(Ge_min_inv, "Inverse of G_min Impulse Response", save_path / "Inverse_Gmin_ImpulseResp.png")
        plot_frequency_responses(Ge_all_freq, Ge_all, "Ge_all Frequency Response", save_path / "G_all_FreqResp.png")
        plot_impulse_responses(Ge_all, "Ge_all Impulse Response", save_path / "G_all_ImpulseResp.png")

    return Fxx_inv, Ge_min_inv, Ge_all
