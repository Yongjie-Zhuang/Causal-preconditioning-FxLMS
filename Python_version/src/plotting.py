"""
plotting.py
Utility functions for plotting frequency and impulse responses.

Author: Yongjie Zhuang
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import signal


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
