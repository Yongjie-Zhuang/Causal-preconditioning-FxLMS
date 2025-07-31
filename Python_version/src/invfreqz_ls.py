"""
invfreqz_ls.py
Custom least squares version of MATLAB invfreqz.

Author: Yongjie Zhuang
"""

import numpy as np


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
