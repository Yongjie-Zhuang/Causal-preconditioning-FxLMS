"""
specfact.py
Spectral factorization for multichannel signals.

Author: Yongjie Zhuang
"""

import numpy as np
from scipy.linalg import solve_discrete_are, cholesky


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
        N_bar[(ii-1)*Nx:ii*Nx, :] = Sxx[Nq + 1 - ii, :, :]

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
        L[:, :, Nq - ii] = g[(ii-1)*Nx:ii*Nx, :]

    # Build F filter
    chol_Re = cholesky(Re, lower=True)
    F = np.zeros((Nq + 1, Nx, Nx))
    F[0, :, :] = chol_Re
    for ii in range(1, Nq + 1):
        F[ii, :, :] = L[:, :, ii - 1] @ chol_Re

    return F, L, Re
