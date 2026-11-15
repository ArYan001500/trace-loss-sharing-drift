"""Wavelet coherence between two temperature records (Morlet, omega0 = 6), with the usual
smoothing in time (Gaussian at each scale) and across scales (boxcar). Returns the squared
coherence, the local complex ratio S(W1 W2*) / S(|W2|^2) and the periods."""
import numpy as np

def cwt(y, dt, periods, w0=6.0):
    n = len(y); y = y - y.mean()
    f = np.fft.fftfreq(n, dt) * 2 * np.pi; Y = np.fft.fft(y)
    scales = periods * (w0 + np.sqrt(2 + w0 ** 2)) / (4 * np.pi)
    W = np.empty((len(scales), n), complex)
    for i, s in enumerate(scales):
        psi = np.pi ** -0.25 * np.exp(-0.5 * (s * f - w0) ** 2) * (f > 0) * np.sqrt(2 * np.pi * s / dt)
        W[i] = np.fft.ifft(Y * psi)
    return W, scales

def _smooth(X, scales, dt):
    out = np.empty_like(X); n = X.shape[1]; f = np.fft.fftfreq(n) * 2 * np.pi
    for i, s in enumerate(scales):
        out[i] = np.fft.ifft(np.fft.fft(X[i]) * np.exp(-0.5 * (s / dt) ** 2 * f ** 2))
    k = np.ones(3) / 3
    for j in range(n): out[:, j] = np.convolve(out[:, j], k, mode="same")
    return out

def coherence(y1, y2, dt=1.0, periods=None):
    periods = np.geomspace(8, 200, 24) if periods is None else periods
    W1, sc = cwt(y1, dt, periods); W2, _ = cwt(y2, dt, periods)
    S = lambda X: _smooth(X / sc[:, None], sc, dt)
    S12 = S(W1 * np.conj(W2)); S11 = S(np.abs(W1) ** 2).real; S22 = S(np.abs(W2) ** 2).real
    R2 = np.abs(S12) ** 2 / (S11 * S22)
    return R2, S12 / S22, periods
