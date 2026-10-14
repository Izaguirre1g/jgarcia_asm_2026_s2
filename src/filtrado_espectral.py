import numpy as np
from scipy.signal import welch, find_peaks
from scipy.ndimage import median_filter

def ventana_raiz_hann(n):
    "Hann que suma 1 con 50 % de solape."
    return np.sqrt(0.5 * (1 - np.cos(2 * np.pi * np.arange(n) / n)))


def estimar_ruido_e_interferencia(x_cal, fs, n_bloque, n_est=8192, kernel_hz=250.0, umbral_pico_db=15.0, tol_armonico_hz=5.0):

    # --- piso de ruido a la resolución del bloque
    _, p_blk = welch(x_cal, fs=fs, window=ventana_raiz_hann(n_bloque), nperseg=n_bloque, noverlap=n_bloque // 2, scaling="spectrum")
    k_blk = max(3, int(kernel_hz / (fs / n_bloque)) | 1)
    pn = median_filter(p_blk, size=k_blk, mode="nearest")

    # --- detección de tonos con mayor resolución
    n_est = min(n_est, len(x_cal))
    f, p = welch(x_cal, fs=fs, window="hann", nperseg=n_est, noverlap=n_est // 2,
                 scaling="spectrum")
    piso = median_filter(p, size=max(3, int(kernel_hz / (fs / n_est)) | 1), mode="nearest")
    idx, _ = find_peaks(p, height=piso * 10 ** (umbral_pico_db / 10))

    # interpolación parabólica del máximo (log-potencia) para precisión sub-bin
    fpk = []
    for i in idx:
        a, b, c = np.log(p[i - 1:i + 2] + 1e-30)
        d = 0.5 * (a - c) / (a - 2 * b + c)
        fpk.append(f[i] + d * (fs / n_est))
    fpk = np.array(sorted(fpk))
    interf = []
    if len(fpk) > 0:
        f0 = fpk[0]
        for fp in fpk:
            n_arm = max(1, round(fp / f0))
            if abs(fp - n_arm * f0) > tol_armonico_hz:
                interf.append(float(fp))
    return pn, interf, {"picos_Hz": fpk.round(1).tolist(),
                        "f0_Hz": float(fpk[0]) if len(fpk) else None}


def reducir_ruido_fft(x, fs, pn, f_interf, n_bloque=1024, alpha=2.0, g_min=0.1, ancho_bins=2):

    hop = n_bloque // 2
    w = ventana_raiz_hann(n_bloque)
    x = np.asarray(x, dtype=np.float64)
    n_pad = (-(len(x) - n_bloque)) % hop if len(x) > n_bloque else n_bloque - len(x)
    xp = np.concatenate([np.zeros(hop), x, np.zeros(n_pad + n_bloque)])  # relleno en bordes
    y = np.zeros_like(xp)
    freqs = np.fft.rfftfreq(n_bloque, 1 / fs)
    escala = np.sum(w) ** 2          # misma normalización que welch(scaling='spectrum')
    mascara_int = np.zeros(len(freqs), dtype=bool)

    for fi in f_interf:
        mascara_int |= np.abs(freqs - fi) <= ancho_bins * fs / n_bloque
    gains = []

    for ini in range(0, len(xp) - n_bloque + 1, hop):
        X = np.fft.rfft(xp[ini:ini + n_bloque] * w)
        P = np.abs(X) ** 2 / escala
        P[1:-1] *= 2                 # potencia de un solo lado, como welch
        g2 = np.maximum(1 - alpha * pn / np.maximum(P, 1e-30), g_min ** 2)
        g = np.sqrt(g2)
        g[mascara_int] = g_min
        gains.append(g)
        y[ini:ini + n_bloque] += np.fft.irfft(X * g, n=n_bloque) * w
        
    return y[hop:hop + len(x)], np.array(gains)