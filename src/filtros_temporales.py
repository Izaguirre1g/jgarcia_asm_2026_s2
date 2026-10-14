import numpy as np
from scipy.signal import (buttord, cheb1ord, ellipord, ellip, sos2tf, sos2zpk, sosfreqz)

def calcular_ordenes(fp, f_stop, rp_db, as_db, fs):
    n_butter, _ = buttord(fp, f_stop, rp_db, as_db, fs=fs)
    n_cheb1, _ = cheb1ord(fp, f_stop, rp_db, as_db, fs=fs)
    n_ellip, _ = ellipord(fp, f_stop, rp_db, as_db, fs=fs)
    return {"butterworth": int(n_butter), "chebyshev1": int(n_cheb1), "eliptico": int(n_ellip)}

def disenar_eliptico_sos(fp, f_stop, rp_db, as_db, fs):
    """Devuelve (sos, b, a, orden, wn_hz); wn es el borde de la banda de paso."""
    orden, wn = ellipord(fp, f_stop, rp_db, as_db, fs=fs)
    sos = ellip(orden, rp_db, as_db, wn, btype="lowpass", fs=fs, output="sos")
    b, a = sos2tf(sos)
    return sos, b, a, int(orden), float(wn)

def respuesta_en_frecuencia(sos, fs, n_puntos=4096):
    """Magnitud (dB), fase (rad) y retardo de grupo (muestras)."""
    w, h = sosfreqz(sos, worN=n_puntos, fs=fs)
    mag_db = 20 * np.log10(np.maximum(np.abs(h), 1e-12))
    fase = np.unwrap(np.angle(h))
    retardo = -np.gradient(fase, 2 * np.pi * w / fs)
    return w, mag_db, fase, retardo

def polos_ceros(sos):
    """Estable si todos los polos cumplen |p| < 1."""
    z, p, k = sos2zpk(sos)
    radio_max = float(np.max(np.abs(p)))
    return z, p, k, radio_max, bool(radio_max < 1.0)

def aplicar_sos_manual(sos, x):
    y = np.asarray(x, dtype=np.float64).copy()
    for b0, b1, b2, a0, a1, a2 in sos:
        s1 = s2 = 0.0
        for n in range(len(y)):
            xn = y[n]
            yn = b0 * xn + s1
            s1 = b1 * xn - a1 * yn + s2
            s2 = b2 * xn - a2 * yn
            y[n] = yn
    return y