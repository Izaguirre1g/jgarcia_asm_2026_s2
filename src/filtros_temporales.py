import numpy as np
from scipy.signal import (buttord, cheb1ord, ellipord, ellip,
                          sos2tf, sos2zpk, sosfreqz)


def calcular_ordenes(fp, f_stop, rp_db, as_db, fs):
    n_butter, _ = buttord(fp, f_stop, rp_db, as_db, fs=fs)
    n_cheb1, _ = cheb1ord(fp, f_stop, rp_db, as_db, fs=fs)
    n_ellip, _ = ellipord(fp, f_stop, rp_db, as_db, fs=fs)
    return {"butterworth": int(n_butter),
            "chebyshev1": int(n_cheb1),
            "eliptico": int(n_ellip)}


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