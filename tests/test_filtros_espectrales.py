import numpy as np

from src.filtrado_espectral import estimar_ruido_e_interferencia, reducir_ruido_fft

FS, N = 16000, 1024


def _senales():
    t = np.arange(6 * FS) / FS
    limpia = sum(a * np.sin(2 * np.pi * f * t) for f, a in [(440, 1.0), (880, 0.55), (1320, 0.3)])
    interf = 0.66 * np.sin(2 * np.pi * 1750 * t)
    return limpia, limpia + interf


def _snr(s, y):
    return 10 * np.log10(np.sum(s ** 2) / np.sum((s - y) ** 2))


def test_reconstruccion_perfecta_sin_modificar():
    _, x = _senales()
    y, _ = reducir_ruido_fft(x, FS, np.zeros(N // 2 + 1), [], N, alpha=0.0, g_min=1.0)
    assert np.max(np.abs(y - x)) < 1e-12


def test_salida_real_y_misma_longitud():
    _, x = _senales()
    pn, fint, _ = estimar_ruido_e_interferencia(x[:2 * FS], FS, N)
    y, _ = reducir_ruido_fft(x, FS, pn, fint, N)
    assert y.shape == x.shape and np.isrealobj(y)


def test_detecta_interferencia_fuera_de_armonicos():
    _, x = _senales()
    _, fint, _ = estimar_ruido_e_interferencia(x[:2 * FS], FS, N)
    assert len(fint) == 1 and abs(fint[0] - 1750) < 5


def test_mejora_snr():
    s, x = _senales()
    pn, fint, _ = estimar_ruido_e_interferencia(x[:2 * FS], FS, N)
    y, _ = reducir_ruido_fft(x, FS, pn, fint, N)
    assert _snr(s, y) > _snr(s, x) + 6.0