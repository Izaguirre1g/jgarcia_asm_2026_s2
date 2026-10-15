import numpy as np
from scipy.signal import sosfilt, lfilter

from src.filtros_temporales import (disenar_eliptico_sos, aplicar_sos_manual, respuesta_en_frecuencia, polos_ceros)

FS, FP, FSTOP, RP, AS = 16000, 1500, 2000, 1.0, 40.0
TOL = 1e-10


def _filtro():
    return disenar_eliptico_sos(FP, FSTOP, RP, AS, FS)


def test_manual_vs_sosfilt_y_lfilter():
    sos, b, a, _, _ = _filtro()
    x = np.random.default_rng(0).standard_normal(20000)
    y = aplicar_sos_manual(sos, x)
    assert np.max(np.abs(y - sosfilt(sos, x))) < TOL
    assert np.max(np.abs(y - lfilter(b, a, x))) < TOL


def test_estabilidad():
    sos, *_ = _filtro()
    radio, estable = polos_ceros(sos)[3:5]
    assert estable and radio < 1.0


def test_especificacion_en_frecuencia():
    sos, *_ = _filtro()
    f, mag, _, _ = respuesta_en_frecuencia(sos, FS)
    assert mag[f <= FP].min() >= -RP - 0.01      # rizado en banda de paso
    assert mag[f >= FSTOP].max() <= -AS + 0.01   # atenuación en rechazo