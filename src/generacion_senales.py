import json
from pathlib import Path

import numpy as np
from scipy.io import wavfile
from scipy.signal import butter, chirp, get_window, sosfiltfilt, welch
from scipy.stats import kurtosis


def rms(x):
    return float(np.sqrt(np.mean(np.square(x, dtype=np.float64))))


# SNR en dB entre la señal limpia y la observada (ecuación 2 de la tarea).
def snr_db(clean, observed):
    error = observed - clean
    return 10 * np.log10((np.sum(clean**2) + 1e-15) / (np.sum(error**2) + 1e-15))


# Entrada y salida suaves (rampa seno^2) para evitar un "clic" al inicio y al final.
def fade(x, fs, ms=20):
    n = min(int(fs * ms / 1000), len(x) // 2)
    if n > 0:
        ramp = np.sin(np.linspace(0, np.pi / 2, n)) ** 2
        x = x.copy()
        x[:n] *= ramp
        x[-n:] *= ramp[::-1]
    return x


# Escala todas las señales con la misma ganancia para que el pico máximo sea `peak`.
# Al usar una ganancia común, las proporciones entre señales y la SNR no cambian.
def normalizar_conjunto(*signals, peak=0.8):
    maximum = max(np.max(np.abs(x)) for x in signals)
    gain = peak / maximum if maximum > 0 else 1.0
    return [x * gain for x in signals], gain


# Señal útil según `kind`: 'tono', 'multitono', 'chirp', 'armonicos' o 'transitoria'.
def generar_limpia(kind, t, freqs, amps, phases_deg, rng):
    if kind == "tono":
        return amps[0] * np.sin(2 * np.pi * freqs[0] * t + np.deg2rad(phases_deg[0]))
    if kind == "multitono":
        return sum(a * np.sin(2 * np.pi * f * t + np.deg2rad(p))
                   for f, a, p in zip(freqs, amps, phases_deg))
    if kind == "chirp":
        return chirp(t, f0=freqs[0], f1=freqs[-1], t1=t[-1], method="linear")
    if kind == "armonicos":
        f0 = freqs[0]
        return sum((1 / k) * np.sin(2 * np.pi * k * f0 * t + rng.uniform(-np.pi, np.pi))
                   for k in range(1, 6))
    if kind == "transitoria":
        env = np.exp(-5 * (t % 1.0)) * (np.sin(2 * np.pi * 2 * t) > 0)
        return env * np.sin(2 * np.pi * freqs[0] * t)
    raise ValueError(f"Tipo de señal no reconocido: {kind}")


# Ruido blanco filtrado con un pasa banda Butterworth (orden 6) entre `band[0]` y `band[1]` Hz.
def ruido_banda(n, fs, band, rng):
    low, high = band
    if not (0 < low < high < fs / 2):
        raise ValueError(f"BANDA_RUIDO debe cumplir 0 < low < high < {fs / 2}")
    sos = butter(6, [low, high], btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, rng.standard_normal(n))


def generar_ruido(kind, t, fs, f_int, band, impulses, rng):
    n = len(t)
    tonal = np.sin(2 * np.pi * f_int * t + rng.uniform(-np.pi, np.pi))
    red = sum((1 / k) * np.sin(2 * np.pi * 60 * k * t + rng.uniform(-np.pi, np.pi))
              for k in range(1, 5))
    white = rng.standard_normal(n)
    band_noise = ruido_banda(n, fs, band, rng)
    # Impulsos: pulsos cortos con forma de ventana Hann, de signo aleatorio
    impulse = np.zeros(n)
    locations = rng.choice(np.arange(int(.1 * fs), n - int(.1 * fs)),
                           size=min(impulses, max(1, n // 100)), replace=False)
    for loc in locations:
        width = max(2, int(rng.uniform(.0005, .004) * fs))
        stop = min(n, loc + width)
        impulse[loc:stop] += rng.choice([-1, 1]) * np.hanning(2 * (stop - loc))[:stop - loc]
    options = {"tonal": tonal, "red_60hz": red, "blanco": white,
               "banda": band_noise, "impulsivo": impulse}
    if kind == "mixto":
        return 0.55 * tonal + 0.30 * band_noise + 0.15 * impulse
    if kind not in options:
        raise ValueError(f"Tipo de ruido no reconocido: {kind}")
    return options[kind]


# Escala el ruido para que la SNR (limpia vs. ruido) sea exactamente `target_db`.
def escalar_ruido_para_snr(clean, noise, target_db):
    target_noise_rms = rms(clean) / (10 ** (target_db / 20))
    return noise * (target_noise_rms / (rms(noise) + 1e-15))


# Magnitud del espectro en dB, con ventana Hann; devuelve (frecuencias en Hz, magnitud en dB).
def espectro(x, fs):
    w = get_window("hann", len(x), fftbins=True)
    X = np.fft.rfft(x * w)
    f = np.fft.rfftfreq(len(x), 1 / fs)
    mag = 20 * np.log10(np.maximum(np.abs(X) / (np.sum(w) / 2), 1e-12))
    return f, mag


def fase_espectro(x, fs, umbral_rel=1e-2):
    X = np.fft.rfft(x)
    f = np.fft.rfftfreq(len(x), 1 / fs)
    fase = np.degrees(np.angle(X))
    fase[np.abs(X) < umbral_rel * np.max(np.abs(X))] = np.nan
    return f, fase


def energia_por_bandas(x, fs, bordes):

    X = np.fft.rfft(x)
    f = np.fft.rfftfreq(len(x), 1 / fs)
    p = np.abs(X) ** 2
    total = np.sum(p)
    return np.array([np.sum(p[(f >= lo) & (f < hi)]) / total
                     for lo, hi in zip(bordes[:-1], bordes[1:])])


def caracteristicas_ruido(x, fs):
    f, pxx = welch(x, fs=fs, nperseg=4096)
    # Piso numérico para no dividir ni tomar log de cero (limita el rango a 120 dB)
    pxx = pxx + 1e-12 * np.max(pxx)
    sfm = float(np.exp(np.mean(np.log(pxx))) / np.mean(pxx))
    picos_db = 10 * np.log10(np.max(pxx) / np.median(pxx))
    acumulada = np.cumsum(pxx) / np.sum(pxx)
    bw95 = float(f[np.searchsorted(acumulada, 0.975)] - f[np.searchsorted(acumulada, 0.025)])
    return {
        "curtosis_exceso": float(kurtosis(x)),
        "planicidad_espectral": sfm,
        "pico_sobre_mediana_dB": float(picos_db),
        "ancho_banda_95_Hz": bw95,
    }


def generar_caso(cfg):

    rng = np.random.default_rng(cfg["semilla"])
    fs = cfg["fs_Hz"]
    t = np.arange(int(fs * cfg["duracion_s"])) / fs
    # Verificación de Nyquist: fs debe ser al menos el doble de la frecuencia más alta
    f_max = max(cfg["frecuencias_utiles_Hz"]
                + [cfg["frecuencia_interferencia_Hz"], cfg["banda_ruido_Hz"][1]])
    if fs < 2 * f_max:
        raise ValueError("La configuración viola Nyquist.")
    clean = generar_limpia(cfg["tipo_senal"], t, cfg["frecuencias_utiles_Hz"],
                           cfg["amplitudes_utiles"], cfg["fases_grados"], rng)
    clean = fade(clean, fs, cfg["fade_ms"])
    noise_raw = generar_ruido(cfg["tipo_ruido"], t, fs, cfg["frecuencia_interferencia_Hz"],
                              tuple(cfg["banda_ruido_Hz"]), cfg["numero_impulsos"], rng)
    noise = escalar_ruido_para_snr(clean, noise_raw, cfg["snr_objetivo_dB"])
    noisy = clean + noise
    (clean, noise, noisy), gain = normalizar_conjunto(clean, noise, noisy, peak=cfg["pico_maximo"])
    return {"t": t, "clean": clean, "noise": noise, "noisy": noisy, "ganancia": gain,
            "snr_obtenida_dB": float(snr_db(clean, noisy))}


# Convierte una señal en [-1, 1] a enteros de 16 bits (formato PCM de los WAV).
def pcm16(x):
    return np.int16(np.clip(x, -1, 1) * 32767)


# Guarda en `carpeta` los tres WAV del caso (limpia, ruido, contaminada) y su JSON de metadatos.
def exportar_caso(cfg, res, carpeta):
    carpeta = Path(carpeta)
    carpeta.mkdir(parents=True, exist_ok=True)
    fs = cfg["fs_Hz"]
    wavfile.write(carpeta / "01_referencia_limpia.wav", fs, pcm16(res["clean"]))
    wavfile.write(carpeta / "02_ruido.wav", fs, pcm16(res["noise"]))
    wavfile.write(carpeta / "03_entrada_contaminada.wav", fs, pcm16(res["noisy"]))
    meta = dict(cfg, N=len(res["t"]), snr_obtenida_dB=res["snr_obtenida_dB"],
                ganancia_comun=float(res["ganancia"]))
    with open(carpeta / "metadatos_prueba.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    return meta
