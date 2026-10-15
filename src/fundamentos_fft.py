"""Fundamentos de la DFT/FFT (Tarea 1, sección 2.1, incisos c y d).

Tono en un bin y entre bins con ventana rectangular y Hann, lóbulos de las
ventanas y simetría conjugada de la DFT de una señal real.
"""
import numpy as np
import matplotlib.pyplot as plt

EPS = np.finfo(float).tiny


# Tono x[n] = sin(2*pi*f*n/fs + fase), amplitud 1 y sin ruido.
def generar_tono(frecuencia_hz, fs, n, fase_rad):
    return np.sin(2.0 * np.pi * frecuencia_hz * np.arange(n) / fs + fase_rad)


def ventana_hann_periodica(n):
    """Hann periódica, la adecuada para una DFT de N puntos: w[n] = 0.5 - 0.5*cos(2*pi*n/N).

    No es igual a np.hanning(N), que es simétrica y divide entre N-1.
    """
    return 0.5 - 0.5 * np.cos(2.0 * np.pi * np.arange(n) / n)


def calcular_espectro(x, w, fs):
    """Espectro de la mitad positiva (0 a fs/2) de la señal con ventana.

    Se divide entre la suma de la ventana para que un tono de amplitud 1 se vea con
    amplitud 1 con cualquier ventana, y se multiplica por 2 (salvo DC y Nyquist) para
    sumar la parte negativa. Esto solo corrige la escala; la fuga espectral se mantiene.
    """
    n = len(x)
    X = np.fft.rfft(x * w)
    f = np.fft.rfftfreq(n, d=1.0 / fs)
    amplitud = np.abs(X) / np.sum(w)
    if n % 2 == 0:
        amplitud[1:-1] *= 2
    else:
        amplitud[1:] *= 2
    return {"X": X, "f": f, "amplitud": amplitud, "fase": np.angle(X)}


def analizar_tonos(cfg):
    """Calcula los cuatro casos (tono en bin / entre bins × rectangular / Hann) y su tabla resumen."""
    n, fs = cfg["n"], cfg["fs_Hz"]
    fase = np.deg2rad(cfg["fase_grados"])
    ventanas = {"rectangular": np.ones(n), "hann": ventana_hann_periodica(n)}
    tonos = {"bin": cfg["f_bin_Hz"], "entre": cfg["f_entre_Hz"]}
    casos, tabla = {}, []
    for nombre_tono, frecuencia in tonos.items():
        x = generar_tono(frecuencia, fs, n, fase)
        for nombre_ventana, w in ventanas.items():
            r = calcular_espectro(x, w, fs)
            r["senal_original"] = x
            casos[(nombre_tono, nombre_ventana)] = r
            a = r["amplitud"]
            i_max = int(np.argmax(a[1:])) + 1
            # Fracción aproximada de la energía (mitad positiva, sin DC) fuera del bin más alto
            energia_total = np.sum(a[1:] ** 2)
            fuera = (energia_total - a[i_max] ** 2) / max(energia_total, EPS)
            tabla.append({
                "tono": nombre_tono, "frecuencia_tono_Hz": frecuencia, "ventana": nombre_ventana,
                "bin_dominante": i_max, "frecuencia_bin_dominante_Hz": float(r["f"][i_max]),
                "amplitud_bin_dominante": float(a[i_max]),
                "fase_bin_dominante_grados": float(np.rad2deg(r["fase"][i_max])),
                "fraccion_energia_fuera_bin_dominante": float(max(0.0, fuera)),
            })
    return casos, tabla


def trazar_bins(ax, datos, frecuencia_tono, fs, n, margen_bins=5, etiquetas=True, k_referencia=None):
    """Dibuja la DFT como puntos discretos (uno por bin), con zoom alrededor del tono."""
    delta_f = fs / n
    k_inferior = int(np.floor(frecuencia_tono / delta_f))
    k_izq = int(k_referencia) if k_referencia is not None else k_inferior
    k_der = k_izq + 1

    f, a = datos["f"], datos["amplitud"]
    inicio = max(0, min(k_izq, k_inferior) - margen_bins)
    final = min(len(f) - 1, max(k_der, k_inferior + 1) + margen_bins)
    ks = np.arange(inicio, final + 1)
    valores = a[ks]
    ybase, ytop = 0, 1.15 * max(1.0, float(np.max(valores)))
    ax.set_ylim(ybase, ytop)
    ax.set_ylabel("Amplitud normalizada")

    ax.vlines(f[ks], ybase, valores, color="#2666a7", alpha=.8, lw=1.6, zorder=2)
    ax.scatter(f[ks], valores, color="#2666a7", s=23, zorder=3, label="Bins de la DFT")
    ax.axvspan(k_izq * delta_f, k_der * delta_f, color="#f6b346", alpha=.12, zorder=0)
    for k, color in [(k_izq, "#d27a0c"), (k_der, "#238b70")]:
        ax.scatter([k * delta_f], [a[k]], color=color, s=71, zorder=5, edgecolors="#202020",
                   linewidths=.5, label=f"Bin k={k} ({k * delta_f:g} Hz)")
    ax.axvline(frecuencia_tono, color="#c73535", ls="--", lw=2.1, zorder=4,
               label=f"Tono real: {frecuencia_tono:g} Hz")
    ax.set_xlim(f[inicio] - delta_f * .4, f[final] + delta_f * .4)
    ax.set_xticks(ks * delta_f)
    ax.set_xticklabels([f"{k * delta_f:g}" for k in ks], fontsize=8)
    ax.grid(alpha=.22, axis="y")
    ax.set_xlabel("Frecuencia (Hz); cada punto es un bin separado Δf = fs/N")
    if etiquetas:
        ax.legend(loc="upper right", fontsize=8, ncol=2)


# Figura 2x2 de magnitud lineal: tono en bin / entre bins × rectangular / Hann.
def figura_comparacion_bins(casos, cfg):
    fs, n, f_bin, f_entre = cfg["fs_Hz"], cfg["n"], cfg["f_bin_Hz"], cfg["f_entre_Hz"]
    fig, axs = plt.subplots(2, 2, figsize=(13, 8.2), layout="constrained")
    for i, (nombre_tono, ftono) in enumerate((("bin", f_bin), ("entre", f_entre))):
        for j, (nombre_ventana, titulo) in enumerate((("rectangular", "Rectangular"), ("hann", "Hann periódica"))):
            trazar_bins(axs[i, j], casos[(nombre_tono, nombre_ventana)], ftono, fs, n,
                        etiquetas=(i == 0 and j == 1), k_referencia=int(round(f_bin / (fs / n))))
            tipo = "bin exacto" if i == 0 else "entre bins"
            axs[i, j].set_title(f"{ftono:g} Hz = {tipo} | {titulo}", fontsize=11, fontweight="bold")
    fig.suptitle("2.1(c) | Magnitud de la DFT: tono alineado frente a tono entre bins", fontweight="bold")
    return fig


def espectros_2x2(casos, campo, cfg):
    """Figura 2x2 en dB (`campo='db'`, fuga espectral) o de fase (`campo='fase'`)."""
    fs, n, fbin, fent = cfg["fs_Hz"], cfg["n"], cfg["f_bin_Hz"], cfg["f_entre_Hz"]
    df = fs / n
    nombres_ventanas = ("Rectangular", "Hann periódica")
    nombres_tonos = (f"Alineado: {fbin:g} Hz", f"Entre bins: {fent:g} Hz")
    fig, ax = plt.subplots(2, 2, figsize=(13, 8.4), sharex=True)
    fig.suptitle("2.1(c) | Magnitud relativa en dB (fuga espectral)" if campo == "db"
                 else "2.1(c) | Fase de la DFT donde la magnitud es significativa")
    for i, tipo_tono in enumerate(("bin", "entre")):
        for j, ventana in enumerate(("rectangular", "hann")):
            s = casos[(tipo_tono, ventana)]
            f, magnitud = s["f"], s["amplitud"]
            axis = ax[i, j]
            if campo == "db":
                rel_db = 20 * np.log10(np.maximum(magnitud / max(magnitud.max(), EPS), 1e-10))
                axis.stem(f, rel_db, linefmt="C0-", markerfmt="C0.", basefmt=" ")
                axis.set_ylim(-95, 4)
                axis.set_ylabel("Magnitud relativa (dB)")
            else:
                # Donde la magnitud es casi 0 la fase es un ángulo arbitrario de redondeo; se omite
                significativa = magnitud > 0.015 * np.max(magnitud)
                axis.scatter(f[significativa], np.rad2deg(s["fase"][significativa]), s=13)
                axis.set_ylim(-195, 195)
                axis.set_yticks([-180, -90, 0, 90, 180])
                axis.set_ylabel("Fase (grados)")
                axis.text(0.02, 0.03, "Solo bins >1.5 % del pico", transform=axis.transAxes, fontsize=8)
            axis.set_title(f"{nombres_tonos[i]} | {nombres_ventanas[j]}")
            axis.set_xlim(max(0, min(fbin, fent) - 10 * df), max(fbin, fent) + 10 * df)
            axis.axvline(fbin if tipo_tono == "bin" else fent, ls="--", color="C3", linewidth=1.3)
            if tipo_tono == "entre":
                abajo = int(np.floor(fent / df)) * df
                axis.axvline(abajo, ls=":", color="0.5", linewidth=1)
                axis.axvline(abajo + df, ls=":", color="0.5", linewidth=1)
            axis.grid(alpha=0.25)
            axis.set_xlabel("Frecuencia (Hz)")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def medir_ventanas(n, fs):
    """Mide el ancho del lóbulo principal y el nivel del primer lóbulo lateral de cada ventana.

    Se rellena con ceros (zero-padding) solo para ver la forma de la ventana con detalle;
    eso interpola la curva, no mejora la resolución real.
    """
    df = fs / n
    ventanas = {"rectangular": np.ones(n), "hann": ventana_hann_periodica(n)}
    sobremuestreo = 256
    nfft = n * sobremuestreo
    frecuencias = np.fft.rfftfreq(nfft, d=1.0 / fs)
    fig, axis = plt.subplots(figsize=(11, 5.7))
    medidas = {}
    for nombre, ventana in ventanas.items():
        respuesta = np.abs(np.fft.rfft(ventana, n=nfft)) / np.sum(ventana)
        db = 20 * np.log10(np.maximum(respuesta, 1e-12))
        # Distancia del centro al primer cero: rectangular 1 bin, Hann 2 bins
        semi_ancho_bins = 1 if nombre == "rectangular" else 2
        ancho_hz = 2.0 * semi_ancho_bins * df
        # Primer máximo local después del primer cero
        idx_inicio = semi_ancho_bins * sobremuestreo + 1
        idx_fin = min(len(db) - 2, int(16 * sobremuestreo))
        maximos = np.where((db[idx_inicio:idx_fin] > db[idx_inicio - 1:idx_fin - 1]) &
                           (db[idx_inicio:idx_fin] >= db[idx_inicio + 1:idx_fin + 1]))[0]
        if len(maximos) == 0:
            raise RuntimeError("No se pudo encontrar el primer lóbulo lateral")
        primer_pico = idx_inicio + int(maximos[0])
        pico_lateral_db = float(db[primer_pico])
        medidas[nombre] = {
            "ancho_lobulo_principal_Hz": ancho_hz,
            "ancho_lobulo_principal_bins": 2 * semi_ancho_bins,
            "primer_lobulo_lateral_dB": round(pico_lateral_db, 3),
            "frecuencia_primer_lobulo_lateral_Hz": float(frecuencias[primer_pico]),
        }
        axis.plot(frecuencias / df, db, linewidth=1.7,
                  label=f"{nombre.capitalize()}: ancho {ancho_hz:g} Hz; 1er lateral {pico_lateral_db:.1f} dB")
        axis.axvline(semi_ancho_bins, color="0.6", ls=":", lw=0.9)
    axis.set(xlim=(0, 15), ylim=(-100, 4), xlabel="Frecuencia relativa al centro (en bins de Δf)",
             ylabel="Magnitud de respuesta de ventana (dB)",
             title="2.1(c) | Lóbulo principal y lóbulos laterales de las ventanas")
    axis.grid(alpha=0.3)
    axis.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    return fig, medidas


def validar_conjugada(X, escala_referencia=None):
    """Comprueba X[N-k] = conj(X[k]); devuelve el error máximo absoluto y el relativo.

    `escala_referencia` evita dividir entre casi cero cuando el espectro editado queda casi vacío.
    """
    n = len(X)
    k = np.arange(n)
    error = float(np.max(np.abs(X[(-k) % n] - np.conj(X[k]))))
    escala = float(np.max(np.abs(X))) if escala_referencia is None else escala_referencia
    return error, error / max(escala, EPS)


# Fase de los dos bins del par conjugado X[kp] y X[kn]: misma magnitud, fases opuestas.
def graficar_fase_conjugada(ax, X, kp, kn):
    indices = np.array([kp, kn])
    fases = np.rad2deg(np.angle(X[indices]))
    for k, fase, color in zip(indices, fases, ("#d27a0c", "#238b70")):
        ax.vlines(k, 0, fase, color=color, lw=2.4)
        ax.scatter([k], [fase], color=color, s=64, zorder=3)
        ax.annotate(f"X[{k}]: {fase:+.1f}°", (k, fase),
                    xytext=(10 if k == kp else -10, 10 if fase >= 0 else -15), textcoords="offset points",
                    ha="left" if k == kp else "right", fontsize=9, fontweight="bold", color=color)
    ax.axhline(0, color="0.45", lw=0.8)
    ax.set(xlim=(-8, len(X) + 7), ylim=(-100, 100), title="Fase del espectro: pares conjugados",
           xlabel="Bin k", ylabel="Fase (grados)")
    ax.set_yticks([-90, -45, 0, 45, 90])
    ax.grid(alpha=.24)


# Texto con los valores del par conjugado y la parte imaginaria de cada IFFT.
def mostrar_resumen_conjugado(ax, X, kp, kn, im_malo, im_bueno):
    ax.axis("off")
    ax.set_title("Comprobación numérica", fontweight="bold")
    mas, menos = X[kp], X[kn]
    fase_mas, fase_menos = np.rad2deg(np.angle([mas, menos]))
    texto = (f"X[{kp}] = {mas.real:.2f} {mas.imag:+.2f}j\n"
             f"X[{kn}] = {menos.real:.2f} {menos.imag:+.2f}j\n\n"
             f"|X[{kp}]| = |X[{kn}]| = {abs(mas):.1f}\n"
             f"Fases: {fase_mas:+.1f}° y {fase_menos:+.1f}°\n\n"
             f"Im(IFFT) al editar solo un bin: {im_malo:.3f}\n"
             f"Im(IFFT) al editar ambos: {im_bueno:.2e}")
    ax.text(.06, .93, texto, transform=ax.transAxes, va="top", ha="left",
            fontsize=10, linespacing=1.55, family="monospace")


def simetria_y_ifft(casos, cfg):
    """Verifica la simetría conjugada en los 4 casos y compara editar el espectro mal y bien.

    Mal: anular solo el bin positivo del tono. Bien: anular también su pareja conjugada.
    Devuelve la figura y un diccionario con los errores y la parte imaginaria de cada IFFT.
    """
    n, fs, fbin = cfg["n"], cfg["fs_Hz"], cfg["f_bin_Hz"]
    errores = {}
    for (tono, ventana), datos in casos.items():
        # rfft solo da la mitad positiva; para ver la pareja conjugada hace falta la fft completa
        w = np.ones(n) if ventana == "rectangular" else ventana_hann_periodica(n)
        e_abs, e_rel = validar_conjugada(np.fft.fft(datos["senal_original"] * w))
        errores[f"{tono}_{ventana}"] = {"error_absoluto": e_abs, "error_relativo": e_rel}
        assert e_rel < 1e-12, "Una señal real debe tener una DFT con simetría conjugada"

    xbin = casos[("bin", "rectangular")]["senal_original"]
    X_original = np.fft.fft(xbin)
    kp = int(round(fbin / (fs / n)))
    kn = (n - kp) % n
    assert kp != kn and 0 < kp < n // 2

    X_incorrecto = X_original.copy()
    X_incorrecto[kp] = 0.0
    X_correcto = X_incorrecto.copy()
    X_correcto[kn] = 0.0

    y_original, y_incorrecto, y_correcto = (np.fft.ifft(X) for X in (X_original, X_incorrecto, X_correcto))
    escala = float(np.max(np.abs(X_original)))
    antes_abs, antes_rel = validar_conjugada(X_original, escala)
    malo_abs, malo_rel = validar_conjugada(X_incorrecto, escala)
    bueno_abs, bueno_rel = validar_conjugada(X_correcto, escala)
    im_malo = float(np.max(np.abs(np.imag(y_incorrecto))))
    im_bueno = float(np.max(np.abs(np.imag(y_correcto))))
    assert antes_rel < 1e-12 and bueno_rel < 1e-12
    assert malo_rel > 1e-2 and im_malo > 1e-2 and im_bueno < 1e-12

    tiempo = np.arange(n) / fs
    vis = slice(0, min(n, 90))
    fig, axs = plt.subplots(3, 2, figsize=(13, 10.6), layout="constrained")
    axs[0, 0].stem(np.arange(n), np.abs(X_original), linefmt="C0-", markerfmt="C0.", basefmt=" ")
    axs[0, 0].set(title=f"DFT completa, picos conjugados k={kp} y k={kn}",
                  xlabel="Índice de bin k", ylabel="|X[k]| (sin normalizar)")
    graficar_fase_conjugada(axs[0, 1], X_original, kp, kn)

    axs[1, 0].plot(np.arange(n), np.abs(X_incorrecto), ".-", markersize=2, label="Incorrecto: anula solo +f")
    axs[1, 0].plot(np.arange(n), np.abs(X_correcto), "--", linewidth=1.1, label="Correcto: anula +f y -f")
    axs[1, 0].set(title="Edición espectral: sin y con pareja conjugada", xlabel="Índice de bin k",
                  ylabel="|X_modificado[k]|")
    axs[1, 0].legend(fontsize=8)
    mostrar_resumen_conjugado(axs[1, 1], X_original, kp, kn, im_malo, im_bueno)

    axs[2, 0].plot(tiempo[vis], np.real(y_original[vis]), label="Original real")
    axs[2, 0].plot(tiempo[vis], np.real(y_incorrecto[vis]), label="Real de IFFT incorrecta", ls="--")
    axs[2, 0].plot(tiempo[vis], np.real(y_correcto[vis]), label="Real de IFFT corregida", ls=":")
    axs[2, 0].set(title="Partes reales reconstruidas", xlabel="Tiempo (s)", ylabel="Amplitud")
    axs[2, 0].legend(fontsize=8)
    axs[2, 1].plot(tiempo[vis], np.imag(y_incorrecto[vis]), label=f"Incorrecto (máx={im_malo:.3f})")
    axs[2, 1].plot(tiempo[vis], np.imag(y_correcto[vis]), label=f"Correcto (máx={im_bueno:.1e})", ls="--")
    axs[2, 1].set(title="La IFFT debe ser real si se preserva la simetría", xlabel="Tiempo (s)",
                  ylabel="Parte imaginaria")
    axs[2, 1].legend(fontsize=8)
    for ax in (axs[0, 0], axs[1, 0], axs[2, 0], axs[2, 1]):
        ax.grid(alpha=0.2)
    fig.suptitle("2.1(d) | Simetría conjugada: X[N-k] = conj(X[k])")

    resultado = {
        "errores_4_casos": errores,
        "k_positivo": kp, "k_negativo": kn,
        "magnitud_bin_positivo": float(abs(X_original[kp])), "magnitud_bin_negativo": float(abs(X_original[kn])),
        "fase_bin_positivo_grados": float(np.rad2deg(np.angle(X_original[kp]))),
        "fase_bin_negativo_grados": float(np.rad2deg(np.angle(X_original[kn]))),
        "error_simetria_original": antes_abs, "error_simetria_edicion_incorrecta": malo_abs,
        "error_simetria_edicion_correcta": bueno_abs,
        "max_imaginaria_edicion_incorrecta": im_malo, "max_imaginaria_edicion_correcta": im_bueno,
    }
    return fig, resultado
