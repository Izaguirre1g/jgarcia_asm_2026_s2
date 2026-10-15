# Informe: Tarea 1, caracterización espectral y prototipado de filtros

Borrador en Markdown; la entrega final es `informe.pdf`.

## Sección 2.1: Fundamentos, muestreo y FFT

Pendiente.

## Sección 2.2: Generación y caracterización de señales

Parámetros de los tres casos: notebook, sección 1, y `metadatos/configuraciones.json`.

### Fase del espectro

en el bin de cada tono la fase de la señal limpia es la fase configurada menos 90° (el seno se mide respecto al coseno de la FFT): -90°, -55° y -140° para 440, 880 y 1320 Hz. Los puntos alrededor de cada tono son bins vecinos con fuga espectral, causada por el fade de 20 ms, y están desfasados 180° respecto al bin central.

- **Caso 1:** los tres tonos conservan su fase. La interferencia aparece como un único punto nuevo en 1750 Hz.
- **Caso 2:** el ruido blanco rellena todo el espectro con fases aleatorias uniformes en ±180°. Los tonos siguen dominando en su bin, por eso conservan la fase.
- **Caso 3:** la fase aleatoria queda confinada a 700-1500 Hz. En los bins exactos de 880 y 1320 Hz la fase del tono sigue visible, porque la energía del tono se concentra en un bin; pero los bins vecinos, con fase distorsionada, ya no son separables del ruido.

### Características para detectar automáticamente el ruido

**Interpretación (sobre el ruido aislado, que es la referencia, y la señal contaminada, que es lo que ve el sistema):**

| Tipo | Curtosis en exceso | Planicidad espectral | Ancho de banda 95 % |
|---|---|---|---|
| Tonal (caso 1) | -1.5 (sinusoide) | ≈ 0 | ≈ 8 Hz |
| Blanco (caso 2) | ≈ 0 (gaussiano) | 0.99 | ≈ 7600 Hz |
| Banda (caso 3) | ≈ 0 (gaussiano) | ≈ 0 | ≈ 740 Hz |

- **Interferencia tonal:** la curtosis de -1.5 es propia de una sinusoide y su ancho de banda es de unos pocos Hz: es un pico aislado fuera de los tonos útiles.
- **Ruido blanco:** es el único con planicidad cercana a 1 y ancho de banda casi igual al disponible. La curtosis no lo distingue del ruido de banda, porque ambos son gaussianos. En la señal contaminada el piso de ruido queda elevado (pico sobre mediana de 34 dB, frente a valores saturados en 120 dB de los otros dos casos).
- **Ruido de banda limitada:** gaussiano pero con ancho de banda acotado (≈ 740 Hz), lo que lo separa del blanco. La planicidad total es baja porque la mayor parte del espectro está vacía.
- **En el sistema real** no hay referencia limpia, así que estos descriptores se estiman sobre bloques de la señal contaminada, donde los tonos útiles también aportan picos. Por eso conviene medir el piso de ruido entre los picos y no solo el máximo.

### ¿Qué no se puede recuperar cuando señal y ruido comparten frecuencias? (caso 3)

En 700-1500 Hz hay tonos útiles (880 y 1320 Hz) y ruido. Un filtro o una máscara espectral solo puede actuar sobre la frecuencia, no distingue qué parte de la energía de un bin es señal y cuál es ruido. Por eso:

1. **Información irrecuperable:** en los bins donde coexisten señal y ruido, el ruido es una suma aditiva con fase aleatoria. Atenuar la banda para quitar ruido también atenúa y distorsiona amplitud y fase de los tonos de 880 y 1320 Hz; dejarla intacta conserva el ruido. No hay forma de separar ambos sin información adicional (modelo de la señal, referencia, etc.).
2. **Compromiso:** cuanto más energía del ruido se elimina, más energía útil se pierde (mayor MSE respecto a la señal limpia). Esto se mide en la sección 2.5.
3. **Qué sí se puede:** estimar y restar el nivel medio de ruido fuera de los picos de los tonos (por ejemplo, sustracción espectral), reduciendo el ruido pero sin eliminar la componente que cae en el mismo bin del tono.
4. **Caso contrario (casos 1 y 2):** con una interferencia separada se puede rechazar solo esa frecuencia sin tocar la señal útil; con ruido blanco, el ruido que cae fuera de los picos puede quitarse, pero el que cae bajo los tonos permanece.

## Sección 2.3: Diseño del filtro temporal

**Caso de diseño: caso 2 (ruido blanco).** El ruido blanco se reparte de 0 a 8 kHz y toda la señal útil está por debajo de 1320 Hz, por lo que un **pasa bajas** corresponde con el espectro observado (sección 2.2).

| Parámetro | Valor | Justificación |
|---|---|---|
| Banda de paso | 0 – 1500 Hz | Contiene 440, 880 y 1320 Hz con 180 Hz de margen |
| Banda de rechazo | ≥ 2000 Hz | Sin señal útil por encima de 1320 Hz |
| Rizado en banda de paso | ≤ 1 dB | Distorsión de amplitud tolerable en los tonos |
| Atenuación en rechazo | ≥ 40 dB | Deja el ruido fuera de banda >100 veces menor en potencia |
| fs | 16 kHz | Igual que las señales |

Se comparan tres prototipos IIR en secciones de segundo orden (SOS): Butterworth, Chebyshev I y elíptico.

**Nota metodológica sobre la SNR.** Un filtro causal retrasa la señal, y la ecuación (2) compara muestra a muestra contra la referencia sin retraso, así que el retardo cuenta como error aunque no sea ruido. Se reportan dos valores: **ΔSNR sin alinear** (ecuación literal) y **ΔSNR alineada**, que compensa el retardo constante (entre 0 y 100 muestras) que maximiza la SNR. Como el retardo de un IIR depende de la frecuencia, la alineación no elimina todo el error de fase. Como referencia se incluye también la salida de **fase cero** (`sosfiltfilt`), que no es causal y no es implementable en tiempo real, pero separa el efecto de la magnitud del efecto de la fase.

### Justificación del tipo y orden 

Para la misma especificación (notebook, sección 8), el orden mínimo es **17** (Butterworth), **8** (Chebyshev I) y **5** (elíptico). Con K secciones, el costo es de 5K coeficientes, 2K estados y 5K multiplicaciones por muestra: 45, 20 y 15 MAC por muestra, respectivamente. El elíptico necesita la tercera parte de los recursos del Butterworth, a cambio de rizado en ambas bandas.

**Evidencia medida en el caso 2** (notebook, sección 8):

- **Especificación:** los tres cumplen. La ganancia en 440, 880 y 1320 Hz está entre 0 y −1 dB, y la atenuación es de al menos 40 dB desde 2000 Hz (figura `2_3_comparacion_prototipos.png`).
- **ΔSNR sin alinear** (ecuación 3 literal): −9.99 dB (Butterworth), −9.43 dB (Chebyshev I) y −6.14 dB (elíptico). Es negativo porque todo filtro causal retrasa la señal respecto a la referencia sin retraso, y ese desfase se contabiliza como error. El MSE del elíptico (0.059) es menos de la mitad que el de los otros dos (0.126 y 0.144).
- **ΔSNR alineado** por el retardo constante que maximiza la SNR (18, 13 y 6 muestras, que coinciden con el retardo de grupo en 440 Hz, el tono dominante): +4.56, +5.14 y +5.34 dB. Sigue por debajo de la referencia de fase cero porque el retardo de un IIR varía con la frecuencia, y un retardo constante solo alinea un tono; los demás quedan desfasados.
- **Referencia de fase cero** (`sosfiltfilt`, no causal, no implementable en tiempo real): +7.27, +7.55 y +6.35 dB, de acuerdo con el valor teórico para ruido blanco, 10·log10(8000/1500..1600) ≈ 7.0-7.3 dB. El elíptico queda ~1 dB por debajo porque su rizado (hasta 1 dB) distorsiona la amplitud de los tonos y `sosfiltfilt` eleva al cuadrado la respuesta.
- **Retardo de grupo:** en 440 Hz, 17.8, 12.3 y 4.7 muestras; en 1320 Hz, 27.3, 23.4 y 11.7 muestras. El pico en el borde de la banda de paso es de 39, 59 y 40 muestras.

**Decisión: filtro elíptico de orden 5.** Tiene el menor costo, el menor MSE, el mayor ΔSNR alineado y el menor retardo en los tonos útiles. Su desventaja medida es el rizado de hasta 1 dB en los tonos y la menor ganancia en la referencia de fase cero. La decisión se basa en estas medidas y no en las propiedades generales de cada prototipo.

### Validación de la ecuación de diferencias

La ecuación de diferencias se implementó a mano (`aplicar_sos_manual`, forma directa II transpuesta por sección) y se comparó con `sosfilt` y con `lfilter(b, a)`, con una **tolerancia declarada de 1e-10**: el `assert` del notebook (sección 8) falla si se excede. Las entradas fueron ruido gaussiano de 96 000 muestras (semilla 0) y la señal contaminada del caso 2.

| Entrada | Error máx. vs `sosfilt` | Error máx. vs `lfilter(b, a)` |
|---|---|---|
| Ruido gaussiano | 0.0 | 1.4e-13 |
| Caso 2 contaminada | 0.0 | 4.7e-14 |

El error frente a `sosfilt` es nulo porque ambas implementan la misma estructura, con las mismas operaciones en float64 y en el mismo orden. Frente a `lfilter(b, a)`, que usa un único polinomio de grado 5, la diferencia es solo de redondeo (≈1e-13). Esa forma se vuelve peor condicionada al aumentar el orden, lo que justifica implementar el filtro como cascada de secciones de segundo orden.

### Aplicación del filtro a los tres casos

El filtro elíptico, diseñado para el caso 2, se aplicó a los tres casos (`resultados/metricas.csv`, figura `2_3_aplicacion_tres_casos.png`):

| Caso | ΔSNR sin alinear | ΔSNR alineado | Energía útil conservada | Observación |
|---|---|---|---|---|
| 1 | −5.97 dB | +9.02 dB | 84.2 % | Interferencia de 1750 Hz atenuada 30.4 dB |
| 2 | −6.14 dB | +5.34 dB | 84.8 % | Caso de diseño |
| 3 | −6.85 dB | −0.08 dB | 86.2 % | Sin mejora |

- **Caso 1:** la interferencia de 1750 Hz cae en la banda de transición (1500-2000 Hz), cerca del cero de transmisión en 1850 Hz, y se atenúa 30.4 dB, menos que los 40 dB de la banda de rechazo. El filtro mejora la SNR aunque no se diseñó para este caso, porque la interferencia está por encima de los tonos útiles.
- **Caso 2:** el ruido por encima de 1700 Hz queda unos 90 dB o más por debajo de la señal. El ΔSNR alineado de +5.34 dB es menor que el valor teórico de ≈ +7 dB, por la distorsión que el filtro causa en los tonos: el rizado (hasta 1 dB) y la dispersión del retardo de grupo, que un retardo constante no corrige.
- **Caso 3:** el ruido ocupa 700-1500 Hz, que contiene los tonos de 880 y 1320 Hz. El filtro solo recorta la cola del ruido por encima de 1500 Hz y no mejora la SNR (−0.08 dB). Un filtro lineal no puede separar señal y ruido que comparten frecuencias (sección 2.2).
- **Energía útil conservada:** 84-86 % en los tres casos. Corresponde al rizado de la banda de paso: las ganancias en los tonos (−0.83, −0.42 y −1.00 dB), ponderadas por la potencia de cada tono, dan unos 84 %. La cifra del caso 3 es algo mayor porque probablemente el ruido dentro de ±20 Hz de cada tono suma energía a la medición.
- **Efecto del retardo:** el ΔSNR sin alinear es negativo en los tres casos porque un filtro causal retrasa la señal respecto a la referencia sin retraso, y ese desfase cuenta como error en la ecuación (3). Con la alineación por 6 muestras, que es el retardo de grupo del filtro en 440 Hz, el resultado es positivo en los casos 1 y 2.

### Cuantización, saturación y memoria

**Microcontrolador seleccionado:** ESP32. Para el procesamiento se propone la biblioteca ESP-DSP de Espressif, que incluye rutinas de filtros biquad en float32 y de FFT.

- **Memoria y cómputo:** el filtro tiene K = 3 secciones, es decir 5K = 15 coeficientes y 2K = 6 estados: **84 bytes** en float32. El costo es de 15 multiplicaciones por muestra, o 0.24 millones de MAC por segundo a 16 kHz, una carga despreciable frente a un procesador de 240 MHz.
- **Formato numérico propuesto: float32.** Como el ESP32 tiene FPU, el filtro se implementa en float32 y no se necesita aritmética de punto fijo. La precisión relativa de float32 (≈6e-8) equivale a más de 20 bits fraccionarios en los coeficientes, por encima de los 14 bits a partir de los cuales la tabla siguiente cumple la especificación. Los resultados de Python en float64 son, por tanto, representativos del filtro en el microcontrolador.
- **Cuantización de coeficientes** (simulada en Python, redondeando a múltiplos de 2⁻ᵇ y recalculando la respuesta). Sirve como referencia por si se usara punto fijo, o para una variante del ESP32 sin FPU:

| Bits fraccionarios | Desviación máx. en banda de paso | Atenuación mín. en rechazo | Estable |
|---|---|---|---|
| 15 | 0.998 dB | 39.995 dB | Sí |
| 14 | 0.997 dB | 39.996 dB | Sí |
| 12 | 0.957 dB | 39.959 dB | Sí |
| 10 | 1.293 dB | 40.263 dB | Sí |
| 8 | 0.776 dB | 38.885 dB | Sí |

  - Con 14 o más bits el filtro es prácticamente idéntico al original (a 0.005 dB de la especificación de 1 dB y 40 dB). El diseño elíptico toca los 40 dB exactos, por lo que cualquier redondeo lo deja unas milésimas por debajo.
  - Con 10 bits el rizado llega a 1.29 dB y excede la especificación; con 8 bits la atenuación baja a 38.9 dB y tampoco cumple.
  - En ningún caso el filtro se vuelve inestable: el radio máximo de los polos pasa de 0.973 a 0.972. La cascada de biquads mantiene los polos casi fijos.
  - Parte de la degradación con pocos bits viene de la ganancia global: el coeficiente b₀ = 0.0141 de la primera sección es muy pequeño y pierde precisión relativa. Con 10 bits, la ganancia cambia en −0.26 dB, y con 8 bits en +0.90 dB, y ese error se suma o se resta al rizado. Se mitigaría sacando la ganancia como un factor aparte.
  - Si se usara punto fijo de 16 bits, el formato sería Q2.14 (rango [−2, 2)), porque los coeficientes a₁ llegan a −1.62 y no caben en Q15 puro.
- **Saturación:** en float32 no hay desbordamiento en el sentido de los enteros, pero la salida debe limitarse al rango del DAC o del formato de salida (por ejemplo, [−1, 1] antes de convertir a 16 bits), para que un pico no produzca recorte con wrap-around. En una variante de punto fijo, los productos se acumularían en un registro más ancho que los datos (32 o 64 bits para entradas de 16) y con saturación aritmética. La simulación de esta sección cuantiza solo los coeficientes; el redondeo de los estados y de los productos no se midió.
- **Latencia:** el retardo de grupo del filtro es de 4.7 muestras en 440 Hz (0.30 ms) y 11.7 muestras en 1320 Hz (0.73 ms), con un pico de 40 muestras (2.5 ms) en el borde de la banda de paso.

## Sección 2.4: Reducción espectral mediante FFT e IFFT

**Cadena de procesamiento** (`src/filtrado_espectral.py`): bloque de N = 1024 muestras (potencia de 2; Δf = 15.6 Hz, 64 ms) → ventana raíz de Hann → `rfft` → ganancia por bin → `irfft` → ventana de síntesis → solapamiento y suma (overlap-add) con 50 % de solape.

- **Solapamiento:** con raíz de Hann en análisis y síntesis, el producto de ventanas es una Hann, que suma 1 con 50 % de solape. Sin modificar el espectro, la reconstrucción es exacta (error ~1e-15, verificado en `tests/`).
- **Simetría conjugada:** se usa `rfft`/`irfft`, que fuerzan $X[N-k] = X^*[k]$. La ganancia es real y por bin, así que la fase del espectro contaminado se conserva y no se introduce distorsión de fase.
- **Ganancia (resta espectral):** $G[k]^2 = \max\left(1 - \alpha\,\dfrac{P_n[k]}{P[k]},\; g_{min}^2\right)$, con $P[k]$ la potencia del bloque y $P_n[k]$ el piso de ruido estimado.
- **Parámetros adaptativos (estimados con los primeros 2 s):** (i) el piso de ruido $P_n[k]$ por bin, con la mediana en frecuencia de la PSD promedio (la mediana ignora los picos angostos de los tonos); (ii) las **frecuencias interferentes**: picos tonales que no pertenecen a la serie armónica de la fundamental (el pico más bajo). En esas frecuencias la ganancia se fija en $g_{min}$.
- **Parámetros de agresividad:** $\alpha = 2$ y $g_{min} = 0.1$.

**Supuestos y limitaciones:**
1. Los tonos útiles son armónicos de una fundamental. El 1750 Hz del caso 1 se distingue porque queda a 10 Hz del cuarto armónico (1760 Hz), con una tolerancia de 5 Hz; una interferencia ubicada sobre un armónico no sería detectada.
2. El ruido es estacionario durante la calibración.
3. Con ruido que comparte frecuencia con los tonos (caso 3) no hay separación posible por ganancia espectral.

### Resultados (α = 2, N = 1024)

| Caso | Interferencia estimada | ΔSNR | SNR de salida | Energía útil conservada |
|---|---|---|---|---|
| 1 | 1750.0 Hz | +18.87 dB | 23.87 dB | 99.8 % |
| 2 | ninguna | +8.77 dB | 13.77 dB | 99.2 % |
| 3 | ninguna | +7.00 dB | 12.00 dB | 95.6 % |

- **Caso 1:** el método detectó la interferencia sin conocer su frecuencia previamente (1750.0 Hz) y la atenuó 19.7 dB, que corresponde al límite impuesto por $g_{min} = 0.1$ (−20 dB). La energía útil se conserva casi intacta (99.8 %). En la figura `2_4_magnitud.png` se ve que el procesamiento deja componentes residuales de entre −55 y −100 dB alrededor de los tonos y de la interferencia; son artefactos de variar la ganancia bloque a bloque, con un nivel despreciable frente a los tonos.
- **Caso 2:** el piso de ruido baja unos 12-15 dB en todo el espectro (de ≈ −60 dB a ≈ −72 dB) y los tres tonos se conservan (99.2 %). El método no elimina el ruido que queda bajo los tonos y en su vecindad inmediata, donde se aprecian pequeños ensanchamientos del piso.
- **Caso 3:** el ruido de banda se reduce unos 15 dB dentro de 700-1500 Hz, pero no desaparece, y los tonos de 880 y 1320 Hz quedan rodeados de ruido residual. La energía útil conservada es la menor (95.6 %), porque en esa banda la ganancia afecta tanto al ruido como a los tonos (sección 2.2). Aun así, el método mejora la SNR en +7.00 dB, mientras que el filtro temporal de la sección 2.3 no logró mejora alguna en este caso (ΔSNR alineado de −0.08 dB).
- **Adaptación al tipo de ruido:** la figura `2_4_piso_ruido.png` muestra que el piso estimado $P_n[k]$ es plano en el caso 2 (≈ −44.5 dB), se concentra en 700-1500 Hz en el caso 3 (≈ −34 dB) y cae por encima de 1500 Hz. En el caso 1 no hay ruido de banda ancha, y el piso estimado (entre −50 y −120 dB) corresponde a la fuga espectral de los tonos, por lo que la resta espectral casi no actúa y la mejora proviene del rechazo de la interferencia detectada.
- **Fase y retardo:** la ganancia es real, así que el procesamiento no desfasa los tonos. El retardo de alineación es de 0 muestras y el ΔSNR sin alinear coincide con el alineado en los tres casos, a diferencia del filtro temporal, que necesitó compensar 6 muestras.

### Compromiso entre ruido eliminado y distorsión 

La figura `2_4_barrido_alpha.png` muestra, para el caso 3, cómo cambian la SNR de salida y la energía útil conservada al variar la agresividad $\alpha$:

- La SNR **sube** hasta un máximo de **13.13 dB en $\alpha = 3.5$** (frente a 12.00 dB en el $\alpha = 2$ usado en las tablas) y luego **baja** de forma sostenida hasta **8.16 dB en $\alpha = 16$**, es decir 5 dB menos que en el máximo.
- La energía útil conservada **solo disminuye** al aumentar $\alpha$: es ≈ 100 % con $\alpha$ pequeño (en esa zona la medida supera ligeramente 100 % porque el ruido dentro de ±20 Hz de cada tono suma energía), 95.6 % en $\alpha = 2$ y **74.1 % en $\alpha = 16$**.
- **Interpretación:** al principio, aumentar $\alpha$ elimina ruido que no está bajo los tonos y la SNR mejora. Pasado el óptimo, la ganancia empieza a recortar los tonos útiles, cuyo bin está contaminado por el ruido de banda, y el error por distorsión supera el beneficio de tener menos ruido. Esto es una consecuencia directa del solapamiento de señal y ruido en frecuencia (sección 2.2).
- Incluso con $\alpha = 16$ la SNR (8.16 dB) sigue por encima de la de entrada (5 dB); lo que se pierde por agresividad excesiva es la mejora, no la señal completa.
- **Elección de $\alpha = 2$:** es un punto intermedio que conserva el 95.6 % de la energía útil; con $\alpha = 3.5$ se ganan 1.1 dB de SNR a costa de unos 4 puntos porcentuales más de pérdida de energía útil (valor leído de la figura). Un valor de $\alpha$ distinto por caso, o adaptado al nivel de ruido, es una línea de mejora.

## Sección 2.5: Comparación y propuesta para el sistema embebido

### Comparación entre señal sin filtrar, filtro temporal y reducción espectral

Resultados medidos con los tres casos de prueba (SNR de entrada de 5 dB). Fuente: `resultados/comparacion_2_5.csv`.

| Caso | Método | ΔSNR sin alinear (dB) | ΔSNR alineada (dB) | Retardo (muestras) | Energía útil conservada (%) |
|---|---|---|---|---|---|
| 1 | Sin filtrar | 0.00 | 0.00 | 0 | 100.0 |
| 1 | Temporal (elíptico) | -5.97 | +9.02 | 6 | 84.2 |
| 1 | Espectral | +18.87 | +18.87 | 0 | 99.8 |
| 2 | Sin filtrar | 0.00 | 0.00 | 0 | 100.8 |
| 2 | Temporal (elíptico) | -6.14 | +5.34 | 6 | 84.8 |
| 2 | Espectral | +8.77 | +8.77 | 0 | 99.2 |
| 3 | Sin filtrar | 0.00 | 0.00 | 0 | 102.5 |
| 3 | Temporal (elíptico) | -6.85 | -0.08 | 6 | 86.2 |
| 3 | Espectral | +7.00 | +7.00 | 0 | 95.6 |

### ΔSNR sin alinear y alineada

La ecuación (3) compara la salida del sistema con la referencia limpia muestra a muestra. Un filtro causal retrasa la señal, y ese retraso se contabiliza como error aunque no sea ruido. Por ejemplo, un tono de 440 Hz que sale retrasado 6 muestras sigue siendo un tono perfecto, pero al restarlo de la referencia sin retraso la diferencia no es cero.

Por eso se reportan dos valores:

- **ΔSNR sin alinear:** aplica la ecuación (3) tal como está escrita. En el filtro temporal es negativo en los tres casos (entre -5.97 y -6.85 dB), es decir, la salida se aleja más de la referencia que la señal sin filtrar, pero la causa principal es el retardo y no una falta de filtrado.
- **ΔSNR alineada:** antes de medir, se desplaza la salida el retardo constante que maximiza la SNR (6 muestras). Mide cuánto ruido se quitó realmente. El retardo de un filtro IIR depende de la frecuencia, así que un retardo constante solo alinea bien un tono y la alineación no elimina todo el error de fase.

La diferencia entre ambos valores es el costo de la latencia del filtro, no de la calidad del filtrado.

### Comparación justa: calidad frente a latencia

El retardo 0 del método espectral en la tabla **no es una propiedad del algoritmo**. Se obtiene porque el procesamiento se ejecutó fuera de línea sobre la señal completa y la salida se recortó para alinearla con la entrada (`y[hop:hop+len(x)]` en `src/filtrado_espectral.py`). En tiempo real, cada muestra de salida requiere haber recibido el bloque completo, así que el método espectral tiene una latencia del orden de un bloque: N = 1024 muestras, es decir, 64 ms a 16 kHz. El filtro temporal tiene 6 muestras de retardo (0.375 ms).

Por lo tanto:

| Criterio | Mejor método | Evidencia |
|---|---|---|
| Reducción de ruido (ΔSNR alineada) | Espectral en los tres casos | +18.87, +8.77 y +7.00 dB frente a +9.02, +5.34 y -0.08 dB |
| Distorsión de la señal útil | Espectral | 95.6-99.8 % de energía útil frente a 84.2-86.2 % |
| Latencia en tiempo real | Temporal | 0.375 ms frente a ≈64 ms |

La decisión del mejor compromiso depende de la latencia que tolere el sistema. Si se admiten decenas de milisegundos, el método espectral ofrece mejor reducción de ruido y menos distorsión. Si se exige una latencia de pocos milisegundos, solo el filtro temporal es viable, y en ese caso el caso 3 (ruido superpuesto a la banda útil) queda sin mejora (-0.08 dB).

### Atenuación del ruido dominante

La tarea pide la atenuación de la interferencia dominante. En el caso 1 es el tono de 1750 Hz. En los casos 2 y 3 no hay un tono, así que se define como la reducción de energía del ruido en la zona donde está:

$$A_{ruido} = 10\log_{10}\left(\frac{E_{ruido,\,entrada}}{E_{ruido,\,salida}}\right)$$

| Caso | Zona donde se mide la energía | Motivo |
|---|---|---|
| 1 | 1750 ± 20 Hz | Interferencia sinusoidal |
| 2 | Todo el espectro | El ruido blanco no tiene una banda dominante |
| 3 | 700 - 1500 Hz | Banda donde está el ruido |

**Cómo se obtiene el ruido a la salida.** Restar la referencia (`y - limpia`) mezcla el ruido residual con el retardo y la distorsión de los tonos. Con el filtro temporal el resultado sería engañoso: su MSE en el caso 2 (0.0594) es mayor que el de la señal sin filtrar (0.0144) solo por el retardo de 6 muestras. En su lugar, la señal limpia y el ruido se pasan por separado por el mismo sistema, y se mide la energía del ruido a la salida. Es válido porque ambos sistemas son lineales con parámetros fijos:

- **Temporal:** es un filtro lineal e invariante en el tiempo, por lo que `filtro(limpia) + filtro(ruido) = filtro(contaminada)`.
- **Espectral:** las ganancias por bloque se calculan una vez a partir de la señal contaminada y se aplican por separado a la limpia y al ruido. Con ganancias fijas el sistema es lineal.

La linealidad se verifica en el notebook con una tolerancia de 1e-10, y existe una prueba automática en `tests/test_filtros_espectrales.py`.

| Caso | Zona | Temporal (dB) | Espectral (dB) |
|---|---|---|---|
| 1 | 1750 ± 20 Hz | 30.36 | 19.61 |
| 2 | Espectro completo | 7.70 | 8.92 |
| 3 | 700 - 1500 Hz | 0.44 | 8.08 |

Fuente: `resultados/metricas.csv`, columna `atenuacion_ruido_dB`.

**Interpretación:**

- **Caso 1:** el filtro temporal atenúa más (30.4 dB) porque la interferencia cae en su banda de transición, cerca de un cero de transmisión. El espectral queda limitado a unos 20 dB por el piso de ganancia `g_min = 0.1` (-20 dB). Aun así el espectral logra mayor ΔSNR (+18.9 frente a +9.0 dB) porque conserva casi toda la energía útil (99.8 % frente a 84.2 %).
- **Caso 2:** el filtro pasa bajas de 1500 Hz solo puede eliminar el ruido por encima de su banda de paso. Para ruido blanco hasta 8 kHz, eso da 10·log10(8000/1500) ≈ 7.3 dB, cercano a los 7.7 dB medidos.
- **Caso 3:** el filtro temporal casi no atenúa el ruido (0.44 dB) porque la banda de 700 a 1500 Hz está dentro de su banda de paso. El método espectral logra 8.1 dB, pero a costa de distorsionar los tonos de 880 y 1320 Hz (95.6 % de energía útil), que es el compromiso descrito en la sección 2.2.

### Justificación de la zona de medición (caso 3)

La zona 700-1500 Hz no se eligió a ojo: es la banda con la que se generó el ruido (`banda_ruido_Hz` en `metadatos/configuraciones.json`), que es el parámetro del filtro pasa banda Butterworth que lo produce. Se comprobó de dos formas (celda de sensibilidad del notebook, sección 10.1):

| Zona | % de la energía del ruido | Temporal (dB) | Espectral (dB) |
|---|---|---|---|
| 700-1500 Hz (elegida) | 98.20 | 0.44 | 8.08 |
| 600-1600 Hz | 99.94 | 0.47 | 8.09 |
| 500-1800 Hz | 99.99 | 0.47 | 8.09 |
| Espectro completo | 100.00 | 0.47 | 8.09 |
| 880 ± 20 Hz (tono) | 5.58 | 0.42 | 2.26 |
| 1320 ± 20 Hz (tono) | 5.63 | 1.00 | 4.35 |

1. **La zona contiene el ruido:** el 98.2 % de su energía cae dentro de 700-1500 Hz. El resto corresponde a las colas del filtro que genera el ruido, que se extienden unos 100 Hz a cada lado de la banda.
2. **El resultado no depende de la zona:** al ampliarla hasta el espectro completo, la atenuación cambia menos de 0.03 dB. Por lo tanto, la elección de los límites no afecta la conclusión.

**Atenuación junto a los tonos.** En las vecindades de 880 y 1320 Hz (±20 Hz), el método espectral atenúa solo 2.3 y 4.4 dB, frente a 8.1 dB en la banda completa. Es la evidencia medida de lo descrito en la sección 2.2: el ruido que comparte frecuencia con un tono no se puede eliminar sin afectar el tono, porque la ganancia por bin actúa sobre la señal y el ruido a la vez. El valor de 8.1 dB de la banda completa es una media dominada por los bins donde no hay tono.

### Alcance y limitaciones de la definición

- **Qué mide:** cuánta energía de ruido sobrevive en la zona definida. No distingue entre el ruido que queda junto a un tono y el del resto de la banda; para eso se reporta aparte la atenuación junto a los tonos.
- **Zona conocida de antemano:** las zonas salen de la configuración de cada caso, que se conoce porque las señales las genera el grupo. Esto es válido para evaluar los métodos con ruido de referencia, pero no es lo que ocurriría en un sistema sin referencia.
- **Sistema real:** sin la señal limpia ni el ruido aislado no se puede medir esta atenuación. La zona del ruido tendría que estimarse sobre la señal contaminada con las características de la sección 2.2 (planicidad espectral, ancho de banda 95 %, pico sobre mediana). El método espectral ya hace una parte de esta estimación para las interferencias tonales (`estimar_ruido_e_interferencia`).
- **Medición en el microcontrolador:** al validar en el ESP32 con la señal de defensa, la zona se tomaría de los metadatos del caso, igual que aquí.

### Tiempo de ejecución por bloque

Se midió con `time.perf_counter()`, tomando la mediana de 15 repeticiones para descartar picos por otras tareas del sistema operativo (notebook, sección 10.2; columna `tiempo_bloque_us` de `resultados/metricas.csv`). El temporal se midió con `sosfilt` por bloques de 1024 muestras conservando el estado; no se usó la implementación manual porque es un lazo en Python y mediría al intérprete, no al algoritmo. El espectral se midió con `reducir_ruido_fft` completo dividido entre el número de bloques; con 50 % de solape cada bloque produce 512 muestras nuevas.

| Método | Muestras por bloque | Tiempo por bloque (µs) | Tiempo por muestra (µs) | Presupuesto de tiempo real (µs) |
|---|---|---|---|---|
| Temporal | 1024 | 41.7 | 0.041 | 64 000 |
| Espectral (N = 1024) | 512 | 48.2 | 0.094 | 32 000 |

Por muestra procesada, el método espectral cuesta unas 2.3 veces más que el temporal, porque calcula una FFT y una IFFT de 1024 puntos para producir 512 muestras nuevas. La calibración del espectral (una sola vez, con los primeros 2 s) tomó unos 5 ms. Son tiempos de PC con NumPy y SciPy compilados: sirven para comparar los métodos entre sí, no como tiempos del ESP32, y varían ligeramente entre ejecuciones.

### Tamaño de bloque del método espectral

Se repitió el método espectral con N = 256, 512, 1024, 2048 y 4096 en los tres casos (notebook, sección 10.4; datos en `resultados/testing/barrido_bloque.csv` y `resultados/testing/barrido_bloque_costos.csv`; figura `resultados/figuras/2_5_barrido_bloque.png`). El parámetro de agresividad se fijó en α = 2, igual que en la sección 2.4, para que N fuera la única variable. En cada N se recalculó la calibración (piso de ruido e interferencia con los primeros 2 s).

| N | Δf (Hz) | Latencia mínima (ms) | ΔSNR casos 1 / 2 / 3 (dB) | Energía útil casos 1 / 2 / 3 (%) | Tiempo por bloque en PC (µs) |
|---|---|---|---|---|---|
| 256 | 62.5 | 16 | 7.13 / 4.78 / 2.70 | 59.3 / 58.8 / 52.5 | 29 |
| 512 | 31.3 | 32 | 17.47 / 8.21 / 5.44 | 97.0 / 95.8 / 88.8 | 36 |
| 1024 | 15.6 | 64 | 18.87 / 8.77 / 7.00 | 99.8 / 99.2 / 95.6 | 49 |
| 2048 | 7.8 | 128 | 18.74 / 8.93 / 7.91 | 100.0 / 99.9 / 97.4 | 72 |
| 4096 | 3.9 | 256 | 18.54 / 8.93 / 8.11 | 100.0 / 100.1 / 98.2 | 126 |

La resolución es Δf = fs/N y, en tiempo real, la latencia es al menos un bloque, N/fs. Los tiempos son de PC con NumPy y SciPy; sirven para comparar tamaños entre sí, no como tiempos del ESP32, y varían un poco entre ejecuciones.

**Interpretación:**

- **N = 256 no es utilizable:** con Δf = 62.5 Hz los bins no separan bien los tonos del ruido, y se pierde entre el 41 y el 47 % de la energía útil.
- **N = 512 a 1024:** subir de 512 a 1024 mejora el ΔSNR entre 0.6 y 1.6 dB y la energía útil entre 2.8 y 6.8 puntos, a cambio de duplicar la latencia de 32 a 64 ms.
- **N ≥ 2048:** en los casos 1 y 2 la mejora es de ±0.2 dB, mientras la latencia vuelve a duplicarse (128 y 256 ms). El caso 3 sí sigue mejorando (de 7.0 a 8.1 dB), porque separar un tono del ruido que lo rodea requiere bins más finos.
- **Propuesta:** N = 1024 como tamaño de bloque base, con N = 512 como alternativa si el sistema exige una latencia menor a 64 ms, aceptando una pérdida de entre 0.6 y 1.6 dB de ΔSNR.

### Comprobación de la dependencia con α

Como el mejor α podría cambiar con N, el barrido se repitió con α = 1 y α = 4. La forma de las curvas es la misma: N = 256 es el peor en los tres valores de α, y la mejora de 512 a 1024 se mantiene. La conclusión sobre N no depende de α.

### Elección de α

El parámetro α es el factor de sustracción de la ganancia por bin, $G^2 = \max\left(1 - \alpha P_n/P,\; g_{min}^2\right)$. Se barrió α con N = 1024 (notebook, sección 10.5; `resultados/testing/barrido_alpha.csv`):

| α | ΔSNR caso 1 (dB) | ΔSNR caso 2 (dB) | ΔSNR caso 3 (dB) | Energía útil caso 3 (%) |
|---|---|---|---|---|
| 1 | 18.83 | 4.61 | 4.14 | 98.8 |
| 2 | 18.87 | 8.77 | 7.00 | 95.6 |
| 3 | 18.91 | 11.91 | 8.09 | 92.8 |
| 4 | 18.95 | 13.83 | 8.03 | 90.3 |
| 6 | 19.05 | 15.10 | 6.91 | 86.0 |
| 8 | 19.16 | 15.16 | 5.78 | 82.4 |
| 12 | 19.24 | 14.77 | 4.26 | 77.3 |
| 32 | 18.17 | 10.69 | 1.38 | 70.3 |

| Caso | α que maximiza el ΔSNR | ΔSNR máximo (dB) | ΔSNR con α = 2 (dB) | Energía útil en el máximo (%) |
|---|---|---|---|---|
| 1 | 12 | 19.24 | 18.87 | 99.1 |
| 2 | 8 | 15.16 | 8.77 | 95.7 |
| 3 | 3 | 8.09 | 7.00 | 92.8 |

**Interpretación:**

- **α = 2 no maximiza el ΔSNR en ningún caso.** Es un valor conservador, y el máximo cae entre α = 3 y α = 12 según el caso.
- **Existe un máximo porque hay dos efectos opuestos.** Con α pequeño queda ruido sin restar. Con α grande se resta de más y se empieza a eliminar señal útil, lo que reduce la energía útil y vuelve a bajar el ΔSNR.
- **El caso 3 es el más sensible.** Su ruido comparte frecuencias con los tonos, por lo que su máximo es el más temprano (α = 3) y cae rápido después: con α = 8, que es el mejor valor para el caso 2, el caso 3 baja de 8.1 a 5.8 dB.
- **El caso 1 casi no depende de α** (entre 18.2 y 19.2 dB en todo el rango), porque la interferencia se elimina con la máscara de ancho fijo y no con la resta.
- **No existe un α óptimo único.** Un valor que maximiza un caso perjudica a otro.
- **El sistema real no puede maximizar el ΔSNR.** Calcularlo requiere la señal limpia, que el microcontrolador no tiene. Por eso α debe fijarse al diseñar, con un valor que funcione bien en general.

Se mantiene α = 2 en toda la comparación para ser coherentes con la sección 2.4, a costa de no obtener el máximo ΔSNR posible. Un α entre 3 y 4 mejoraría el resultado en los casos 2 y 3, a cambio de más distorsión de la señal útil.

### Memoria aproximada en el ESP32

Se estimó con un **modelo analítico**: se cuentan los buffers que necesita cada algoritmo (notebook, sección 10.6; datos en `resultados/testing/memoria_espectral_por_N.csv`, `resultados/testing/memoria_calibracion.csv` y la columna `memoria_algoritmo_bytes` de `resultados/metricas.csv`). No es una medición en el microcontrolador; debe confirmarse con el mapa de memoria del compilador (`idf.py size`) al portar el código.

**Datos externos usados:**

- **Memoria del ESP32:** 520 KB de SRAM interna, de los cuales 320 KB son DRAM y 200 KB IRAM. El máximo de DRAM asignada estáticamente es 160 KB; los otros 160 KB solo están disponibles como heap. Si se usa Bluetooth, se pierden 64 KB de DRAM ([ESP-IDF, Memory Types](https://docs.espressif.com/projects/esp-idf/en/stable/esp32/api-guides/memory-types.html), versión 6.1; conviene verificar la cifra en la versión de ESP-IDF que se use).
- **FFT de ESP-DSP:** la FFT compleja se calcula *in situ* sobre un arreglo intercalado de N puntos complejos (2N floats). El tamaño de la tabla de factores se pasa a `dsps_fft2r_init_fc32` en palabras float, y su máximo por defecto es 4096 ([dsps_fft2r.h](https://github.com/espressif/esp-dsp/blob/master/modules/fft/include/dsps_fft2r.h)). Para señales reales, el ejemplo oficial usa una FFT compleja de N/2 puntos más una conversión a espectro real ([fft4real](https://github.com/espressif/esp-dsp/blob/master/examples/fft4real/main/dsps_fft4real_main.c)).

**Supuestos:** formato float32; muestras de E/S en int16 con doble buffer; las ganancias se calculan bin por bin sobre el mismo buffer, sin guardar la matriz de ganancias que sí guarda la versión de Python. El tamaño de la tabla de factores (N/2 floats para N/2 puntos complejos) es una **inferencia** compatible con la documentación y con un reporte externo de 16 KB para 4096 puntos, pero la documentación oficial no da la fórmula explícita.

**Desglose del método espectral con N = 1024:**

| Componente | Bytes |
|---|---|
| Ventana raíz de Hann | 4096 |
| Buffer de trabajo FFT (512 complejos, in situ) | 4096 |
| Tabla de factores (twiddle) | 4096 |
| Piso de ruido Pn | 2052 |
| Solape de salida (overlap-add) | 2048 |
| Muestras previas de entrada | 2048 |
| Lista de interferencias (8 frecuencias) | 32 |
| **Algoritmo** | **18 468** |
| E/S con doble buffer, entrada y salida (int16) | 4096 |
| **Total en régimen** | **22 564** |

**Comparación de métodos** (sin buffers de E/S): el filtro temporal necesita 15 coeficientes y 6 estados, es decir **84 bytes**, frente a **18 468 bytes** del método espectral. El espectral usa unas 220 veces más memoria.

**Efecto del tamaño de bloque** (en régimen, con E/S):

| N | Total (KB) | % de la DRAM (320 KB) |
|---|---|---|
| 256 | 8.5 | 2.7 |
| 512 | 13.0 | 4.1 |
| 1024 | 22.0 | 6.9 |
| 2048 | 40.0 | 12.5 |
| 4096 | 80.0 | 25.0 |

**Calibración.** El método estima el piso de ruido y las interferencias una sola vez, con los primeros 2 s de señal. En Python esa estimación usa una FFT de 8192 puntos. Se comprobó que con una FFT de calibración de 1024, 2048, 4096 u 8192 puntos el resultado es idéntico en los tres casos (misma interferencia detectada y mismo ΔSNR), por lo que se puede reducir. La memoria extra, sin reutilizar buffers, es:

| FFT de calibración | Pico (KB) | % de DRAM (320 KB) | % del máximo estático (160 KB) |
|---|---|---|---|
| 1024 | 32.0 | 10.0 | 20.0 |
| 2048 | 46.0 | 14.4 | 28.8 |
| 4096 | 74.0 | 23.1 | 46.3 |
| 8192 | 130.0 | 40.6 | 81.3 |

**Interpretación:**

- **La memoria no es el factor limitante.** Con N = 1024 el método espectral usa 22 KB, el 6.9 % de la DRAM. Incluso N = 4096 cabe (25 %); lo que limita el tamaño de bloque es la latencia, no la memoria.
- **La calibración con la FFT de 8192 puntos es la parte más cara:** su pico de 130 KB ocuparía el 81 % del máximo de DRAM estática. Se propone calibrar con 2048 puntos (pico de 46 KB), con el mismo resultado en los tres casos. Además, los buffers de calibración se pueden liberar o reutilizar después, así que el costo permanente es el de régimen (22 KB).
- **Limitaciones:** solo se probó con tres casos, con tonos útiles armónicos y una interferencia a 1750 Hz; con otros ruidos habría que repetir la comprobación del tamaño de calibración. El modelo no incluye la pila de FreeRTOS, el código, ni la pila de comunicación entre nodos (WiFi o Bluetooth), que puede consumir una parte importante de la DRAM.
