# Método de Momentos: Parche Conductor en Interfaz Dieléctrica

## Descripción

Implementación del **Método de Momentos (MoM)** para resolver el problema electrostático de un parche conductor rectangular ubicado en la interfaz entre dos medios dieléctricos.

**Configuración del problema:**
- Parche conductor rectangular en $z = 0$
- Medio 1 ($z > 0$): Aire con $\varepsilon_1 = \varepsilon_0$
- Medio 2 ($z < 0$): Dieléctrico con $\varepsilon_2 = 4\varepsilon_0$
- Voltaje aplicado: $V_0 = 1$ V
- Objetivo: Determinar $\sigma(x,y)$ y campos eléctricos

---

## Formulación Matemática

### Ecuación Integral

$$\iint_S \sigma(x',y') \, G(x,y,z|x',y',0) \, dS' = V_0$$

donde:
- $\sigma(x',y')$: densidad de carga superficial (incógnita)
- $G$: función de Green para interfaz dieléctrica
- $S$: superficie del conductor

### Función de Green (Método de Carga Imagen)

**Región 1** ($z > 0$):
$$G = \frac{1}{4\pi\varepsilon_1} \left[\frac{1}{R_{\text{directo}}} + \frac{q_{\text{imagen}}}{R_{\text{imagen}}}\right]$$

**Región 2** ($z < 0$):
$$G = \frac{1}{4\pi\varepsilon_2} \left[\frac{q_{\text{trans}}}{R_{\text{directo}}}\right]$$

donde:
$$q_{\text{imagen}} = \frac{\varepsilon_1 - \varepsilon_2}{\varepsilon_1 + \varepsilon_2}, \quad q_{\text{trans}} = \frac{2\varepsilon_1}{\varepsilon_1 + \varepsilon_2}$$

### Discretización MoM

**Funciones base:** Constantes por tramos (funciones pulso)

**Sistema matricial:**
$$[Z]\{\sigma\} = \{V_0\}$$

**Término propio (diagonal):** Para celda cuadrada de lado $L$:
$$Z_{ii} = \frac{1}{\varepsilon_1 + \varepsilon_2} \cdot \frac{L}{\pi} \cdot \left[\ln(1+\sqrt{2}) + \sqrt{2}\right]$$

**Términos mutuos:** Cuadratura gaussiana 3×3

---

## Archivos del Proyecto

```
/
├── solver.c              # Solver MoM en C
├── visualizador.py           # Visualización Python
├── charge_distribution.csv   # Resultados: distribución de carga
├── field_data.csv           # Resultados: campos eléctricos
└── README.md                # Este archivo
```

--

## Compilación y Ejecución

### 1. Compilar el Solver

```bash
gcc -o mom_solver mom_solver.c -lm -O3
```

### 2. Ejecutar Simulación

```bash
./mom_solver
```

**Salida esperada:**
```
Total charge: 9.776e-12 C
Capacitance (MoM): 9.776e-12 F
Max charge density: 3.391e-09 C/m²
Min charge density: 5.124e-10 C/m²
Charge density ratio (max/min): 6.62
```

### 3. Generar Visualizaciones

```bash
python visualizador.py
```

**Requisitos Python:**
```bash
pip install numpy pandas matplotlib seaborn
```

---

## Resultados

### Parámetros de la Simulación

| Parámetro | Valor |
|-----------|-------|
| Tamaño del parche | 10 × 10 cm |
| Voltaje aplicado | 1.0 V |
| Discretización | 20 × 20 celdas |
| $\varepsilon_1$ (aire) | $8.854 \times 10^{-12}$ F/m |
| $\varepsilon_2$ (dieléctrico) | $3.542 \times 10^{-11}$ F/m |
| Contraste dieléctrico | $\varepsilon_2/\varepsilon_1 = 4.0$ |

### Resultados Numéricos

- **Carga total:** $Q = 9.78$ pC
- **Capacitancia:** $C = Q/V_0 = 9.78$ pF
- **Densidad de carga máxima:** $\sigma_{\max} = 3.39$ nC/m² (bordes)
- **Densidad de carga mínima:** $\sigma_{\min} = 0.51$ nC/m² (centro)
- **Razón máx/mín:** 6.62 (evidencia de singularidad de borde)
- **Residuo numérico:** $6.52 \times 10^{-16}$ (excelente convergencia)

### Visualizaciones Generadas

El script de Python genera 7 gráficas:

1. **`distribucion_carga.png`**
   - Mapa de calor 2D
   - Superficie 3D
   - Cortes por líneas centrales

2. **`estadisticas_carga.png`**
   - Histograma de distribución
   - Diagrama de caja con estadísticas

3. **`campo_potencial.png`**
   - Potencial en tres planos ($z$ = 1, 3, 5 cm)

4. **`potencial_xz.png`**
   - Corte transversal del potencial
   - Perfil vertical

5. **`campo_electrico_magnitud.png`**
   - Magnitud de $|\vec{E}|$ en múltiples planos

6. **`campo_electrico_vectores.png`**
   - Campo vectorial 2D
   - Campo vectorial en corte $xz$

7. **`campo_electrico_3d.png`**
   - Visualización 3D de magnitud
   - Líneas de campo

---

## Características Físicas Observadas

### Singularidad de Borde
La densidad de carga aumenta hacia los bordes y esquinas del parche (factor ~6.6×), comportamiento esperado en conductores.

### Continuidad del Potencial
El potencial es continuo en la interfaz $z=0$: $V_1(z=0^+) = V_2(z=0^-)$

### Discontinuidad del Campo Normal
La componente normal del campo eléctrico presenta discontinuidad en la interfaz:
$$\varepsilon_1 E_{1z}(z=0^+) = \varepsilon_2 E_{2z}(z=0^-)$$

### Simetría
Para un parche cuadrado, la distribución de carga es simétrica respecto a $x=0$ y $y=0$.

---

## Parámetros Configurables

En `mom_solver.c` (líneas 443-465):

```c
/* Permitividades */
p.eps1 = 8.854e-12;        // Aire
p.eps2 = 4.0 * 8.854e-12;  // Dieléctrico

/* Geometría */
p.patch_width = 0.1;       // 10 cm
p.patch_height = 0.1;      // 10 cm

/* Discretización */
p.N_x = 20;                // Celdas en x
p.N_y = 20;                // Celdas en y

/* Dominio de evaluación */
p.x_min = -0.15;  p.x_max = 0.15;  // ±15 cm
p.y_min = -0.15;  p.y_max = 0.15;
p.z_min = -0.1;   p.z_max = 0.1;   // -10 a +10 cm
```

### Estudio de Convergencia

| Malla | Incógnitas | Tiempo | Memoria |
|-------|-----------|--------|---------|
| 10×10 | 100 | ~2 s | ~10 KB |
| 20×20 | 400 | ~8 s | ~13 KB |
| 40×40 | 1600 | ~3 min | ~200 KB |
| 60×60 | 3600 | ~15 min | ~1 MB |

---

## Verificación Numérica

### Checks de Auto-Consistencia

1. **Conservación de carga:** $\iint \sigma \, dS = Q$ constante
2. **Condición de frontera:** $V = V_0$ en el conductor 
3. **Simetría:** Solución simétrica para parche cuadrado
4. **Residuo:** $|[Z]\{\sigma\} - \{V_0\}| < 10^{-15}$

### Validación de Orden de Magnitud

La capacitancia obtenida ($C \approx 10$ pF) es consistente con la estimación de orden:

$$C \sim (\varepsilon_1 + \varepsilon_2) \cdot L \sim 5\varepsilon_0 \cdot 0.1 \text{ m} \sim 4.4 \text{ pF}$$

El factor ~2 de diferencia se debe a efectos geométricos (cuadrado vs. circular).

---

## Implementación Técnica

### Aspectos Clave del Código

**1. Término Propio Analítico**
- Evita singularidad en integración numérica
- Fórmula exacta para celda cuadrada
- Crucial para precisión y estabilidad

**2. Cuadratura Gaussiana**
- 3×3 puntos para términos mutuos
- Balance entre precisión y eficiencia

**3. Eliminación Gaussiana con Pivoteo**
- Solver directo estable
- Adecuado para $N < 2000$

**4. Gradiente Numérico**
- Campo eléctrico: $\vec{E} = -\nabla V$
- Diferencias finitas centradas

---

## Referencias

1. **Harrington, R.F. (1968)**: "Field Computation by Moment Methods", IEEE Press.

2. **Jackson, J.D. (1999)**: "Classical Electrodynamics", 3rd Ed., Wiley, Cap. 4.

3. **Mosig, J.R. (2024)**: "Roger F. Harrington and the Method of Moments", *IEEE Antennas and Propagation Magazine*, Vol. 66, No. 2.

---

## Autores y Licencia

**Proyecto académico** - Curso de Electrodinamica I

Código disponible bajo licencia MIT para fines educativos.

---

**Última actualización:** Diciembre 2024
