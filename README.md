# Replicación Campbell (1994) — Secciones 2, 3a, 3b y 4

Este repositorio contiene varios archivos Dynare (.mod) que implementan versiones del modelo RBC de Campbell (1994):

- `campbel_sec2_fixed_labor.mod`  — Sección 2: oferta laboral fija (L = 1), variables en log, escrito en forma no lineal usando `exp(...)`.
- `campbel_sec3a_separable_labor.mod` — Sección 3a: oferta laboral variable, utilidad separable `u(c,n)=log c + phi_n log(1-n)`.
- `campbel_sec3b_nonseparable_labor.mod` — Sección 3b: oferta laboral variable, utilidad no separable (Cobb–Douglas consumo/ocio compuesto).
- `campbel_sec4_gov_spending.mod` — Sección 4: introduce gasto público `g` y financiación por impuesto proporcional `tau` y shock de gasto.


Resumen rápido de ecuaciones comunes
- Producción (Cobb–Douglas, variables en log y escrito en niveles con `exp`):
  $ Y_t = A_t K_{t-1}^\alpha N_t^{1-\alpha}$  
  (en los archivos: `exp(y) = exp(a) * exp(k(-1))^alpha * exp(n)^(1-alpha)`)
- Restricción de recursos:    
  $Y_t = C_t + I_t + G_t$   
  (`exp(y) = exp(c) + exp(i) + exp(g)`).
- Acumulación de capital:  
  $ K_t = (1-\delta)K_{t-1} + I_t$     
  (`exp(k) = (1-delta) * exp(k(-1)) + exp(i)`).
- Euler (CRRA/log caso sigma=1):  
  $$C_t^{-\sigma} = \beta E_t\left[C_{t+1}^{-\sigma}\left((1-\tau)\frac{\alpha Y_{t+1}}{K_t} + 1-\delta\right)\right]$$
- FOC trabajo‑ocio (separable): $$\dfrac{\phi_n}{1-N_t} = \dfrac{1}{C_t}(1-\tau)\dfrac{(1-\alpha)Y_t}{N_t}$$

Notas importantes
- Convención de logs: en estos .mod muchas variables (`y,c,k,i,n,a,g`) están definidas en logaritmos. Para trabajar en niveles las ecuaciones usan `exp(var)`; esto es una elección del modelador, no algo que Dynare determine automáticamente.
- `g_y` en algunos archivos es la fracción de `Y` usada en `initval` para fijar el estado estacionario; en la implementación original `g` puede ser exógeno (AR(1) en logs). Si quieres que `G` sea estrictamente igual a `\tau Y` usar la versión explícita o añadir la ecuación correspondiente.

Salida y ficheros generados
- Dynare genera automáticamente archivos `.log`, `.mat` y estructuras en memoria como `oo_` y `M_`.
- Los scripts de Campbell han sido modificados para guardar figuras en `plots/` (PNG). Revisa `plots/` tras ejecutar los .mod.
- Los IRFs están en `oo_.irfs` (por ejemplo `oo_.irfs.y_e_a` para la respuesta de `y` al shock `e_a`).

Cómo ejecutar los scripts con Dynare

Requisitos previos
- Tener Dynare instalado y accesible desde MATLAB o Octave.
- Abrir MATLAB o Octave y situar el directorio de trabajo en esta carpeta: la que contiene los `.mod`.

Ejecución desde MATLAB / Octave (prompt)
1. Abrir MATLAB o Octave.
2. Cambiar al directorio `Final` (o la carpeta donde están los .mod):

```matlab
cd '/ruta/a/Final'
```

3. Ejecutar Dynare con el nombre del archivo (sin extensión):

```matlab
dynare campbel_sec2_fixed_labor
dynare campbel_sec3a_separable_labor
dynare campbel_sec3b_nonseparable_labor
dynare campbel_sec4_gov_spending
```

Nota: puedes ejecutar cualquiera de los `.mod` según qué sección quieras correr. Dynare cargará, calculará el estado estacionario (o usará `initval`), hará `check` y correrá `stoch_simul` tal como está programado en el archivo.

Ejecución desde la línea de comandos (Octave + Dynare)

Si prefieres lanzar desde terminal usando Octave (suponiendo Dynare en el path de Octave):

```bash
octave --eval "dynare campbel_sec2_fixed_labor"
```

Recomendaciones y resolución de errores comunes
- Si `steady` no converge: comprobar `initval`, elegir valores iniciales razonables (p. ej. `n ~ log(0.3)`), o usar `steady(solve_algo=4)` u opciones de `steady` que prueben métodos numéricos distintos.
- Si `check` muestra problemas de determinacy o eigenvalues fuera de rango: revisar parámetros (`beta`, `alpha`, `delta`) y restricciones fiscales (si `tau>0` asegurar coherencia de `G` y recaudación).
- Para experimentar con `tau` (impuesto proporcional) y comparar efectos: modifica `tau` en el bloque `parameters` y lanza el `.mod`. Si quieres que `G` se ajuste automáticamente a la recaudación, usa `campbel_sec4_gov_spending_explicit.mod` que implementa `exp(g) = tau*exp(y) + exp(T)`.
