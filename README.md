# Método de Newton-Raphson para encontrar raíces en MATLAB y GNU Octave

[![CI](https://github.com/KevinAndrey96/matlab-newton-raphson-root-finder/actions/workflows/ci.yml/badge.svg)](https://github.com/KevinAndrey96/matlab-newton-raphson-root-finder/actions/workflows/ci.yml)
![MATLAB](https://img.shields.io/badge/MATLAB-R2026a-0076A8)
![GNU Octave](https://img.shields.io/badge/GNU%20Octave-11-0790C0)
[![Licencia: MIT](https://img.shields.io/badge/licencia-MIT-yellow.svg)](LICENSE)

Funciones para resolver ecuaciones no lineales $f(x) = 0$ en MATLAB y GNU Octave: **Newton-Raphson** con derivada analítica o numérica, límite de iteraciones, corrección por **multiplicidad** para raíces repetidas, un método **Newton con bisección de respaldo** que converge siempre que exista un intervalo con cambio de signo, y un estimador del **orden de convergencia** observado. Sin toolboxes, probado en Octave 11 y en MATLAB en cada cambio.

## Contenido

- [Fundamento matemático](#fundamento-matemático)
- [Uso](#uso)
- [Resultados](#resultados)
- [Estructura](#estructura)
- [Pruebas y calidad](#pruebas-y-calidad)
- [Historia del proyecto](#historia-del-proyecto)

## Fundamento matemático

Linealizando $f$ alrededor de la estimación $x_k$ con el polinomio de Taylor de primer orden e igualando a cero se obtiene la iteración

$$
x_{k+1} = x_k - m\,\frac{f(x_k)}{f'(x_k)},
$$

donde $m = 1$ en el método clásico. Cerca de una raíz simple el error cumple $e_{k+1} \approx \frac{f''(r)}{2 f'(r)}\,e_k^2$: el número de cifras correctas se duplica en cada paso (convergencia cuadrática). Si además $f''(r) = 0$, como ocurre con $\sin x$ en $\pi$, la convergencia es cúbica.

**Raíces repetidas.** Si $r$ tiene multiplicidad $m > 1$, el método clásico solo converge linealmente con razón $(m - 1)/m$. Multiplicar el paso por $m$ recupera la convergencia cuadrática.

**Derivada numérica.** Cuando no se entrega $f'$, se usa la diferencia centrada

$$
f'(x) \approx \frac{f(x + h) - f(x - h)}{2h}, \qquad h = \varepsilon^{1/3}\max(1, |x|),
$$

con el paso que equilibra el error de truncamiento $O(h^2)$ y el de redondeo $O(\varepsilon/h)$.

**Newton con bisección de respaldo.** Dado $[a, b]$ con $f(a) f(b) < 0$, en cada paso se intenta Newton y se recurre a bisección si el punto sale del intervalo o si el residuo no se reduce lo suficiente. El intervalo siempre contiene una raíz, de modo que el método no diverge y conserva la velocidad de Newton cerca de la solución.

**Orden observado.** Con $e_k = |x_k - r|$, el orden se estima como $p_k = \log(e_{k+1}/e_k) / \log(e_k/e_{k-1})$.

## Uso

```matlab
addpath('src');

% Derivada numérica automática
[raiz, info] = newton_raphson(@sin, 3);

% Derivada analítica y tolerancia más estricta
opciones = struct('derivative', @(x) 3*x.^2 - 2, 'tolerance', 1e-15);
[raiz, info] = newton_raphson(@(x) x.^3 - 2*x - 5, 2, opciones);
convergence_order(info.history, raiz)

% Raíz triple con corrección de multiplicidad
opciones = struct('derivative', @(x) 3*(x - 2).^2, 'multiplicity', 3);
raiz = newton_raphson(@(x) (x - 2).^3, 5, opciones);

% Método robusto sobre un intervalo con cambio de signo
[raiz, info] = newton_bisection(@atan, @(x) 1 ./ (1 + x.^2), -10, 3);
```

| Función | Descripción |
| :--- | :--- |
| `newton_raphson(f, x0, opciones)` | Newton clásico. Opciones: `derivative`, `tolerance` (residuo, `1e-10`), `step_tolerance` (paso relativo, `1e-14`), `max_iterations` (`50`) y `multiplicity` (`1`). Devuelve `converged`, `iterations`, `residual` y `history`. |
| `newton_bisection(f, df, a, b, tol, max_iter)` | Newton protegido por bisección. Devuelve el número de pasos de cada tipo y el intervalo final. |
| `convergence_order(history, raiz)` | Orden de convergencia observado en cada paso. |

Los errores usan identificadores `newton:*` (`newton:zeroDerivative`, `newton:nonFinite`, `newton:invalidBracket`, entre otros) para poder capturarlos con `try/catch`.

## Resultados

`examples/newton_demo.m` reproduce los ejercicios originales y compara los casos de raíz simple y repetida:

```text
sin(x) from x0 = 3 (original exercise)
   k                    x(k)    |x(k) - r|    order
   0       3.000000000000000     1.416e-01      NaN
   1       3.142546543082244     9.539e-04      NaN
   2       3.141592653300413     2.894e-10    3.002
   3       3.141592653589793     0.000e+00      NaN

x^2.01 - 13x + 30 from x0 = 2
   k                    x(k)    |x(k) - r|    order
   0       2.000000000000000     1.014e+00      NaN
   1       2.896759104768069     1.177e-01      NaN
   2       3.012435835994965     1.997e-03    1.892
   3       3.014432027919935     5.950e-07    1.992
   4       3.014432622879442     0.000e+00      NaN

(x - 1)^2 without multiplicity: 25 iterations, order 1.000
(x - 1)^2 with multiplicity 2: 1 iteration

atan(x) on [-10, 3] with the safeguarded method: root 2.555e-19
  3 Newton steps, 2 bisection steps (plain Newton from 3 diverges)
```

```bash
octave-cli examples/newton_demo.m
matlab -batch "run('examples/newton_demo.m')"
```

## Estructura

```text
src/
  newton_raphson.m       Newton-Raphson con opciones y multiplicidad
  newton_bisection.m     Newton protegido con bisección
  convergence_order.m    orden de convergencia observado
examples/
  newton_demo.m          ejercicios originales y raíces repetidas
  print_newton_run.m     tabla de iteraciones del ejemplo
tests/
  run_tests.m            ejecutor compatible con MATLAB y Octave
  test_*.m               suites de pruebas
```

## Pruebas y calidad

```bash
octave-cli --eval "cd tests; run_tests"
matlab -batch "cd tests; run_tests"
```

17 pruebas verifican convergencia cuadrática en raíces simples, convergencia lineal en raíces dobles, la corrección por multiplicidad, el agotamiento del presupuesto de iteraciones, derivadas nulas, divergencia, el método protegido en funciones donde Newton diverge y la validación de entradas. El pipeline de GitHub Actions ejecuta MISS_HIT (`mh_style`, `mh_lint`, `mh_metric`), las pruebas en el contenedor oficial de GNU Octave 11 y las mismas pruebas en MATLAB con las acciones de MathWorks.

## Historia del proyecto

El repositorio empezó en 2019 como un ejercicio de clase con tres archivos: una función fija, su derivada numérica y un ciclo `while` sin límite de iteraciones. En 2026 se reescribió como una librería reutilizable con funciones como parámetro, validación de entradas, manejo de derivadas nulas y divergencia, el método protegido con bisección, pruebas automatizadas y CI en MATLAB y Octave. El código original se conserva en el tag [`v0.1-matlab-coursework-2019`](https://github.com/KevinAndrey96/matlab-newton-raphson-root-finder/tree/v0.1-matlab-coursework-2019).

## Autor

**Kevin Andrey Herrera Silva** · [GitHub](https://github.com/KevinAndrey96)

## Licencia

Distribuido bajo la [licencia MIT](LICENSE).
