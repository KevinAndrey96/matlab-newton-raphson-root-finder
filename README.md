# Método de Newton-Raphson para Búsqueda de Raíces en MATLAB y GNU Octave

[![MATLAB](https://img.shields.io/badge/Language-MATLAB%20%2F%20GNU%20Octave-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Implementación modular del algoritmo iterativo de **Newton-Raphson** para la localización y aproximación de raíces en funciones continuas y no lineales ($f(x) = 0$). Incorpora diferenciación numérica mediante diferencias centradas y control estricto de convergencia con tolerancia residual de alta precisión ($|f(x)| \le 10^{-10}$).

---

## 📐 Fundamento Teórico y Matemático

El método de Newton-Raphson es uno de los algoritmos de resolución numérica de ecuaciones escalares no lineales más eficientes debido a su tasa de convergencia cuadrática ($p = 2$) en las proximidades de una raíz simple.

### 1. Formulación Iterativa

A partir del desarrollo en serie de Taylor de primer orden de $f(x)$ alrededor de una estimación inicial $x_k$:

$$f(x) \approx f(x_k) + f'(x_k)(x - x_k)$$

Imponiendo la condición de raíz $f(x) = 0$ y despejando $x = x_{k+1}$, se obtiene la regla de actualización:

$$x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)}$$

### 2. Diferenciación Numérica Centrada

En escenarios prácticos donde la derivada analítica $f'(x)$ es costosa de derivar o evaluar simbólicamente, se aproxima mediante la fórmula de diferencias finitas centradas con paso $h = 10^{-5}$:

$$f'(x) \approx \frac{f(x + h) - f(x - h)}{2h} + O(h^2)$$

### 3. Criterio de Parada y Tolerancia

El algoritmo itera continuamente hasta satisfacer el criterio de tolerancia de residuo funcional:

$$|f(x_k)| \le 10^{-10}$$

Garantizando que el punto $x^*$ encontrado aproxima el cero de la función con precisión de al menos diez cifras decimales.

---

## 📂 Estructura del Proyecto

```text
newton-raphson/
├── f.m           # Funcion continua objetivo: f(x) = sin(x)
├── df.m          # Derivada numerica por diferencias finitas centradas con h = 10^-5
└── newt.m        # Bucle iterativo de Newton-Raphson y reporte de raiz encontrada
```

---

## 🚀 Requisitos y Ejecución

### Prerrequisitos

- **MATLAB** (R2016a o superior) o **GNU Octave** (4.0 o superior).

### Ejecución en GNU Octave (Terminal)

```bash
# Iniciar Octave en la carpeta del repositorio
octave -q

# Invocar la función con un punto de inicio x0 cercano a pi (ej: 3.0)
octave:1> newt(3.0)
```

### Ejecución en MATLAB

1. Abre MATLAB y coloca la ruta de trabajo en la carpeta `newton-raphson`.
2. En la consola de comandos (`Command Window`), ejecuta la función proporcionando una estimación inicial:

```matlab
% Aproximar la raíz de sin(x) comenzando en x0 = 3.0
raiz = newt(3.0);
```

### Salida de Ejemplo

```text
La raiz es:       3.141592653590 con un error de:       0.000000000000 
```

---

## 💡 Personalización de la Función

Para resolver una función no lineal diferente, simplemente edita el archivo `f.m`:

```matlab
function y = f(x)
  % Ejemplo: f(x) = x^2.01 - 13*x + 30
  y = x^2.01 - 13*x + 30;
endfunction
```

Al utilizar diferencias numéricas en `df.m`, no es necesario recalcular la derivada manualmente.

---

## 📄 Licencia

Distribuido bajo la Licencia MIT. Consulta el archivo `LICENSE` para más información.
