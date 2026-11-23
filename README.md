[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=emanuelalejandro9-web/MMP2)

# Práctica 2 — Caos en sistemas biológicos y control de células patológicas

Repositorio académico de la asignatura **Modelado Matemático**, correspondiente a la **Maestría en Ciencias de la Ingeniería** del Tecnológico Nacional de México / Instituto Tecnológico de Tijuana. La práctica integra modelado mecanicista de sistemas dinámicos no lineales en tres dimensiones, comparación de métodos numéricos de diferencias finitas, cálculo simbólico y numérico de puntos de equilibrio, análisis de estabilidad asintótica local mediante la matriz Jacobiana, localización de conjuntos compactos invariantes y análisis de estabilidad asintótica global mediante el método directo de Lyapunov[cite: 5, 16, 17, 18].

El caso de estudio utiliza el modelo tridimensional de cáncer propuesto por **Itik y Banks (2010)** para representar la dinámica no lineal entre tres poblaciones celulares: células patológicas o tumorales, células normales o sanas y células efectoras del sistema inmune (citotóxicas o *Natural Killer*), evaluando tanto su comportamiento caótico natural como la dosis crítica de inmunoterapia $\rho_i$ necesaria para erradicar la población patológica[cite: 16, 17].

**Docente:** Dr. Paul Antonio Valle Trujillo  
**Departamento:** Ingeniería Eléctrica y Electrónica  
**Institución:** Tecnológico Nacional de México / Instituto Tecnológico de Tijuana  
**Asignatura:** Modelado Matemático  
**Programa:** Maestría en Ciencias de la Ingeniería  

**Alumno:** Flores Gonzalez Emmanuel Alejandro. M26210085. m26210085@tectijuana.edu.mx

---

<a id="contenido"></a>
## Contenido

1. [Resumen](#resumen)
2. [Objetivos](#objetivos)
3. [Modelo matemático](#modelo)
4. [Regímenes dinámicos y parámetros](#regimenes)
5. [Análisis matemático: Positividad y equilibrios](#analisis-equilibrios)
6. [Estabilidad asintótica local](#estabilidad-local)
7. [Conjuntos compactos invariantes y estabilidad global](#estabilidad-global)
8. [Métodos numéricos](#metodos)
9. [Flujo de trabajo computacional](#flujo)
10. [Actividades y secciones del cuaderno MATLAB](#actividades)
11. [Funciones auxiliares](#funciones)
12. [Estructura del repositorio](#estructura)
13. [Requisitos](#requisitos)
14. [Reproducibilidad](#reproducibilidad)
15. [Referencias](#referencias)

---

<a id="resumen"></a>
## Resumen

Esta práctica estudia la dinámica no lineal y el control de un sistema biológico tridimensional descrito por ecuaciones diferenciales ordinarias de primer orden que exhibe bifurcaciones, órbitas periódicas, ciclos límite y comportamiento caótico[cite: 16, 17]. Computacionalmente, se implementan y comparan tres métodos de integración numérica de paso fijo (**Euler**, **Heun** y **Runge–Kutta de cuarto orden**) mediante la función `odesolver`[cite: 19, 20, 21, 22, 44, 45].

Posteriormente, utilizando cálculo simbólico y numérico en MATLAB, se determinan los seis puntos de equilibrio del sistema y se evalúa su estabilidad local a partir de los valores propios de la matriz Jacobiana[cite: 26, 27, 28, 29]. Finalmente, se deducen analíticamente las cotas de localización en el octante no negativo $\mathbb{R}_{+,0}^3$ y se aplican tanto la linealización local como el **método directo de Lyapunov** para calcular los umbrales del parámetro de control $\rho_i$ (tratamiento de inmunoterapia) que garantizan la eliminación local y global de las células patológicas, validando ambos resultados mediante experimentación *in silico*[cite: 5, 6, 8, 9, 12, 13, 16, 36, 40].

**Palabras clave:** atractor caótico; bifurcación; conjuntos compactos invariantes; ecuaciones diferenciales ordinarias; estabilidad de Lyapunov; inmunoterapia; MATLAB; métodos numéricos; sistemas biológicos no lineales.

[Volver al contenido](#contenido)

---

<a id="objetivos"></a>
## Objetivos

### Objetivo general

Comparar métodos de diferencias finitas para solucionar sistemas de ecuaciones diferenciales ordinarias no lineales de primer orden que exhiben diferentes conjuntos compactos invariantes y comportamiento caótico, estableciendo condiciones analíticas y computacionales para el control y erradicación de poblaciones celulares patológicas[cite: 6, 16].

### Objetivos específicos

- Implementar y comparar los métodos numéricos de Euler, Heun y Runge–Kutta de cuarto orden (RK4) sobre el régimen de atractor caótico.
- Simular las soluciones temporales y las trayectorias tridimensionales en el espacio de fase para cinco regímenes dinámicos distintos determinados por el parámetro de bifurcación $a_{12}$.
- Calcular de forma simbólica y numérica los seis puntos de equilibrio del sistema y determinar su estabilidad local mediante los valores propios de la matriz Jacobiana[cite: 26, 27, 29].
- Demostrar matemáticamente la positividad e invariancia del sistema en el octante no negativo $\mathbb{R}_{+,0}^3$[cite: 5, 6].
- Determinar los conjuntos compactos invariantes $K_x$, $K_y$ y $K_z$ mediante funciones localizadoras y derivadas de Lie[cite: 9, 10, 11, 12].
- Deducir y simular las condiciones críticas de tratamiento $\rho_i$ para la **estabilidad asintótica local** y la **estabilidad asintótica global** (mediante el método directo de Lyapunov) del equilibrio libre de células patológicas[cite: 8, 12, 13, 36, 40].

[Volver al contenido](#contenido)

---

<a id="modelo"></a>
## Modelo matemático

El sistema biológico que describe la interacción entre tres poblaciones celulares se formula mediante tres ecuaciones diferenciales ordinarias no lineales acopladas de primer orden[cite: 16]:

$$\frac{dx}{dt} = r_1 x (1 - b_1 x) - a_{12} x y - a_{13} x z,$$

$$\frac{dy}{dt} = r_2 y (1 - b_2 y) - a_{21} x y,$$

$$\frac{dz}{dt} = (r_3 - a_{31}) x z - d_3 z + \rho_i(t),$$

donde todos los parámetros cinéticos son estrictamente positivos, las condiciones iniciales son no negativas $(x(0), y(0), z(0)) \in \mathbb{R}_{+,0}^3$ y el parámetro de control satisface $\rho_i \ge 0$[cite: 5, 16, 17].

### Variables de estado y parámetros

| Símbolo | Valor base | Descripción |
|---|:---:|---|
| $x(t)$ | Variable | Población de células anormales o patológicas (tumorales)[cite: 16] |
| $y(t)$ | Variable | Población de células normales o sanas[cite: 16] |
| $z(t)$ | Variable | Población de células efectoras del sistema inmune (citotóxicas o NK)[cite: 16] |
| $r_1$ | `1.0` | Tasa de crecimiento intrínseco de las células patológicas[cite: 16] |
| $b_1$ | `1.0` | Inverso de la capacidad de carga máxima de las células patológicas[cite: 16] |
| $a_{12}$ | `[0, 1]` | Fracción de células patológicas eliminadas por competencia con células sanas[cite: 16] |
| $a_{13}$ | `2.5` | Fracción de células patológicas eliminadas por acción de células efectoras[cite: 16] |
| $r_2$ | `0.6` | Tasa de crecimiento intrínseco de las células sanas[cite: 16] |
| $b_2$ | `1.0` | Inverso de la capacidad de carga máxima de las células sanas[cite: 16] |
| $a_{21}$ | `1.5` | Fracción de células sanas eliminadas o desplazadas por células patológicas[cite: 16] |
| $r_3$ | `3.0` | Tasa de reclutamiento máximo de células efectoras[cite: 16] |
| $a_{31}$ | `0.3` | Fracción de células efectoras inactivadas por las células patológicas[cite: 16] |
| $d_3$ | `0.4` | Tasa de mortalidad natural de las células efectoras[cite: 16] |
| $\rho_i(t)$ | Control | Parámetro de entrada correspondiente al tratamiento de inmunoterapia[cite: 16] |

[Volver al contenido](#contenido)

---

<a id="regimenes"></a>
## Regímenes dinámicos y parámetros

Al variar el parámetro de bifurcación $a_{12} \in [0,1]$ en ausencia de tratamiento ($\rho_i = 0$), el sistema exhibe cinco comportamientos cualitativos distintos configurados mediante la función `getdata(dynamics)`[cite: 16, 17, 19, 48]:

| Dinámica (`dynamics`) | Comportamiento observado | $a_{12}$ | Tiempo final (`tend`) | Condiciones iniciales $(x_0, y_0, z_0)$ |
|---|---|:---:|:---:|:---:|
| `'chaos'` | Atractor caótico[cite: 17, 48] | `1.00` | `1000`[cite: 17, 19] | `(0.555, 0.148, 0.172)`[cite: 17, 48] |
| `'orbit'` | Órbita periódica[cite: 17, 48] | `0.90` | `500`[cite: 17, 22] | `(0.254781, 0.115167, 0.600527)`[cite: 17, 48] |
| `'cycleI'` | Ciclo límite interno[cite: 17, 48] | `0.80` | `500`[cite: 17, 23] | `(0.2, 0.6, 0.1)`[cite: 17, 48] |
| `'cycleE'` | Ciclo límite externo[cite: 17, 48] | `0.75` | `500`[cite: 17, 24] | `(0.4, 0.2, 0.4)`[cite: 17, 48] |
| `'equilibrium'` | Convergencia a punto de equilibrio[cite: 17, 48] | `0.50` | `500`[cite: 17, 25] | `(0.5, 0.5, 0.4)`[cite: 17, 48] |

[Volver al contenido](#contenido)

---

<a id="analisis-equilibrios"></a>
## Análisis matemático: Positividad y equilibrios

### 1. Positividad e invariancia del sistema

Evaluando el campo vectorial sobre las fronteras del octante no negativo $\mathbb{R}_{+,0}^3 = \{(x,y,z) \in \mathbb{R}^3 \mid x \ge 0,\; y \ge 0,\; z \ge 0\}$[cite: 6, 17]:

$$\left.\dot{x}\right\vert{}_{x=0} = 0, \qquad \left.\dot{y}\right\vert{}_{y=0} = 0, \qquad \left.\dot{z}\right\vert{}_{z=0} = \rho_i \ge 0.$$

Por el criterio de De Leenheer y Aeyels, ninguna trayectoria con condiciones iniciales no negativas puede abandonar $\mathbb{R}_{+,0}^3$, garantizando que el dominio es **positivamente invariante** para todo $t \ge 0$ y manteniendo el sentido biológico de las poblaciones celulares[cite: 5, 6].

### 2. Puntos de equilibrio libres de células patológicas

Igualando las derivadas a cero y fijando $x = 0$ para buscar los estados libres de enfermedad[cite: 6]:

$$r_2 y (1 - b_2 y) = 0, \qquad -d_3 z + \rho_i = 0,$$

se obtienen dos puntos de equilibrio de interés biológico en la frontera libre de tumor[cite: 6, 7]:

$$E_0 = \left(x_0^*, y_0^*, z_0^*\right) = \left(0,\; 0,\; \frac{\rho_i}{d_3}\right),$$

$$E_1 = \left(x_1^*, y_1^*, z_1^*\right) = \left(0,\; \frac{1}{b_2},\; \frac{\rho_i}{d_3}\right).$$

El punto $E_0$ representa la extinción total de células patológicas y sanas, mientras que $E_1$ corresponde al **equilibrio saludable deseado**, donde las células sanas alcanzan su capacidad de carga $1/b_2$ y las células patológicas son erradicadas[cite: 7]. En total, la resolución simbólica con `solve` identifica **6 puntos de equilibrio** algebraicos en el sistema[cite: 26, 48].

[Volver al contenido](#contenido)

---

<a id="estabilidad-local"></a>
## Estabilidad asintótica local

La matriz Jacobiana del sistema está dada por[cite: 7, 29]:

$$J(x,y,z) = \begin{bmatrix} r_1 - 2 b_1 r_1 x - a_{12} y - a_{13} z & -a_{12} x & -a_{13} x \\ -a_{21} y & r_2 - 2 b_2 r_2 y - a_{21} x & 0 \\ (r_3 - a_{31}) z & 0 & (r_3 - a_{31}) x - d_3 \end{bmatrix}.$$

### Evaluación en $E_0 = (0, 0, \rho_i/d_3)$

Al evaluar $J$ en el origen libre de células sanas y tumorales se obtiene una matriz triangular inferior con valores propios[cite: 9]:

$$\lambda_1 = r_1 - \frac{a_{13}\rho_i}{d_3}, \qquad \lambda_2 = r_2 > 0, \qquad \lambda_3 = -d_3 < 0.$$

Dado que $\lambda_2 = r_2 > 0$ para cualquier combinación positiva de parámetros, el equilibrio $E_0$ es siempre un **punto silla inestable**[cite: 9].

### Evaluación en $E_1 = (0, 1/b_2, \rho_i/d_3)$

Al evaluar $J$ en el equilibrio saludable $E_1$ se obtiene[cite: 7, 8]:

$$J(E_1) = \begin{bmatrix} r_1 - \dfrac{a_{12}}{b_2} - \dfrac{a_{13}\rho_i}{d_3} & 0 & 0 \\ -\dfrac{a_{21}}{b_2} & -r_2 & 0 \\ \dfrac{\rho_i}{d_3}(r_3 - a_{31}) & 0 & -d_3 \end{bmatrix},$$

cuyos valores propios son $\lambda_1 = r_1 - \frac{a_{12}}{b_2} - \frac{a_{13}\rho_i}{d_3}$, $\lambda_2 = -r_2 < 0$ y $\lambda_3 = -d_3 < 0$[cite: 8]. Para garantizar $\lambda_1 < 0$ y lograr la **estabilidad asintótica local** de $E_1$, el parámetro de inmunoterapia debe satisfacer[cite: 8]:

$$\rho_i > \frac{d_3}{a_{13}}\left(r_1 - \frac{a_{12}}{b_2}\right).$$

Si $\rho_i < \frac{d_3}{a_{13}}\left(r_1 - \frac{a_{12}}{b_2}\right)$, las células patológicas persisten en el sistema[cite: 8, 9]. En el cuaderno de MATLAB se verifica esta condición utilizando un margen del 1 % sobre el umbral crítico (`rho_L = (d3/a13)*(r1 - a12/b2)*1.01`, y `+ 0.01` para el caso caótico donde $r_1 - a_{12}/b_2 = 0$).

[Volver al contenido](#contenido)

---

<a id="estabilidad-global"></a>
## Conjuntos compactos invariantes y estabilidad global

### 1. Localización de conjuntos compactos invariantes

Mediante el método de localización de conjuntos compactos invariantes usando funciones localizadoras $h_1 = x$, $h_2 = y$ y $h_3 = z$ y sus respectivas derivadas de Lie $L_f h_i = \dot{x}_i$, se determinan las cotas superiores e inferiores de las trayectorias en $\mathbb{R}_{+,0}^3$ (bajo la condición $r_3 > a_{31}$)[cite: 9, 10, 11, 12]:

$$K_x = \left\{ 0 \le x(t) \le \frac{1}{b_1} \right\}, \qquad K_y = \left\{ 0 \le y(t) \le \frac{1}{b_2} \right\}, \qquad K_z = \left\{ z(t) \ge \frac{\rho_i}{d_3} \right\}.$$

### 2. Método directo de Lyapunov

Para garantizar la erradicación global de las células patológicas desde cualquier condición inicial positiva, se propone la función candidata de Lyapunov definida positiva y radialmente desacotada en $\mathbb{R}_{+,0}^3$[cite: 12]:

$$V(x) = x.$$

Calculando su derivada temporal a lo largo de las trayectorias del sistema y acotando sobre el dominio $K_x \cap K_y \cap K_z$ ($x \ge 0$, $y \ge 0$, $z \ge \rho_i/d_3$)[cite: 12, 13]:

$$\dot{V} = \dot{x} = \left(r_1 - b_1 r_1 x - a_{12} y - a_{13} z\right) x \le \left(r_1 - \frac{a_{13}}{d_3}\rho_i\right) x.$$

Para asegurar que $\dot{V} < 0$ para todo $x > 0$ (derivada definida negativa), se obtiene la **condición suficiente de eliminación global**[cite: 13]:

$$\rho_i > \frac{r_1 d_3}{a_{13}}.$$

Con los parámetros del modelo ($r_1 = 1$, $d_3 = 0.4$, $a_{13} = 2.5$), el umbral global es $\rho_{\text{global}} = 0.16$[cite: 13, 16]. Al aplicar $\rho_i = 1.01 \cdot \frac{r_1 d_3}{a_{13}} = 0.1616$, todas las trayectorias convergen globalmente al equilibrio saludable $(x^*, y^*, z^*) = (0, 1, 0.4040)$ independientemente de la dinámica original del sistema.

[Volver al contenido](#contenido)

---

<a id="metodos"></a>
## Métodos numéricos

La función `odesolver` implementa tres esquemas de integración explícitos de paso fijo (`dt = 1e-3`, y `dt = 1e-4` en equilibrio) para resolver $\dot{\mathbf{X}} = \mathbf{F}(\mathbf{X}, \rho_i(t))$[cite: 19, 39, 44, 45]:

### 1. Método de Euler explícito (`'euler'`)

$$\mathbf{X}_{n+1} = \mathbf{X}_n + h\,\mathbf{F}(\mathbf{X}_n, \rho_n).$$

### 2. Método de Heun o Euler mejorado (`'heun'`)

$$\widetilde{\mathbf{X}}_{n+1} = \mathbf{X}_n + h\,\mathbf{F}(\mathbf{X}_n, \rho_n), \qquad \mathbf{X}_{n+1} = \mathbf{X}_n + \frac{h}{2}\left[\mathbf{F}(\mathbf{X}_n, \rho_n) + \mathbf{F}(\widetilde{\mathbf{X}}_{n+1}, \rho_n)\right].$$

### 3. Método de Runge–Kutta de cuarto orden (`'runge-kutta'`)

$$\begin{aligned} \mathbf{k}_1 &= \mathbf{F}(\mathbf{X}_n, \rho_n), \\ \mathbf{k}_2 &= \mathbf{F}\left(\mathbf{X}_n + \frac{h}{2}\mathbf{k}_1,\; \rho_n\right), \\ \mathbf{k}_3 &= \mathbf{F}\left(\mathbf{X}_n + \frac{h}{2}\mathbf{k}_2,\; \rho_n\right), \\ \mathbf{k}_4 &= \mathbf{F}\left(\mathbf{X}_n + h\mathbf{k}_3,\; \rho_n\right), \\ \mathbf{X}_{n+1} &= \mathbf{X}_n + \frac{h}{6}\left(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4\right). \end{aligned}$$

Al finalizar la integración, `odesolver` submuestrea la salida con un factor `h = 10` (`1/(h*dt)`) para optimizar el renderizado de las series temporales y trayectorias en $\mathbb{R}^3$[cite: 45].

[Volver al contenido](#contenido)

---

<a id="flujo"></a>
## Flujo de trabajo computacional

```mermaid
flowchart LR
    A[getdata: Selección de dinámica] --> B[Parámetros y CI]
    B --> C{Comparación de solvers}
    C --> D1[Euler]
    C --> D2[Heun]
    C --> D3[Runge–Kutta RK4]
    D1 --> E[plotEDOs: Series de tiempo y espacio de fase 3D]
    D2 --> E
    D3 --> E
    B --> F[Análisis simbólico: solve y jacobian]
    F --> G[equilibria: 6 equilibrios y valores propios]
    G --> H[plotCCI: Trayectorias y puntos de equilibrio]
    G --> I1[Control local: rho_L > d3/a13 * r1 - a12/b2]
    G --> I2[Control global Lyapunov: rho_L > r1*d3 / a13]
    I1 --> J[Simulación con inmunoterapia constante]
    I2 --> J
    J --> K[Exportación de figuras vectoriales PDF]