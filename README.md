[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=emanuelalejandro9-web/MMP3)

# Práctica 2 — Caos en sistemas biológicos y control de células patológicas

Repositorio académico de la asignatura **Modelado Matemático**, correspondiente a la **Maestría en Ciencias de la Ingeniería** del Tecnológico Nacional de México / Instituto Tecnológico de Tijuana. La práctica integra modelado mecanicista de sistemas dinámicos no lineales en tres dimensiones, comparación de métodos numéricos de diferencias finitas, cálculo simbólico y numérico de puntos de equilibrio, análisis de estabilidad asintótica local mediante la matriz Jacobiana, localización de conjuntos compactos invariantes y análisis de estabilidad asintótica global mediante el método directo de Lyapunov.

El caso de estudio utiliza el modelo tridimensional de cáncer propuesto por **Itik y Banks (2010)** para representar la dinámica no lineal entre tres poblaciones celulares: células patológicas o tumorales, células normales o sanas y células efectoras del sistema inmune (citotóxicas o *Natural Killer*), evaluando tanto su comportamiento caótico natural como la dosis crítica de inmunoterapia $\rho_i$ necesaria para erradicar la población patológica.

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
3. [Fuente y preparación de los datos](#datos)
4. [Modelo matemático](#modelo)
5. [Regímenes dinámicos y parámetros](#regimenes)
6. [Análisis matemático: Positividad y equilibrios](#analisis-equilibrios)
7. [Estabilidad asintótica local](#estabilidad-local)
8. [Conjuntos compactos invariantes y estabilidad global](#estabilidad-global)
9. [Métodos numéricos](#metodos)
10. [Flujo de trabajo computacional](#flujo)
11. [Actividades y secciones del cuaderno MATLAB](#actividades)
12. [Funciones auxiliares](#funciones)
13. [Estructura del repositorio](#estructura)
14. [Requisitos](#requisitos)
15. [Reproducibilidad](#reproducibilidad)
16. [Referencias](#referencias)

---

<a id="resumen"></a>
## Resumen

Esta práctica estudia la dinámica no lineal y el control de un sistema biológico tridimensional descrito por ecuaciones diferenciales ordinarias de primer orden que exhibe bifurcaciones, órbitas periódicas, ciclos límite y comportamiento caótico. Computacionalmente, se implementan y comparan tres métodos de integración numérica de paso fijo (**Euler**, **Heun** y **Runge–Kutta de cuarto orden**) mediante la función `odesolver`.

Posteriormente, utilizando cálculo simbólico y numérico en MATLAB, se determinan los seis puntos de equilibrio del sistema y se evalúa su estabilidad local a partir de los valores propios de la matriz Jacobiana. Finalmente, se deducen analíticamente las cotas de localización en el octante no negativo $\mathbb{R}_{+,0}^{3}$ y se aplican tanto la linealización local como el **método directo de Lyapunov** para calcular los umbrales del parámetro de control $\rho_i$ (tratamiento de inmunoterapia) que garantizan la eliminación local y global de las células patológicas, validando ambos resultados mediante experimentación *in silico*.

**Palabras clave:** atractor caótico; bifurcación; conjuntos compactos invariantes; ecuaciones diferenciales ordinarias; estabilidad de Lyapunov; inmunoterapia; MATLAB; métodos numéricos; sistemas biológicos no lineales.

[Volver al contenido](#contenido)

---

<a id="objetivos"></a>
## Objetivos

### Objetivo general

Comparar métodos de diferencias finitas para solucionar sistemas de ecuaciones diferenciales ordinarias no lineales de primer orden que exhiben diferentes conjuntos compactos invariantes y comportamiento caótico, estableciendo condiciones analíticas y computacionales para el control y erradicación de poblaciones celulares patológicas.

### Objetivos específicos

- Implementar y comparar los métodos numéricos de Euler, Heun y Runge–Kutta de cuarto orden (RK4) sobre el régimen de atractor caótico.
- Simular las soluciones temporales y las trayectorias tridimensionales en el espacio de fase para cinco regímenes dinámicos distintos determinados por el parámetro de bifurcación $a_{12}$.
- Calcular de forma simbólica y numérica los seis puntos de equilibrio del sistema y determinar su estabilidad local mediante los valores propios de la matriz Jacobiana.
- Demostrar matemáticamente la positividad e invariancia del sistema en el octante no negativo $\mathbb{R}_{+,0}^{3}$.
- Determinar los conjuntos compactos invariantes $K_x$, $K_y$ y $K_z$ mediante funciones localizadoras y derivadas de Lie.
- Deducir y simular las condiciones críticas de tratamiento $\rho_i$ para la **estabilidad asintótica local** y la **estabilidad asintótica global** (mediante el método directo de Lyapunov) del equilibrio libre de células patológicas.

[Volver al contenido](#contenido)

---

<a id="datos"></a>
## Fuente y preparación de los datos

A diferencia de un problema inverso con series experimentales externas, en esta práctica los parámetros cinéticos base y las condiciones iniciales provienen directamente de la parametrización reportada por **Itik y Banks (2010)** en *Chaos in a three-dimensional cancer model* [2]. Todos los valores están adimensionalizados en el octante no negativo $\mathbb{R}_{+,0}^{3}$.

### Inicialización y estructuración en MATLAB

La preparación de los parámetros y condiciones iniciales para cada experimento numérico se centraliza mediante la función auxiliar `getdata(dynamics)`, la cual empaqueta el vector de parámetros `par` y el vector de estado inicial `f0` dentro de una estructura `D`:

```matlab
par = [r1, b1, a12, a13, r2, b2, a21, r3, a31, d3];
D.par = par;
D.f0  = f0;