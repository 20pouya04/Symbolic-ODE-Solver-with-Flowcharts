# Symbolic ODE Solver with Auto-Generated Decision Flowcharts

A Python tool (**SymPy + Graphviz**) that symbolically solves **second-order
linear ODEs** — both **constant-coefficient** and **Cauchy–Euler** forms —
and automatically renders the corresponding solution-method flowchart.

**Status:** Course Project, Differential Equations (2024)

## Background

This project comes from my **Differential Equations** course during
university. It was originally a script I wrote to automate solving these
standard ODE forms (rather than doing the characteristic-equation case
analysis by hand) and to draw the accompanying method flowchart.

## Supported equation types

**1. Constant-coefficient equations**

```
y'' + K*y' + L*y = 0
```
Characteristic equation: `m^2 + K*m + L = 0`

**2. Cauchy-Euler (equidimensional) equations**

```
a2*x^2*y'' + a1*x*y' + a0*y = 0
```
Characteristic equation: `a2*m^2 + (a1 - a2)*m + a0 = 0`

In both cases, the general solution depends on the discriminant (`Δ`) of
the characteristic equation:

| Discriminant | Roots | General solution form |
|---|---|---|
| `Δ > 0` | Two distinct real roots `m1, m2` | Exponential / power terms for each root |
| `Δ = 0` | One repeated real root `m` | Exponential / power term, plus a second term multiplied by `x` (or `ln x`) |
| `Δ < 0` | Complex conjugate roots `α ± βi` | Oscillatory solution with `sin`/`cos` |

## Files

| File | Description |
|---|---|
| `ode_solver.py` | Clean, documented Python module — solves both ODE types and generates their flowcharts. |
| `Flowchart.ipynb` | Original Jupyter notebook version of the project. |
| `LICENSE` | MIT License. |

## Requirements

- Python 3.8+
- [sympy](https://www.sympy.org/) — symbolic math
- [graphviz](https://pypi.org/project/graphviz/) (Python package) — for flowchart generation
- [Graphviz](https://graphviz.org/download/) (system package, provides the `dot` executable) — required only if you want to render flowcharts

Install the Python dependencies with:

```bash
pip install -r requirements.txt
```

On most Linux systems, the Graphviz system package can be installed with:

```bash
sudo apt-get install graphviz
```

On Windows/macOS, download it from [graphviz.org](https://graphviz.org/download/) and make sure the `dot` executable is on your `PATH`.

## Usage

Run the script directly to see both example ODEs solved, and to generate
their flowcharts as PNG files in the current directory:

```bash
python ode_solver.py
```

### Using your own equations

```python
from ode_solver import solve_constant_coefficient_ode, solve_cauchy_euler_ode

# y'' - 4y' + 4y = 0
result = solve_constant_coefficient_ode(K=-4, L=4)
print(result["y_gh"])

# x^2 y'' + 3x y' + y = 0
result = solve_cauchy_euler_ode(a2=1, a1=3, a0=1)
print(result["y_gh"])
```

Each call returns a dictionary containing the discriminant, the roots (or
`alpha`/`beta` for the complex case), and the symbolic general solution
`y_gh`.

To generate the flowcharts as PNG images:

```python
from ode_solver import generate_constant_coefficient_flowchart, generate_cauchy_euler_flowchart

generate_constant_coefficient_flowchart("my_flowchart_1")
generate_cauchy_euler_flowchart("my_flowchart_2")
```

Alternatively, open and run `Flowchart.ipynb` in Jupyter for the original
notebook walkthrough.

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE)
for details.
