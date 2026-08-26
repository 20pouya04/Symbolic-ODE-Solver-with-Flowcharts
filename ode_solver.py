"""
Second-Order Linear ODE Solver & Flowchart Generator
======================================================

Solves two classic families of second-order linear homogeneous ODEs by
building and solving their characteristic equation, and generates a
flowchart (using Graphviz) that visualizes the solution method.

Supported equation types
-------------------------
1. Constant-coefficient equations:
       y'' + K*y' + L*y = 0
   Characteristic equation:  m^2 + K*m + L = 0

2. Cauchy-Euler (equidimensional) equations:
       a2*x^2*y'' + a1*x*y' + a0*y = 0
   Characteristic equation:  a2*m^2 + (a1 - a2)*m + a0 = 0

In both cases, the nature of the general solution depends on the
discriminant (Delta) of the characteristic equation:

    Delta > 0  -> two distinct real roots      -> exponential / power solution
    Delta = 0  -> one repeated real root       -> solution with an extra x / ln(x) term
    Delta < 0  -> complex conjugate roots      -> oscillatory (sin/cos) solution

Usage
-----
    python ode_solver.py

Or import the functions to solve your own coefficients:

    from ode_solver import solve_constant_coefficient_ode, solve_cauchy_euler_ode

    solve_constant_coefficient_ode(K=-4, L=4)
    solve_cauchy_euler_ode(a2=1, a1=3, a0=1)

Flowchart generation requires the Graphviz system package (the `dot`
executable) in addition to the `graphviz` Python package. If Graphviz is
not installed, the solver functions above still work — only the
`generate_*_flowchart` functions need it.
"""

from sympy import Function, Eq, Derivative, symbols, solve, exp, cos, sin, ln

try:
    from graphviz import Digraph

    GRAPHVIZ_AVAILABLE = True
except ImportError:
    GRAPHVIZ_AVAILABLE = False


# ---------------------------------------------------------------------------
# Solvers
# ---------------------------------------------------------------------------

def solve_constant_coefficient_ode(K, L):
    """
    Solve  y'' + K*y' + L*y = 0  by finding the roots of its characteristic
    equation  m^2 + K*m + L = 0  and building the general solution y_gh.

    Returns a dict with the discriminant, the roots (or alpha/beta for the
    complex case), and the symbolic general solution `y_gh`.
    """
    x = symbols("x")
    C1, C2 = symbols("C1 C2")
    r = symbols("r")

    char_eq = Eq(r**2 + K * r + L, 0)
    roots = solve(char_eq, r)

    delta = K**2 - 4 * L

    if delta > 0:
        m1, m2 = roots
        y_gh = C1 * exp(m1 * x) + C2 * exp(m2 * x)
        return {"case": "distinct_real_roots", "delta": delta, "m1": m1, "m2": m2, "y_gh": y_gh}

    elif delta == 0:
        m = roots[0]
        y_gh = (C1 + C2 * x) * exp(m * x)
        return {"case": "repeated_real_root", "delta": delta, "m": m, "y_gh": y_gh}

    else:
        alpha = roots[0].as_real_imag()[0]
        beta = abs(roots[0].as_real_imag()[1])
        y_gh = C1 * exp(alpha * x) * cos(beta * x) + C2 * exp(alpha * x) * sin(beta * x)
        return {"case": "complex_roots", "delta": delta, "alpha": alpha, "beta": beta, "y_gh": y_gh}


def solve_cauchy_euler_ode(a2, a1, a0):
    """
    Solve the Cauchy-Euler equation  a2*x^2*y'' + a1*x*y' + a0*y = 0  by
    finding the roots of its characteristic equation
    a2*m^2 + (a1 - a2)*m + a0 = 0  and building the general solution y_gh.

    Returns a dict with the discriminant, the roots (or alpha/beta for the
    complex case), and the symbolic general solution `y_gh`.
    """
    x = symbols("x")
    C1, C2 = symbols("C1 C2")
    r = symbols("r")

    char_eq = Eq(a2 * r**2 + (a1 - a2) * r + a0, 0)
    roots = solve(char_eq, r)

    delta = (a1 - a2) ** 2 - 4 * a2 * a0

    if delta > 0:
        m1, m2 = roots
        y_gh = C1 * x**m1 + C2 * x**m2
        return {"case": "distinct_real_roots", "delta": delta, "m1": m1, "m2": m2, "y_gh": y_gh}

    elif delta == 0:
        m = roots[0]
        y_gh = (C1 + C2 * ln(x)) * x**m
        return {"case": "repeated_real_root", "delta": delta, "m": m, "y_gh": y_gh}

    else:
        alpha = roots[0].as_real_imag()[0]
        beta = abs(roots[0].as_real_imag()[1])
        y_gh = x**alpha * (C1 * cos(beta * ln(x)) + C2 * sin(beta * ln(x)))
        return {"case": "complex_roots", "delta": delta, "alpha": alpha, "beta": beta, "y_gh": y_gh}


# ---------------------------------------------------------------------------
# Flowcharts
# ---------------------------------------------------------------------------

def _build_solver_flowchart(title, node_a_label, label_left, label_right, m_formula_g, m_formula_h):
    """Shared flowchart structure for both solver types."""
    dot = Digraph(engine="dot")
    dot.attr(size="10,4", ratio="fill")

    dot.node("Title", title, shape="plaintext", fontsize="20")
    dot.node("A", node_a_label, shape="rectangle")
    dot.node("B", "Δ = b^2 - 4ac", shape="rectangle")
    dot.node("C", "Δ > 0", shape="triangle")
    dot.node("D", "Δ = 0", shape="triangle")
    dot.node("E", "β = √Δ / 2a\nα = -b / 2a", shape="rectangle")
    dot.node("F", "END", shape="circle")
    dot.node("G", m_formula_g, shape="rectangle")
    dot.node("H", m_formula_h, shape="rectangle")

    dot.edge("A", "B")
    dot.edge("B", "C")
    dot.edge("C", "D", label="NO")
    dot.edge("D", "E", label="NO")
    dot.edge("E", "F")
    dot.edge("C", "G", label="yes")
    dot.edge("D", "H", label="yes")
    dot.edge("G", "F")
    dot.edge("H", "F")

    dot.node("label_right", "y\u2033 \u2192 m^2\n y\u2032 \u2192 m\n y \u2192 1", shape="plaintext", width="0.1", height="0.1")
    dot.node("label_left", label_left, shape="plaintext", width="0.1", height="0.1")

    return dot


def generate_constant_coefficient_flowchart(output_path="flowchart_constant_coefficient"):
    """
    Render a flowchart (PNG) illustrating the solution method for
    y'' + K*y' + L*y = 0. Requires Graphviz (`dot`) to be installed.
    """
    if not GRAPHVIZ_AVAILABLE:
        raise RuntimeError("The 'graphviz' Python package is not installed. Run: pip install graphviz")

    dot = _build_solver_flowchart(
        title="Characteristic Equation Solver",
        node_a_label="m^2 + Km + L = 0",
        label_left="y\u2033 + Ky\u2032 + Ly = 0",
        label_right=None,
        m_formula_g="m1, m2 =\n(-b ± √Δ) / 2\n\ny_gh = C1 e^(m1 x) + C2 e^(m2 x)",
        m_formula_h="m =\n- b / 2a\n\ny_gh = C1 e^(mx) + C2 x e^(mx)",
    )
    dot.node("E", "β = √Δ / 2a\nα = -b / 2a\n\ny_gh = e^(α x) [ C₁ cos(β x) + C₂ sin(β x) ]", shape="rectangle")
    dot.render(output_path, format="png", cleanup=True)
    return f"{output_path}.png"


def generate_cauchy_euler_flowchart(output_path="flowchart_cauchy_euler"):
    """
    Render a flowchart (PNG) illustrating the solution method for
    a2*x^2*y'' + a1*x*y' + a0*y = 0. Requires Graphviz (`dot`) to be installed.
    """
    if not GRAPHVIZ_AVAILABLE:
        raise RuntimeError("The 'graphviz' Python package is not installed. Run: pip install graphviz")

    dot = _build_solver_flowchart(
        title="Cauchy-Euler Equation Solver",
        node_a_label="a1m^2 + (a1+a2)m + a0 = 0",
        label_left="a2x^2y\u2033 + a1xy\u2032 + a0y = 0",
        label_right=None,
        m_formula_g="m1, m2 =\n(-b ± √Δ) / 2a\n\ny_gh = C1 x^(m1) + C2 x^(m2)",
        m_formula_h="m =\n- b / 2a\n\ny_gh = C1 x^(m) + C2 x^(m)Lnx",
    )
    dot.node("E", "β = √Δ / 2a\nα = -b / 2a\n\ny_gh = x^α [ C₁ cos(β Lnx) + C₂ sin(β Lnx) ]", shape="rectangle")
    dot.render(output_path, format="png", cleanup=True)
    return f"{output_path}.png"


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------

def main():
    print("=== Constant-coefficient ODE: y'' - 4y' + 4y = 0 ===")
    result = solve_constant_coefficient_ode(K=-4, L=4)
    for key, value in result.items():
        print(f"{key}: {value}")

    print("\n=== Cauchy-Euler ODE: x^2 y'' + 3x y' + y = 0 ===")
    result = solve_cauchy_euler_ode(a2=1, a1=3, a0=1)
    for key, value in result.items():
        print(f"{key}: {value}")

    if GRAPHVIZ_AVAILABLE:
        print("\nGenerating flowcharts...")
        try:
            path1 = generate_constant_coefficient_flowchart()
            path2 = generate_cauchy_euler_flowchart()
            print(f"Saved: {path1}")
            print(f"Saved: {path2}")
        except Exception as exc:
            print(f"Could not render flowcharts (is the Graphviz 'dot' executable installed?): {exc}")
    else:
        print("\n'graphviz' package not installed — skipping flowchart generation.")


if __name__ == "__main__":
    main()
