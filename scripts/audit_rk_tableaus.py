#!/usr/bin/env python3
"""Audit the order conditions of the FOODIE Runge-Kutta tableaus.

The coefficients are parsed from the Fortran sources (decimal literals or ``real(p_I8P)/real(q_I8P)`` rationals) and the
order conditions of all the rooted trees up to the claimed order are evaluated in 60 digits arithmetic: a tableau passes if
all the residuals are below the double precision exactness threshold and the abscissae are the row sums of the matrix.
The low storage (Williamson 2N) schemes are converted to their Butcher form. Only the Python standard library is used.

Exit status is the number of failed tableaus.
"""

from __future__ import annotations

import re
import sys
from decimal import Decimal, getcontext
from functools import lru_cache
from pathlib import Path

getcontext().prec = 60

LIB = Path(__file__).resolve().parent.parent / "src" / "lib"
TOLERANCE = Decimal("1e-14")
NUM = (
    r"(-?\s*(?:real\(\s*-?\d+_I8P\s*,\s*kind=R_P\)\s*/\s*real\(\s*\d+_I8P\s*,\s*kind=R_P\)"
    r"|-?\d+\.\d*(?:[eE][-+]?\d+)?_R_P(?:\s*/\s*\d+\.\d*_R_P)?))"
)
# embedded pairs: (order of the solution advanced, beta(:,1); order of the error estimator, beta(:,2))
EMBEDDED_ORDERS = {
    "runge_kutta_emd_stages_2_order_2": (2, 1),
    "runge_kutta_emd_stages_6_order_5": (5, 4),
    "runge_kutta_emd_stages_7_order_4": (5, 4),
    "runge_kutta_emd_stages_9_order_6": (6, 5),
    "runge_kutta_emd_stages_17_order_10": (10, 8),
}

Tree = tuple  # rooted tree: sorted tuple of its children trees
Matrix = list[list[Decimal]]


def value(token: str) -> Decimal:
    """Return the value of a Fortran real literal or of a ``real(p_I8P)/real(q_I8P)`` rational.

    Parameters
    ----------
    token : str
        Fortran expression.

    Returns
    -------
    Decimal
        Value of the expression.
    """
    token = token.replace(" ", "")
    negative = token.startswith("-")
    token = token.lstrip("-")
    rational = re.match(r"real\((-?\d+)_I8P,kind=R_P\)/real\((\d+)_I8P,kind=R_P\)", token)
    if rational:
        result = Decimal(rational.group(1)) / Decimal(rational.group(2))
    else:
        parts = token.replace("_R_P", "").split("/")
        result = Decimal(parts[0]) / (Decimal(parts[1]) if len(parts) > 1 else Decimal(1))
    return -result if negative else result


def schemes_coefficients(path: Path) -> dict[str, dict[tuple[str, tuple[int, ...]], Decimal]]:
    """Return the coefficients of all the schemes defined in a FOODIE integrator source.

    Parameters
    ----------
    path : Path
        Fortran source file.

    Returns
    -------
    dict
        For each scheme name the map ``(array name, indexes) -> value``.
    """
    source = path.read_text()
    schemes = {}
    for block in re.split(r"\n\s*case\('", source)[1:]:
        name = block.split("'")[0]
        body = block.split("case default")[0]
        coefficients = {}
        for match in re.finditer(r"self%(\w+)\(\s*([\d\s,]+)\)\s*=\s*" + NUM, body):
            key = (match.group(1), tuple(int(i) for i in match.group(2).split(",")))
            coefficients[key] = value(match.group(3))
        if coefficients:
            schemes[name] = coefficients
    return schemes


@lru_cache(None)
def trees(order: int) -> list[Tree]:
    """Return all the rooted trees with the given number of nodes.

    Parameters
    ----------
    order : int
        Number of nodes.

    Returns
    -------
    list of Tree
        Rooted trees, each one as the sorted tuple of its children.
    """
    if order == 1:
        return [()]
    found = set()

    def forests(nodes: int, bound: tuple[int, str] | None):
        if nodes == 0:
            yield ()
            return
        for k in range(nodes, 0, -1):
            for tree in trees(k):
                key = (k, repr(tree))
                if bound is not None and key > bound:
                    continue
                for rest in forests(nodes - k, key):
                    yield (tree, *rest)

    for forest in forests(order - 1, None):
        found.add(tuple(sorted(forest, key=repr)))
    return sorted(found, key=repr)


def tree_order(tree: Tree) -> int:
    """Return the number of nodes of a tree."""
    return 1 + sum(tree_order(child) for child in tree)


def tree_density(tree: Tree) -> Decimal:
    """Return the density (gamma) of a tree."""
    density = Decimal(tree_order(tree))
    for child in tree:
        density *= tree_density(child)
    return density


def max_residual(a: Matrix, b: list[Decimal], order: int) -> Decimal:
    """Return the max residual of the order conditions of all the trees up to the given order.

    Parameters
    ----------
    a : list of list of Decimal
        Butcher matrix.
    b : list of Decimal
        Weights.
    order : int
        Order checked.

    Returns
    -------
    Decimal
        Max absolute residual ``b . Phi(t) - 1 / gamma(t)``.
    """
    stages = len(b)
    cache: dict[Tree, list[Decimal]] = {}

    def weights(tree: Tree) -> list[Decimal]:
        if tree not in cache:
            vector = [Decimal(1)] * stages
            for child in tree:
                w = weights(child)
                vector = [vector[i] * sum(a[i][j] * w[j] for j in range(stages)) for i in range(stages)]
            cache[tree] = vector
        return cache[tree]

    return max(
        abs(sum(bi * wi for bi, wi in zip(b, weights(tree), strict=True)) - 1 / tree_density(tree))
        for p in range(1, order + 1)
        for tree in trees(p)
    )


def low_storage_to_butcher(a_ls: list[Decimal], b_ls: list[Decimal]) -> tuple[Matrix, list[Decimal]]:
    """Convert a Williamson 2N low storage scheme into its Butcher form.

    Parameters
    ----------
    a_ls, b_ls : list of Decimal
        Low storage coefficients, ``dU_j = A_j dU_{j-1} + Dt R(U_{j-1})``, ``U_j = U_{j-1} + B_j dU_j``.

    Returns
    -------
    tuple
        Butcher matrix and weights.
    """
    stages = len(a_ls)
    u = [Decimal(0)] * stages
    du = [Decimal(0)] * stages
    a = []
    for j in range(stages):
        a.append(list(u))
        du = [a_ls[j] * x for x in du]
        du[j] += 1
        u = [ui + b_ls[j] * di for ui, di in zip(u, du, strict=True)]
    return a, u


def claimed_order(name: str) -> int:
    """Return the order encoded in the scheme name."""
    return int(re.search(r"order_(\d+)", name).group(1))


def check(label: str, a: Matrix, b: list[Decimal], c: list[Decimal], order: int) -> bool:
    """Check a tableau and print the outcome.

    Parameters
    ----------
    label : str
        Tableau label.
    a : list of list of Decimal
        Butcher matrix.
    b, c : list of Decimal
        Weights and abscissae (as coded).
    order : int
        Order checked.

    Returns
    -------
    bool
        True if the tableau passes.
    """
    residual = max_residual(a, b, order)
    c_error = max(abs(ci - sum(row)) for ci, row in zip(c, a, strict=True))
    passed = residual <= TOLERANCE and c_error <= TOLERANCE
    print(
        f"{'PASS' if passed else 'FAIL'} {label:48s} order {order:2d}: "
        f"max |order residual| = {float(residual):.2e}, max |c - row sum| = {float(c_error):.2e}"
    )
    return passed


def main() -> int:
    """Audit all the Runge-Kutta tableaus.

    Returns
    -------
    int
        Number of failed tableaus.
    """
    failures = 0
    for name, co in schemes_coefficients(LIB / "foodie_integrator_runge_kutta_ssp.F90").items():
        s = max(k[1][0] for k in co if k[0] == "beta")
        a = [[co.get(("alph", (i, j)), Decimal(0)) for j in range(1, s + 1)] for i in range(1, s + 1)]
        b = [co.get(("beta", (i,)), Decimal(0)) for i in range(1, s + 1)]
        c = [co.get(("gamm", (i,)), Decimal(0)) for i in range(1, s + 1)]
        failures += not check(name, a, b, c, claimed_order(name))
    for name, co in schemes_coefficients(LIB / "foodie_integrator_runge_kutta_embedded.F90").items():
        s = max(k[1][0] for k in co if k[0] == "beta")
        a = [[co.get(("alph", (i, j)), Decimal(0)) for j in range(1, s + 1)] for i in range(1, s + 1)]
        c = [co.get(("gamm", (i,)), Decimal(0)) for i in range(1, s + 1)]
        for column, order in enumerate(EMBEDDED_ORDERS[name], start=1):
            b = [co.get(("beta", (i, column)), Decimal(0)) for i in range(1, s + 1)]
            failures += not check(f"{name} [beta(:,{column})]", a, b, c, order)
    for name, co in schemes_coefficients(LIB / "foodie_integrator_runge_kutta_low_storage.F90").items():
        if not any(k[0] == "A" for k in co):
            continue
        s = max(k[1][0] for k in co if k[0] == "A")
        a_ls = [co.get(("A", (i,)), Decimal(0)) for i in range(1, s + 1)]
        b_ls = [co.get(("B", (i,)), Decimal(0)) for i in range(1, s + 1)]
        a, b = low_storage_to_butcher(a_ls, b_ls)
        failures += not check(name, a, b, [co.get(("C", (i,)), Decimal(0)) for i in range(1, s + 1)], claimed_order(name))
    print(f"{failures} failures")
    return failures


if __name__ == "__main__":
    sys.exit(main())
