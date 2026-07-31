#!/usr/bin/env python3
"""Dependency-free sanity check for the registered MLS mathematics.

This does not test CubeCL. It independently checks the Wendland weighting,
Jacobi eigensolver, normal orientation and five projection iterations on a
plane. It is useful before spending time compiling the ROCm dependency tree.
"""

from __future__ import annotations

import math


def rotate(a: list[list[float]], v: list[list[float]], p: int, q: int) -> None:
    apq = a[p][q]
    if abs(apq) <= 1e-12:
        return
    tau = (a[q][q] - a[p][p]) / (2.0 * apq)
    t = (1.0 if tau >= 0 else -1.0) / (abs(tau) + math.sqrt(1.0 + tau * tau))
    c = 1.0 / math.sqrt(1.0 + t * t)
    s = t * c
    app, aqq = a[p][p], a[q][q]
    a[p][p] = app - t * apq
    a[q][q] = aqq + t * apq
    a[p][q] = a[q][p] = 0.0
    for k in range(3):
        if k not in (p, q):
            akp, akq = a[k][p], a[k][q]
            a[k][p] = a[p][k] = c * akp - s * akq
            a[k][q] = a[q][k] = s * akp + c * akq
    for row in range(3):
        vip, viq = v[row][p], v[row][q]
        v[row][p] = c * vip - s * viq
        v[row][q] = s * vip + c * viq


def smallest_eigenvector(cov: list[list[float]]) -> tuple[float, float, float]:
    a = [row[:] for row in cov]
    v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    for _ in range(8):
        rotate(a, v, 0, 1)
        rotate(a, v, 0, 2)
        rotate(a, v, 1, 2)
    col = min(range(3), key=lambda i: a[i][i])
    n = [v[0][col], v[1][col], v[2][col]]
    length = math.sqrt(sum(x * x for x in n))
    n = [x / length for x in n]
    if sum(n) < 0:
        n = [-x for x in n]
    return tuple(n)


def weight(d2: float, radius: float) -> float:
    q = math.sqrt(d2) / radius
    if q >= 1.0:
        return 0.0
    x = 1.0 - q
    return x**4 * (4.0 * q + 1.0)


def project_once(points: list[tuple[float, float, float]], p: tuple[float, float, float], radius: float):
    neighbours = []
    for q in points:
        d2 = sum((q[i] - p[i]) ** 2 for i in range(3))
        w = weight(d2, radius)
        if w > 0:
            neighbours.append((q, w))
    sw = sum(w for _, w in neighbours)
    c = [sum(w * q[i] for q, w in neighbours) / sw for i in range(3)]
    cov = [[0.0] * 3 for _ in range(3)]
    for q, w in neighbours:
        d = [q[i] - c[i] for i in range(3)]
        for i in range(3):
            for j in range(3):
                cov[i][j] += w * d[i] * d[j] / sw
    n = smallest_eigenvector(cov)
    signed = sum((p[i] - c[i]) * n[i] for i in range(3))
    return tuple(p[i] - signed * n[i] for i in range(3)), n


def main() -> None:
    points = [(float(x), float(y), 0.0) for y in range(-6, 7) for x in range(-6, 7)]
    p = (0.25, -0.5, 1.75)
    n = (0.0, 0.0, 0.0)
    for _ in range(5):
        p, n = project_once(points, p, 4.0)
    assert abs(p[2]) < 1e-10, p
    assert n[2] > 0.999999, n
    print("PASS plane projection")
    print(f"vertex={p}")
    print(f"normal={n}")


if __name__ == "__main__":
    main()
