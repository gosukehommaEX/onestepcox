"""Independent reference implementation for the onestepcox tests.

Written from the definitions, not from the R code: risk sets are counted by
brute force (one comparison per subject and event time), the Cox estimate is
the root of the score found by Brent's method, and the derivatives are written
per event time without the weight-table form used in the package.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm


def risk_table(time, event, trt, strata=None):
    """Rows (stratum, t, nT, nC, OT, O) for each distinct event time per stratum."""
    time = np.asarray(time, float); event = np.asarray(event, int); trt = np.asarray(trt, bool)
    if strata is None:
        strata = np.zeros(len(time), int)
    strata = np.asarray(strata)
    rows = []
    for s in sorted(set(strata.tolist()), key=str):
        m = strata == s
        t_s, e_s, g_s = time[m], event[m], trt[m]
        for t in sorted(set(t_s[e_s == 1].tolist())):
            at = t_s >= t
            dead = (t_s == t) & (e_s == 1)
            rows.append((s, t, int(np.sum(at & g_s)), int(np.sum(at & ~g_s)),
                         int(np.sum(dead & g_s)), int(np.sum(dead))))
    return rows


def loglik(eta, rows, ties):
    ll = 0.0
    th = np.exp(eta)
    for (_, _, nT, nC, OT, O) in rows:
        OC = O - OT
        ll += OT * eta
        if ties == "breslow":
            ll -= O * np.log(nC + nT * th)
        else:
            for r in range(O):
                ll -= np.log((nC - r * OC / O) + (nT - r * OT / O) * th)
    return ll


def derivs(eta, rows, ties):
    """Score U and I = -U', J = -U'', L = -U''' of the log partial likelihood."""
    th = np.exp(eta)
    U = I = J = L = 0.0
    for (_, _, nT, nC, OT, O) in rows:
        OC = O - OT
        U += OT
        subsets = [(nT, nC, O)] if ties == "breslow" else \
            [(nT - r * OT / O, nC - r * OC / O, 1) for r in range(O)]
        for (a, b, w) in subsets:
            p = a * th / (b + a * th)
            q = 1 - p
            U -= w * p
            I += w * p * q
            J += w * p * q * (q - p)
            L += w * p * q * (1 - 6 * p * q)
    return U, I, J, L


def cox(rows, ties):
    f = lambda e: derivs(e, rows, ties)[0]
    eta = brentq(f, -20, 20, xtol=1e-15, rtol=4 * np.finfo(float).eps, maxiter=500)
    I = derivs(eta, rows, ties)[1]
    return eta, 1 / np.sqrt(I)


def anchors(rows):
    OT = sum(r[4] for r in rows); O = sum(r[5] for r in rows); OC = O - OT
    ET = sum(r[2] * r[5] / (r[2] + r[3]) for r in rows)
    EC = sum(r[3] * r[5] / (r[2] + r[3]) for r in rows)
    V = 0.0
    for (_, _, nT, nC, oT, o) in rows:
        n = nT + nC
        if n > 1:
            V += nT * nC * o * (n - o) / (n * n * (n - 1))
    R = sum(r[4] * r[3] / (r[2] + r[3]) for r in rows)
    S = sum((r[5] - r[4]) * r[2] / (r[2] + r[3]) for r in rows)
    vMH = sum(r[5] * r[2] * r[3] / (r[2] + r[3]) ** 2 for r in rows) / (R * S)
    return {"MH": (np.log(R / S), np.sqrt(vMH)),
            "Pike": (np.log(OT * EC / (OC * ET)), np.sqrt(1 / ET + 1 / EC)),
            "Peto": ((OT - ET) / V, 1 / np.sqrt(V))}


def one_step(eta0, rows, ties, kind):
    U, I, J, L = derivs(eta0, rows, ties)
    if kind == "Newton":
        d = U / I
    elif kind == "Halley":
        d = 2 * U * I / (2 * I * I + U * J)
    else:
        d = U / I - J * U * U / (2 * I ** 3)
    info = I + J * d + L * d * d / 2
    return eta0 + d, 1 / np.sqrt(info)


METHODS = ["MH_Halley", "MH_Newton", "Pike_Halley", "Pike_Newton", "MH_Chebyshev",
           "Pike_Chebyshev", "MH", "Pike", "Peto", "Cox_NR"]


def all_methods(time, event, trt, strata=None, ties="efron"):
    rows = risk_table(time, event, trt, strata)
    a = anchors(rows)
    out = {}
    for m in METHODS:
        if m in a:
            out[m] = a[m]
        elif m == "Cox_NR":
            out[m] = cox(rows, ties)
        else:
            anc, kind = m.split("_")
            out[m] = one_step(a[anc][0], rows, ties, kind)
    return out
