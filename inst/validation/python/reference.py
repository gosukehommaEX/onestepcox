"""Independent Python implementation of the estimators of onestepcox::EstHR().

Written from the definitions in the documentation of EstHR(), without the
weight-table form of the package: the risk sets are counted by comparing
every subject with every event time, the Efron terms are written per tie
block, and the Cox estimate is the root of the score found by Brent's
method on a bracket (not by Newton-Raphson). The rules for missing values
follow the documentation: no estimate when there is no event in one of the
groups, when the anchor is not finite, or (Halley) when the denominator is
not positive.
"""
import numpy as np
from scipy.optimize import brentq

METHODS = ["MH_Halley", "MH_Newton", "Pike_Halley", "Pike_Newton", "MH_Chebyshev",
           "Pike_Chebyshev", "MH", "Pike", "Peto", "Cox_NR"]
NA = (np.nan, np.nan)


def risk_table(time, event, trt, stratum):
    """Arrays n_T, n_C, O_T, O over the distinct event times within strata."""
    nT, nC, OT, O = [], [], [], []
    for s in np.unique(stratum):
        m = stratum == s
        t, e, g = time[m], event[m], trt[m]
        ev_times = np.unique(t[e == 1])
        if len(ev_times) == 0:
            continue
        at_risk = t[:, None] >= ev_times[None, :]            # subjects x times
        dies = (t[:, None] == ev_times[None, :]) & (e[:, None] == 1)
        nT.append((at_risk & g[:, None]).sum(0))
        nC.append((at_risk & ~g[:, None]).sum(0))
        OT.append((dies & g[:, None]).sum(0))
        O.append(dies.sum(0))
    cat = lambda v: np.concatenate(v).astype(float) if v else np.zeros(0)
    return cat(nT), cat(nC), cat(OT), cat(O)


def terms(nT, nC, OT, O, ties):
    """Numbers at risk and weights of the terms of the log partial likelihood."""
    if ties == "breslow":
        return nT, nC, O
    a, b = [], []
    for k in range(len(O)):
        Ok = int(O[k])
        r = np.arange(Ok) / Ok
        a.append(nT[k] - r * OT[k])
        b.append(nC[k] - r * (O[k] - OT[k]))
    return np.concatenate(a), np.concatenate(b), np.ones(int(O.sum()))


def derivs(eta, tt, OT_total):
    a, b, w = tt
    p = a * np.exp(eta) / (b + a * np.exp(eta))
    q = 1 - p
    U = OT_total - np.sum(w * p)
    I = np.sum(w * p * q)
    J = np.sum(w * p * q * (q - p))
    L = np.sum(w * p * q * (1 - 6 * p * q))
    return U, I, J, L


def estimates(time, event, trt, stratum, ties):
    """Dictionary method -> (coef, se) with NaN for no estimate."""
    out = {m: NA for m in METHODS}
    if len(time) == 0 or event.sum() == 0:
        return out
    nT, nC, OT, O = risk_table(time, event, trt, stratum)
    OT_tot, O_tot = OT.sum(), O.sum()
    if OT_tot == 0 or OT_tot == O_tot:
        return out
    n = nT + nC
    OC = O - OT
    with np.errstate(divide="ignore", invalid="ignore"):
        R = np.sum(OT * nC / n)
        S = np.sum(OC * nT / n)
        mh = (np.log(R / S), np.sqrt(np.sum(O * nT * nC / n ** 2) / (R * S)))
        ET = np.sum(nT * O / n)
        EC = np.sum(nC * O / n)
        pike = (np.log(OT_tot * EC / ((O_tot - OT_tot) * ET)), np.sqrt(1 / ET + 1 / EC))
        ok = n > 1
        V = np.sum(nT[ok] * nC[ok] * O[ok] * (n[ok] - O[ok]) / (n[ok] ** 2 * (n[ok] - 1)))
        peto = ((OT_tot - ET) / V, 1 / np.sqrt(V))
    for name, val in (("MH", mh), ("Pike", pike), ("Peto", peto)):
        if np.isfinite(val[0]):
            out[name] = (val[0], val[1] if np.isfinite(val[1]) else np.nan)
    tt = terms(nT, nC, OT, O, ties)
    for anchor, val in (("MH", mh), ("Pike", pike)):
        eta0 = val[0]
        if not np.isfinite(eta0):
            continue
        U, I, J, L = derivs(eta0, tt, OT_tot)
        if not I > 0:
            continue
        steps = {"Newton": U / I,
                 "Halley": 2 * U * I / (2 * I * I + U * J) if 2 * I * I + U * J > 0 else np.nan,
                 "Chebyshev": U / I - J * U * U / (2 * I ** 3)}
        for kind, d in steps.items():
            if not np.isfinite(d):
                continue
            info = I + J * d + L * d * d / 2
            out[anchor + "_" + kind] = (eta0 + d, 1 / np.sqrt(info) if info > 0 else np.nan)
    # Cox: root of the score on a bracket; it exists when the MH estimate is finite
    if np.isfinite(mh[0]):
        f = lambda e: derivs(e, tt, OT_tot)[0]
        lo, hi = mh[0] - 1, mh[0] + 1
        while f(lo) <= 0:
            lo -= 1
        while f(hi) >= 0:
            hi += 1
        root = brentq(f, lo, hi, xtol=1e-14, rtol=1e-15, maxiter=500)
        out["Cox_NR"] = (root, 1 / np.sqrt(derivs(root, tt, OT_tot)[1]))
    return out
