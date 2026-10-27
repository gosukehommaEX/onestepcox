"""Independent reference for LimitHR(): large-sample limits of the estimators.

Written from the definitions with scipy: the analysis time solves the event
fraction computed by numerical integration over the accrual time (not the
closed form used in the package), each limiting sum is a separate quad
integral, and the score is O_T minus the integral of p dN.
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

METHODS = ["MH_Halley", "MH_Newton", "Pike_Halley", "Pike_Newton", "MH_Chebyshev",
           "Pike_Chebyshev", "MH", "Pike", "Peto", "Cox"]
QOPT = dict(epsabs=0, epsrel=1e-13, limit=500)


def limits(theta, hazard_C, n_T, n_C, a, E, dropout_hazard):
    lam = {"T": theta * hazard_C, "C": hazard_C}
    pi = {"T": n_T / (n_T + n_C), "C": n_C / (n_T + n_C)}
    lamD = dropout_hazard
    e_frac = E / (n_T + n_C)

    def F_event(j, f):
        r = lam[j] + lamD
        return lam[j] / r * (1 - np.exp(-r * f))

    def frac_events(tau):
        w = min(tau, a)
        return sum(pi[j] * quad(lambda u: F_event(j, tau - u), 0, w, **QOPT)[0] / a for j in "TC")

    hi = a
    while frac_events(hi) <= e_frac:
        hi *= 2
    tau = brentq(lambda s: frac_events(s) - e_frac, 1e-12, hi, xtol=1e-13, rtol=1e-15, maxiter=500)

    def G(t):
        return min(1.0, (tau - t) / a)

    def y(j, t):
        return pi[j] * G(t) * np.exp(-(lam[j] + lamD) * t)

    brk = [0.0, tau - a, tau] if tau > a else [0.0, tau]

    def integ(f):
        return sum(quad(f, lo, hi_, **QOPT)[0] for lo, hi_ in zip(brk[:-1], brk[1:]))

    def dN(t):
        return lam["T"] * y("T", t) + lam["C"] * y("C", t)

    def p(t, eta):
        a_ = y("T", t) * np.exp(eta)
        return a_ / (y("C", t) + a_)

    O_T = integ(lambda t: lam["T"] * y("T", t))
    O_C = integ(lambda t: lam["C"] * y("C", t))
    E_T = integ(lambda t: y("T", t) / (y("T", t) + y("C", t)) * dN(t))
    E_C = integ(lambda t: y("C", t) / (y("T", t) + y("C", t)) * dN(t))
    V = integ(lambda t: y("T", t) * y("C", t) / (y("T", t) + y("C", t)) ** 2 * dN(t))
    R = integ(lambda t: lam["T"] * y("T", t) * y("C", t) / (y("T", t) + y("C", t)))
    S = integ(lambda t: lam["C"] * y("C", t) * y("T", t) / (y("T", t) + y("C", t)))

    def derivs(eta):
        U = O_T - integ(lambda t: p(t, eta) * dN(t))
        I = integ(lambda t: p(t, eta) * (1 - p(t, eta)) * dN(t))
        J = integ(lambda t: p(t, eta) * (1 - p(t, eta)) * (1 - 2 * p(t, eta)) * dN(t))
        return U, I, J

    out = {"MH": np.log(R / S), "Pike": np.log(O_T * E_C / (O_C * E_T)), "Peto": (O_T - E_T) / V}
    for anc in ["MH", "Pike"]:
        e0 = out[anc]
        U, I, J = derivs(e0)
        out[anc + "_Newton"] = e0 + U / I
        out[anc + "_Halley"] = e0 + 2 * U * I / (2 * I * I + U * J)
        out[anc + "_Chebyshev"] = e0 + U / I - J * U * U / (2 * I ** 3)
    out["Cox"] = brentq(lambda e: derivs(e)[0], np.log(theta) - 1, np.log(theta) + 1, xtol=1e-14, maxiter=500)
    return tau, out


DESIGNS = [
    # name, hazard_C, n_T, n_C, a, E, dropout_hazard, HRs
    ("A", np.log(2) / 18, 100, 100, 12, 160, 0.0, [0.1, 0.5, 0.7, 1.0, 2.0, 10.0]),
    ("B", np.log(2) / 6, 200, 100, 24, 90, 0.02, [0.5, 3.0]),
]
