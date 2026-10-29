"""Expected values of tests/testthat/test-summarize_oc.R.

Independent implementation, in numpy, of the summary of one method in one
scenario of CompareOC(), applied to estimates written by hand (the same
numbers as oc_inputs() in the test file). Quantiles are numpy's default
linear interpolation, which is type 7 of R's quantile(). The hazard ratios
with two decimals are formatted with '%.2f', as R's sprintf("%.2f").

Usage: python3 tools/test_expected/summarize_oc_expected.py
"""
import numpy as np

nan = float('nan')
# Method: one replicate without an estimate (the third)
x_coef = np.log(np.array([0.50, 0.70, nan, 1.25, 0.92, 2.10]))
x_se = np.array([0.20, 0.25, nan, 0.30, 0.22, 0.35])
x_lower = np.array([0.35, 0.45, nan, 0.70, 0.60, 1.20])
x_upper = np.array([0.72, 1.10, nan, 2.20, 1.40, 3.60])
# Reference: one replicate without an estimate (the fifth)
r_coef = np.array([np.log(0.503), np.log(0.70) + 1e-4, np.log(0.80),
                   np.log(1.25) - 0.02, nan, np.log(2.104)])
r_se = np.array([0.21, 0.24, 0.26, 0.31, nan, 0.33])
hr_true = 0.8
time_per_sim = 0.0012

ok = ~np.isnan(x_coef)
pair = ok & ~np.isnan(r_coef)
dev = np.abs(x_coef[pair] - r_coef[pair])
devse = dev / r_se[pair]
ratio = x_se[pair] / r_se[pair]
mismatch = [('%.2f' % np.exp(a)) != ('%.2f' % np.exp(b))
            for a, b in zip(x_coef[pair], r_coef[pair])]
err = x_coef[ok] - np.log(hr_true)
covered = (x_lower[ok] <= hr_true) & (hr_true <= x_upper[ok])

values = [
    ('n_NA', float(np.sum(~ok))),
    ('dev_median', np.quantile(dev, 0.5)),
    ('dev_q90', np.quantile(dev, 0.9)),
    ('dev_q99', np.quantile(dev, 0.99)),
    ('dev_max', np.max(dev)),
    ('devse_median', np.quantile(devse, 0.5)),
    ('devse_q99', np.quantile(devse, 0.99)),
    ('devse_max', np.max(devse)),
    ('mismatch_2dp', float(np.mean(mismatch))),
    ('se_ratio_median', np.quantile(ratio, 0.5)),
    ('se_ratio_min', np.min(ratio)),
    ('se_ratio_max', np.max(ratio)),
    ('bias', np.mean(err)),
    ('sd', np.std(x_coef[ok], ddof=1)),
    ('rmse', np.sqrt(np.mean(err ** 2))),
    ('se_mean', np.mean(x_se[ok])),
    ('coverage', float(np.mean(covered))),
    ('time_per_sim', time_per_sim),
]
print('pairs with both estimates:', int(np.sum(pair)), '; mismatches:', mismatch)
print('Expected values (17 significant digits):')
for name, v in values:
    print(f'  {name} = {v:.17g},')
# Check of the relation between bias, sd, and rmse (divisor n - 1 for sd)
n = int(np.sum(ok))
lhs = np.mean(err ** 2)
rhs = np.mean(err) ** 2 + np.var(x_coef[ok], ddof=1) * (n - 1) / n
print(f'rmse^2 = {lhs:.17g}, bias^2 + sd^2 (n - 1) / n = {rhs:.17g}')
