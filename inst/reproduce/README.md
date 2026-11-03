# Reproduction of published results

Each script compares results of `onestepcox` with results published by
others and writes a table with a status per check (PASS or FAIL) to
`results/`. Run the scripts from the root of the package after
`devtools::load_all()` (or with the installed package);
`inst/validation/summarize_validation.R` summarizes all tables.

| Script | Published result | Checks | Time |
|---|---|---|---|
| `reproduce_Berry1991.R` | Berry, Kitchin and Mock (1991, Statistics in Medicine 10: 749-755), Figure 1: mean Peto and Pike estimates of the hazard ratio over 1000 simulated trials, for 50, 60, and 70 per cent censoring | 26 | a few minutes |
| `reproduce_TherneauGrambsch2000.R` | Therneau and Grambsch (2000, Modeling Survival Data, Appendix E), test data 1: closed-form score, information, and Cox estimate for the Breslow and Efron handling of ties | 24 | seconds |

## Berry, Kitchin and Mock (1991)

The article gives its results in figures only. The design of its simulation
(100 subjects per group, the group of each death drawn with probability
`n_iA theta / (n_iA theta + n_iB)`, the study stopped after 100, 80, or 60
deaths) is simulated 10,000 times. The values of Figure 1 were read from the
PDF (`Berry1991_figure1_digitized.csv`): each plotted line was followed on
both sides of its marker, and the line of the true values, which passes
through `(theta, theta)`, calibrated the axes (largest residual 0.033). The
points are those where the lines of the three series are separate: true
hazard ratios 8, 12, 15, and 16 for each censoring level and 4 for 70 per
cent; at 4 the lines coincide for 50 and 60 per cent, and at 16 only the left
side of the marker can be followed. The reading script is kept with the
project records (`Handover/checks_20261028/digitize_Berry1991_fig1.py`)
because it needs the PDF of the article. A check passes when the simulated
mean differs from the read value by at most 3 times the Monte Carlo standard
error of the difference (1000 trials in the article, 10,000 here) plus the
reading error (0.1, or 0.15 at a true hazard ratio of 16). Figure 2
(coverage of the confidence intervals) is not reproduced: its four series
overlap too much to be read.

## Therneau and Grambsch (2000)

The formulas in the script are those of the tests `book1.R` and `book2.R` of
the 'survival' package, which take them from Appendix E of the book. They
agree with the book: for the Breslow approximation, the printed estimate
1.475285 and information at the estimate 0.634168; for the Efron
approximation, the closed forms of the score and the information on
p. 327. The value printed there for the Efron information at the
estimate, 0.652077, disagrees with the closed form on the same page, which
gives 0.612632; the script checks the closed form (checked on 2026-11-03
for the first independent review of the article).

## Not yet done

`reproduce_GreenlandRobins1985.R` (the numerical example of Greenland and
Robins, 1985, Biometrics 41: 55-68, for the Mantel-Haenszel estimator and its
variance) will be written when the article is available.
