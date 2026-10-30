# Scripts of the article

The figures, tables, and numbers of the article are computed by the scripts
of this folder (main text) and of `inst/sup_info/` (Supporting
Information). Each part has two steps: a data script, which runs the
simulations and saves the results in `results/data/`, and a table and figure
script, which reads the saved results and writes `results/tables/` and
`results/figures/`. The settings shared by all scripts (design, methods,
scenarios, seeds, the reference, and the quick test run) are in
`settings.R`, and the helpers of the two table and figure scripts in
`figure_helpers.R`.

## How to run

Open `onestepcox.Rproj`, so that the working directory is the root of the
package, and run

```r
devtools::load_all()
# Quick test run: 1/100 of the replicates, results in results/quick/
options(onestepcox.quick = TRUE)
source("inst/scripts/data_generate_main.R")
source("inst/sup_info/data_generate_supplement.R")
source("inst/scripts/table_and_figure_manuscript.R")
source("inst/sup_info/table_and_figure_supplement.R")
# Full run
options(onestepcox.quick = FALSE)
source("inst/scripts/data_generate_main.R")
source("inst/sup_info/data_generate_supplement.R")
source("inst/scripts/table_and_figure_manuscript.R")
source("inst/sup_info/table_and_figure_supplement.R")
```

The packages 'asaur' (pharmacoSmoking data) and 'microbenchmark' (computing
times) are needed besides the imports of the package; the data scripts check
them first. The data scripts write `run_log_main.csv` and
`run_log_supplement.csv` (elapsed minutes of each section and the versions
of R and the packages) and the session information to `results/data/`. The
largest scenario of grid A (900 subjects, 10,000 replicates) generates 9
million rows at once and needs about 2 GB of memory.

## Design

All simulations use `GenSimData()`: exponential survival times with control
hazard log(2)/18 per month, uniform accrual over 12 months, and the analysis
at the E-th event with E = 0.4 (n_T + n_C), rounded. The methods are `Cox`
(`survival::coxph()` with the Efron handling of ties and `eps = 1e-12`),
`MH_Halley`, `MH_Newton`, `Pike_Halley`, `Pike_Newton`, `MH`, `Pike`,
`Peto`, and `Pike_Chebyshev`, and the reference `Cox_NR`. No estimate is
trimmed.

Reference of the deviations: `Cox_NR`, the root of the Cox score found by
Newton-Raphson iteration on the risk-set summaries, which stops when the
step is below 1e-12 and is accurate to rounding. `survival::coxph()` stops
on the relative change of the log partial likelihood, so that its estimate
can differ from the root by about 1e-12 to 1e-11, which would hide the
accuracy of the proposed estimator at large n; its difference from `Cox_NR`
is reported in the row of `Cox`.

Replicates in which the reference has no estimate are left out of the
summaries of all methods, so that every method is summarized over the same
replicates; their number is `n_undefined`. `Cox_NR` has no estimate when the
Cox estimate is not finite: when a group has no event (which happens in
small samples with extreme hazard ratios), or when every event of a group
occurs at a time when no subject of the other group is at risk in the
stratum, so that the Mantel-Haenszel anchor is not finite. It has none
either, with a warning, when the Newton-Raphson iteration does not converge.
When the Cox estimate is not finite, `survival::coxph()` warns that the
coefficient may be infinite. These warnings are counted (`n_warning`) and not
printed. Because of its tolerance (`toler.inf`, the square root of
`eps` = 1e-12, so 1e-6), `survival::coxph()` can give the same warning in
replicates with a finite Cox estimate, which stay in the summaries, so that
`n_warning` can exceed `n_undefined`.

| Section | Script | Content | Replicates | Base seed |
|---|---|---|---|---|
| `limit` | main | `LimitHR()` at 81 hazard ratios from 0.1 to 10 (equally spaced on the log scale) and at the hazard ratios of `fixed`, 1:1 | none | none |
| `fixed` | main | 1:1, HR 0.1, 0.2, 0.5, 2, 5, 10, n_T + n_C = 100, 300, 1,000, 3,000, 10,000, 30,000, 100,000 (42 scenarios; n varies fastest) | 1,000 for n up to 3,000, 200 for 10,000 and 30,000, 100 for 100,000 | 100001 |
| `local` | main | HR exp(c / sqrt(n)) with c = -10 log 2 and 10 log 2 (HR 0.5 and 2 at n = 100), the same n (14 scenarios) | as `fixed` | 200001 |
| `gridA` | main | Continuous times; n_C = 30, 100, 300; n_T = r n_C with r = 1, 2; HR 0.2, 0.5, 0.7, 1, 2, 5 (36 scenarios; HR fastest, then r, then n_C) | 10,000 | 300001 |
| `gridB` | main | n_T = n_C = 100, HR 0.2, 0.5, 1, 2, 5, with B1 dropout of 10 per cent per year, B2 times rounded up to months (ties), B3 2 strata (probabilities 0.5 and 0.5, hazard multipliers 1 and 3), B4 20 strata (equal probabilities, multipliers exp(-1) to exp(1)), B5 ties and 20 strata (25 scenarios) | 10,000 | 400001 |
| `application` | main | gbsg of 'survival' (`hormon`, control 0) and pharmacoSmoking of 'asaur' (`grp`, control `"patchOnly"`), without and with strata (`grade`, `ageGroup2`), Efron and Breslow | none | none |
| `limit_sup` | sup_info | `limit` with 2:1 allocation and with dropout of 10 per cent per year | none | none |
| `gridB_breslow` | sup_info | B2 and B5 with the Breslow handling of ties for all methods | 10,000 | 450001 |
| `gridC` | sup_info | n_T + n_C = 30, 60; r = 1, 2; HR 0.1, 0.2, 0.5, 2, 5, 10 (24 scenarios) | 10,000 | 500001 |
| `timing` | sup_info | One dataset each with n = 200, 1,000, 5,000, 20,000, HR 0.5: aggregation of the risk sets, the Mantel-Haenszel anchor and Halley step, Newton-Raphson iteration from the same anchor (tolerance 1e-8 and 1e-12), `EstHR()` with `"MH_Halley"`, `"Cox_NR"`, and `"Cox"`, and `survival::coxph()` | 200 timings each (`microbenchmark`) | 600001 |

Seeds: scenario i of a section uses base + 100 (i - 1), as `CompareOC()`
does; for n of 30,000 or more the replicates are run in chunks of 10, and
chunk k uses base + 100 (i - 1) + (k - 1). The bases are at least 50,000
apart and a section uses at most 4,200 seeds above its base, so that no
seed is used twice. The estimates of every replicate are saved for `fixed`
and `local`; for the grids only the summaries are saved, because the
estimates of `gridA` alone would take about 100 MB.

## Results

Main text (`inst/scripts/results/`):

| File | Content |
|---|---|
| `data/limit.rds` | `LimitHR()` results: `grid` and `fixed` |
| `data/fixed.rds`, `local.rds` | Lists with `scenarios`, `estimates` (every replicate and method, unchanged), and `summary` (the columns of the summary of `CompareOC()` against `Cox_NR`, with `nsim`, `n_undefined`, and `n_warning`) |
| `data/gridA.rds`, `gridB.rds` | Lists with `scenarios` and `summary` (as above, with `HR_true`; attribute `nsim`) |
| `data/application.rds`, `application.csv` | Estimates of the applications (the CSV with 17 significant digits) |
| `data/fig4_envelope.rds` | Constants of the dashed lines of Figure 4 |
| `tables/table1_operating.tex` | Table 1: coverage and bias, n_T = n_C = 100 |
| `tables/table2_application.tex` | Table 2: applications, Efron |
| `figures/Fig1.pdf`, `.eps` | Large-sample bias of the Pike-anchored and closed-form estimators against the hazard ratio; absolute biases below 1e-10 (near HR 1, at the accuracy of the numerical integration) are not drawn |
| `figures/Fig2.pdf`, `.eps` | Deviation from the Cox estimate against n (median and 99th percentile; dotted lines: large-sample limits of the Pike-anchored estimators) |
| `figures/Fig3.pdf`, `.eps` | The deviation divided by the standard error of the Cox estimate |
| `figures/Fig4.pdf`, `.eps` | Error after one Newton or Halley step against the error of the anchor; in each scenario 200 random replicates and the 20 with the largest anchor errors are drawn; dashed lines: K e^2 and K e^3 with K the largest ratio over all replicates with an anchor error of at least 1e-3 (below it the error after one step is near the rounding error) |

Supporting Information (`inst/sup_info/results/`):

| File | Content |
|---|---|
| `data/limit_sup.rds`, `gridB_breslow.rds`, `gridC.rds`, `timing.rds` | Results of the sections of the supplement |
| `tables/tableS1_scenarios.tex` | Table S1: the scenarios |
| `tables/tableS2_gridA.tex`, `tableS3_gridB.tex` | Tables S2 and S3: all summaries of grids A and B (B2 and B5 also with Breslow ties); long tables |
| `tables/tableS4_gridC.tex` | Table S4: small samples and extreme hazard ratios; long table |
| `tables/tableS5_timing.tex` | Table S5: computing time |
| `tables/tableS6_application.tex` | Table S6: applications, all methods, Efron and Breslow; long table |
| `tables/tableS7_validation.tex` | Table S7: results of `inst/reproduce/` and `inst/validation/` (from `inst/validation/results/summary.csv`) |
| `figures/FigS1.pdf`, `.eps` | Figure S1: Figure 1 with 2:1 allocation and with dropout |
| `figures/FigS2.pdf`, `.eps` | Figure S2: Figures 2 and 3 for HR 0.2 and 5 |

The tables need the LaTeX packages booktabs and array (the ragged-right
paragraph columns of Tables S1 and S7). The long tables (Tables S2 to S4 and
S6) need longtable, and Tables S2 to S4, which are set on landscape pages,
need pdflscape.

## Numbers of the article and copies to the manuscript

Every number written in the main text, the Supporting Information, or the
cover letter is entered in the ledger `number_locations.csv` before it is
written. The ledger is a UTF-8 file, edited as text (Excel's plain CSV is in
another encoding). Its columns are `id`, `file` (`main.tex`,
`supporting_information.tex`, or `cover_letter.tex`), `anchor` (text of the
.tex around the number, found exactly once; avoid `\"` in it), `text_value`
(the number as written), `source_file` (an .rds or .csv file of the results,
relative to the root of the package), `source_expr` (one R expression),
`rounding` (the rule that writes the value as in the text), and `type`:

- `computed`: a result; `source_expr` uses `x`, the content of
  `source_file`;
- `input`: a setting of the simulations; `source_expr` uses an object of
  `settings.R` (such as `n_rep(10000L)` or `event_frac`); a design value that
  is not in `settings.R` is read from the results as `computed`;
- `literature:<bibkey>`: a number taken from a cited work; only the
  citation of `<bibkey>` in the same file is checked by the scripts, and the
  value is checked against the cited work by hand.

The rounding rules are described at the top of `number_check.R` and tested
on `number_rules_test.csv`. The values are rounded with `round()`, as the
tables are (`fmt_fixed()` in `figure_helpers.R`), or with `signif()` for the
rules `signif` and `approx`, so that at a decimal tie (such as 94.45 written
with one decimal) the text follows R. The rules `signif` and `approx` write
no trailing zeros (0.10 needs `fixed:2`). A negative number is written
`$-0.012$`, and thousands are separated by `,` or `{,}` (not `\,`).

After any change of the results, of `settings.R`, of the ledger, or of the
text, run, from the root of the package in a fresh R session (in RStudio,
Session > Restart R),

```r
devtools::load_all()
source("inst/scripts/number_check.R")
source("inst/scripts/copy_results.R")
```

and then, from the folder `Manuscript`,

```
python3 tools/verify_copies.py
python3 tools/check_numbers.py
python3 tools/classify_numbers.py
```

Each Python script passes with `RESULT: PASS` (exit status 0); the reports
are in `Manuscript/tools/logs/`.

- `number_check.R` recomputes each number of type `computed` or `input` and
  writes `results/number_check.csv`, with the MD5 of each source file (of
  `settings.R` for an input without a source file).
- `copy_results.R` copies the `.tex` tables and the `.pdf` and `.eps`
  figures of `results/` (not those of `results/quick/`) to
  `Manuscript/Main_body/` and `Manuscript/Supp_info/`, and records their MD5
  in `Manuscript/tools/logs/copy_results.csv`. The manuscript folder is
  `../../Manuscript` from the root of the package, or the one given by
  `options(onestepcox.manuscript = "<path>")`. Copies of earlier runs are not
  removed; `verify_copies.py` reports them, and they are moved by hand to
  `_to_delete/`.
- `verify_copies.py` checks the copies, and that every file that `main.tex`
  and `supporting_information.tex` read with `\input`, `\include`,
  `\subfile`, or `\includegraphics` is one of them (so the text is written in
  `main.tex` and `supporting_information.tex`, not in section files).
- `check_numbers.py` checks each row of the ledger: the anchor in the text,
  the text value in it with its sign, and the result of `number_check.R`,
  which must be up to date with the ledger and the source files.
- `classify_numbers.py` lists every number of the text, outside the
  arguments of commands such as `\cite`, `\ref`, `\label`, `\input`,
  `\includegraphics`, and `\url` and outside lengths, and fails on any
  number that is not inside an anchor of the ledger, inside an anchor of
  `Manuscript/tools/number_whitelist.csv`, the number of a table, figure,
  section, equation, or the like (`Table~1`, `Figures 2--4`, `Section 2.1`;
  write `Tables 1--3` rather than `Tables 1, 2 and 3`), or a single digit in
  mathematics. A single digit in inline mathematics is a value, and fails,
  when it is all of the formula (`$2$`) or follows `=`, `<`, `>`, `\le`,
  `\neq`, or the like outside a subscript or superscript (`$\theta = 5$`,
  but not `$\sum_{i=1}^{n}$`); any single digit of a display equation is a
  formula, while a number of two or more digits is checked as in the text.
  A macro whose definition contains a number fails too, apart from the
  macros of the layout (`\arraystretch` and the like) and the definitions
  inside an anchor of the whitelist. Not read: the
  preamble apart from macro definitions (`\title`, `\author`), a number
  without a leading zero (`.5`), and a number right after `\\`.

## Computing time

In the full run of 2026-10-29 (R 4.6.1 on Windows; `run_log_main.csv` and
`run_log_supplement.csv`), `data_generate_main.R` took 90 minutes (the grids
67 minutes, `fixed` and `local` 23 minutes) and `data_generate_supplement.R`
32 minutes; the two table and figure scripts took under a minute each
(times of their output files). In the quick test run, the two data scripts
took 6 minutes.
