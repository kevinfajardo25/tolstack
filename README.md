# tolstack

1D tolerance stack-up analysis for mechanical assemblies: worst case, RSS (root sum square) and Monte Carlo, with yield, Cpk and a ranking of which part tolerances matter most.

## Why

A worst-case stack tells you whether an assembly can *ever* fail. A statistical stack tells you *how often* it will. Most real design decisions sit between the two, and the useful question is usually "which tolerance should I tighten?" This tool answers all three from a simple CSV file.

## Install and run

```bash
git clone https://github.com/kevinfajardo25/tolstack.git
cd tolstack
pip install -e ".[test]"
python -m tolstack examples/hinge_gap.csv --lsl 0.10 --usl 0.70
pytest -q
```

## Defining a stack

Each row is one dimension in the loop. `direction` is `1` if the dimension opens the gap and `-1` if it closes it. Tolerances are magnitudes, so `+0.15 / -0.05` is written `0.15, 0.05`. `dist` is `normal` (tolerance equals plus or minus 3 sigma) or `uniform`.

```csv
name,nominal,plus,minus,direction,dist
housing_cavity_width,40.00,0.10,0.10,1,normal
left_wall_spacer,2.00,0.05,0.05,-1,normal
pcb_width,30.00,0.10,0.10,-1,normal
connector_overhang,1.50,0.15,0.05,-1,normal
right_wall_spacer,2.00,0.05,0.05,-1,normal
foam_pad_compressed,4.10,0.20,0.20,-1,uniform
```

## Worked example: PCB in a housing

The clearance gap must stay between 0.10 and 0.70 mm.

```
6 dimensions, nominal gap 0.4000, mean gap 0.3500

worst_case   mean    0.3500   range [  -0.2500,    0.9500]
rss          mean    0.3500   range [  -0.0437,    0.7437]   yield 96.778%   Cpk  0.64
monte_carlo  mean    0.3498   range [   0.0234,    0.6745]   yield 98.077%   Cpk  0.63

Variance contribution:
  foam_pad_compressed       77.4%  #######################################
  housing_cavity_width       6.5%  ###
  pcb_width                  6.5%  ###
  connector_overhang         6.5%  ###
  left_wall_spacer           1.6%  #
  right_wall_spacer          1.6%  #
```

Reading the result:

* The worst case goes negative, so parts can interfere.
* About 2 to 3% of assemblies fall outside spec, and Cpk is 0.64, far below the usual 1.33 target.
* The asymmetric connector tolerance shifts the mean gap from 0.40 to 0.35 mm.
* The foam pad drives 77% of the variation. Tightening every machined part would barely help.

Changing only the foam pad to a 4.05 +/- 0.10 mm controlled part (`examples/hinge_gap_improved.csv`):

```
rss          mean    0.4000   range [   0.1879,    0.6121]   yield 99.998%   Cpk  1.41
monte_carlo  mean    0.3999   range [   0.1863,    0.6110]   yield 99.996%   Cpk  1.41
```

One part change takes Cpk from 0.64 to 1.41.

## Method notes

* **Asymmetric tolerances** become a shifted mean plus a symmetric half tolerance: `mean = nominal + (plus - minus)/2`, `half = (plus + minus)/2`.
* **RSS** assumes independent parts, with sigma equal to half tolerance divided by 3 for normal parts and divided by the square root of 3 for uniform parts.
* **Monte Carlo** samples each part from its distribution. The range reported is the 0.135 to 99.865 percentile band, the same coverage as plus or minus 3 sigma, so it can be compared directly with RSS.
* **Cpk** is the distance from the mean to the nearest spec limit divided by 3 sigma.
* The tool is 1D. Geometric tolerances (position, flatness, perpendicularity) must be converted into equivalent linear contributions before they are entered.

## Tests

13 tests cover the analytic cases (three equal tolerances give an RSS of tol times the square root of 3), asymmetric tolerances, uniform sigma, agreement between Monte Carlo and RSS for yield and Cpk, and input validation.

## License

MIT
