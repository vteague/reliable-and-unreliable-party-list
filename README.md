Auditing Free List Hamiltonian elections.
-----------------------------------------

This repository contains code for determining the set of assertions required to audit the distribution of seats in a Free List Hamiltonian election. 
The actual sampling code is copied from (an older version of) Philip Stark's [SHANGRLA audit tools](https://github.com/pbstark/SHANGRLA).

Data has been provided for a local election in Hesse in 2016, summarised from [Hesse online sources](https://www.statistik-hessen.de/k2016/html/index.htm).

A script `run_kk_hesse.sh` has been provided for running the experiments in the paper [Assertion-based approaches to auditing complex elections](https://arxiv.org/abs/2107.11903).
It expects directories called `output_[n]pc` for n = 5, 10, 20.

Usage (example with 0 error rate and 5% risk limit): 

```
python3 audit.py -d data/Local_Hesse_2016/Bergstrasse.csv -r 0.05 -g 0.1 -e 0 -reps 1 
```

Command-line options:
```
-d data file
-e error rate (default 0.002)
-r risk limit (default 0.05)
-g g parameter for risk calculation in some risk functions (default 0.1)
-t threshold for testing (default 1/2)
-rf risk function (default kaplan_kolmogorov; only that one is implemented)
-reps number of repeats of random test (default 20)
-s seed (default 9368663)
-scf social choice function (default free-list-hamiltonian; sainte-lague also accepted)
```


Currently, kaplan-kolmogorov is the only supported risk function.