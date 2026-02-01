Auditing Free List Hamiltonian elections.
-----------------------------------------

This repository contains code for determining the set of assertions required to audit the distribution of seats in a Free List Hamiltonian election. 
The actual sampling code is copied from (an older version of) Philip Stark's [SHANGRLA audit tools](https://github.com/pbstark/SHANGRLA).

An alternative version directly uses (the current version of) the SHANGRLA code. Use
```
pip3 install git+https://github.com/pbstark/SHANGRLA.git@main
```
to install it. If you don't want to use this version, you need to comment out the relevant imports.

Data has been provided for a local election in Hesse in 2016, summarised from [Hesse online sources](https://www.statistik-hessen.de/k2016/html/index.htm).

A script `run_kk_hesse.sh` has been provided for running the experiments in the paper [Assertion-based approaches to auditing complex elections](https://arxiv.org/abs/2107.11903).
It expects directories called `output_[n]pc` for n = 5, 10, 20.

Usage (example with 0 error rate and 5% risk limit): 

```
python3 audit.py -d data/Local_Hesse_2016/Bergstrasse.csv -r 0.05 -g 0.1 -e 0 -reps 1 
```

Currently, kaplan-kolmogorov is the only supported risk function.