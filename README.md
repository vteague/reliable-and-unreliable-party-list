Auditing Free List Hamiltonian elections, with extensions to Sainte Laguë.
-----------------------------------------

The main branch of this repository contains code for determining the set of assertions required to audit the distribution 
of seats in a  Sainte-Laguë party-list election where voters may cast complex ballots that mix their
votes among different party lists.

This code was used to compute the estimates in the paper "Assorter based RLAs for Complex
 Sainte-Laguë elections: Reducing sample sizes by trusting manually counted ballots," by Budurushi,
Hilt, Mack, Teague and Volkamer. Currently only a fixed average fraction of untrusted ballots is supported
in the sample size estimation.

The actual sampling code uses Philip Stark's [SHANGRLA audit tools](https://github.com/pbstark/SHANGRLA). 
Run
```
pip3 install git+https://github.com/pbstark/SHANGRLA.git@main
```
to install SHANGRLA. 

Data has been provided for 
- local elections in Karlsruhe, 2019 and 2024,
- other German local elections, 2024, though not in a format currently read by the code,
- an artificial contest based on Karlsruhe data, but with a very small margin. 

Usage (Sainte Laguë example with 0 error rate and 5% risk limit):
```
python3 audit.py -d data/Karlsruhe/2024_local-council_Karlsruhe.csv -r 0.05 -g 0.1 -e 0 -reps 1 -scf sainte-lague 
```

Complete command-line options:
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

To replicate the simulation summarised in Table 1 of the paper, run
```
python3 audit.py -d data/Karlsruhe/2024_local-council_Karlsruhe.csv -r 0.05 -g 0.1 -e 0 -reps 1 -scf sainte-lague 
python3 audit.py -d data/Karlsruhe/2019_local-council_Karlsruhe.csv -r 0.05 -g 0.1 -e 0 -reps 1 -scf sainte-lague 
python3 audit.py -d data/Karlsruhe/artificial_close_contest.csv -r 0.05 -g 0.1 -e 0 -reps 1 -scf sainte-lague 
```