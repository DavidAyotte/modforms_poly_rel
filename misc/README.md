# Counting weights of generators

This script count the number of occurrences of the weights of the generators for a modular forms ring for `Gamma0` or `Gamma1`.

To use this script, load the file `weight_occurrences.sage` in your SageMath session.

```sage
    sage: load("weights_occurrences.sage")
```

A text file will be created in the directory with the ouput of the results. The files `Gamma0_max_30.txt` and `Gamma1_max_20.txt` already contains some results.

You may also open the file `weight_occurrences.sage` and change the initial parameters (ex: the group and the maximum level).

## Explanations

The script runs through all rings of modular forms for `Gamma0` or `Gamma1` up to a chosen maximum level and use SageMath to compute a generating set for each rings.

Then, the script will output strings of the form:
```
Gamma0(N):
    [(k1, N1), (k2, N2)]
    [[F1, k1],
     [G1, k1],
     [H1, k1],
     [F2, k2],
     [G2, k2],
     [H3, k2]]
```
where:

* `[(k1, N1), (k2, N2)]` represents the list of unique weights `ki` and their occurences `Ni`;
* `[[F1, K1],...,[H3, k2]]` is the list of all generators with their weight.
