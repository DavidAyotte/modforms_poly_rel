# Polynomials relations between modular forms rings generators

This package contains certain utility functions that finds polynomials relations between
the generators of the graded ring of modular forms for Gamma0(N).

## EXAMPLE:

```sage
    sage: load("relations.sage")
    sage: M = ModularFormsRing(Gamma0(6))
    sage: relations(M, 4)
    [    1  2/11  1/11 -1/11     0     0]
    sage: f1,f2,f3,f4,f5,f6 = homogeneous_monomials_of_weight(M, 4)
    sage: (f1 + (2/11)*f2 + (1/11)*f3 - (1/11)*f4).is_zero()
    True
```

More relations can be found in the file [result.txt](https://github.com/DavidAyotte/modforms_poly_rel/blob/main/result.txt).

## Content of the repo:

* [relations.sage](https://github.com/DavidAyotte/modforms_poly_rel/blob/main/relations.sage): contains all the utility functions;
* [script.sage](https://github.com/DavidAyotte/modforms_poly_rel/blob/main/script.sage): contains the script that execute the functions defined in [relations.sage](https://github.com/DavidAyotte/modforms_poly_rel/blob/main/relations.sage);
* [result.txt](https://github.com/DavidAyotte/modforms_poly_rel/blob/main/result.txt): contains the result of the script computations.

After cloning this repo, if you only want to use the utility functions (and don't run the script), you should execute the following command in your SageMath session:

```sage
    sage: load('relations.sage')
```

This will only load the functions defined in `relations.sage` in your current session. To see the docstring of a given function, you can write `name_of_function?`.
