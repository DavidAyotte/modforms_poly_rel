# Polynomials relations between modular forms rings generators

This package contains certain utility functions that finds polynomials relations between
the generators of the graded ring of modular forms for Gamma0(N).

Note that this implementation works only if the generators all have the same weight.

EXAMPLE:

```python
    sage: M = ModularFormsRing(Gamma0(6))
    sage: relations(M, 4)
    [ 0  0  1 -1  2 11]
    sage: g1, g2, g3 = M.gen_forms()
    sage: (g1*g3 - g2*g2 + 2*g2*g3 + 11*g3*g3).is_zero()
    True
```

More relations can be found in the file [result.txt](https://github.com/DavidAyotte/modforms_poly_rel/blob/main/result.txt).

## Content of the repo:

* [relation.sage](https://github.com/DavidAyotte/modforms_poly_rel/blob/main/relations.sage): contains all the utility functions;
* [script.sage](https://github.com/DavidAyotte/modforms_poly_rel/blob/main/script.sage): contains the script that execute the functions defined in [relation.sage](https://github.com/DavidAyotte/modforms_poly_rel/blob/main/relations.sage);
* [result.txt](https://github.com/DavidAyotte/modforms_poly_rel/blob/main/result.txt): contains the result of the script computations.

If you only need to have access to the utility function, after cloning this repo, you should write in your SageMath session:

```sage
    sage: load('relation.sage')
```

Note: this code was tested for [SageMath](https://www.sagemath.org/) version 9.0 and above.