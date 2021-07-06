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

Note: this code was tested for SageMath version 9.0 and above.