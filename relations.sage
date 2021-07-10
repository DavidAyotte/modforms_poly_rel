
r"""
This package contains certain utility functions that finds polynomials relations between
the generators of the graded ring of modular forms for Gamma0(N).

Note that this implementation works only if the generators all have the same weight.

EXAMPLES::

    sage: M = ModularFormsRing(Gamma0(6))
    sage: relations(M, 4)
    [ 0  0  1 -1  2 11]
    sage: g1, g2, g3 = M.gen_forms()
    sage: (g1*g3 - g2*g2 + 2*g2*g3 + 11*g3*g3).is_zero()
    True

AUTHORS:

- David Ayotte (2021-07-05): initial version

"""

# ****************************************************************************
#       Copyright (C) 2021 David Ayotte
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from itertools import combinations_with_replacement

def nb_of_generators(M):
    r"""
    Return the number of generator for given ring of modular forms.

    EXAMPLES::

        sage: nb_of_generators(ModularFormsRing(1))
        2
        sage: nb_of_generators(ModularFormsRing(Gamma0(6)))
        3
    """
    return len(M.gen_forms())

def weights_of_generators(M):
    r"""
    Return a list of the weights of generators of the given `ModularFormsRing`

    EXAMPLES::

        sage: M = ModularFormsRing(1)
        sage: weights_of_generators(M)
        [4, 6]
        sage: M = ModularFormsRing(Gamma0(6))
        [2, 2, 2]
    """
    return [ M.generators()[f][0] for f in range(0, nb_of_generators(M)) ]

def check_generators_weights(M):
    r"""
    Check whether the generators of `M` all have the same weights or not.

    EXAMPLE::

        sage: M = ModularFormsRing(1)
        sage: check_generators_weights(M)
        False
        sage: M = ModularFormsRing(Gamma0(6))
        sage: check_generators_weights(M)
        True
    """
    if len(set(weights_of_generators(M))) == 1:
        return True
    else:
        return False

def homogeneous_monomials_of_weight(M, weight):
    r"""
    Returns the list of all homogeneous monomials of weight `weight`.

    EXAMPLES::

        sage: M = ModularFormsRing(Gamma0(6))
        sage: M.gen_forms()
        [1 + 24*q^3 + O(q^6),
        q + 5*q^3 - 2*q^4 + 6*q^5 + O(q^6),
        q^2 - 2*q^3 + 3*q^4 + O(q^6)]
        sage: g0, g1, g2 = M.gen_forms()
        sage: homogeneous_monomials_of_weight(M, 4)
        [q^4 - 4*q^5 + O(q^6),
        q^3 - 2*q^4 + 8*q^5 + O(q^6),
        q^2 - 2*q^3 + 3*q^4 + 24*q^5 + O(q^6),
        q^2 + 10*q^4 - 4*q^5 + O(q^6),
        q + 5*q^3 + 22*q^4 + 6*q^5 + O(q^6),
        1 + 48*q^3 + O(q^6)]
        sage: g0*g0
        1 + 48*q^3 + O(q^6)
        sage: g0*g1
        q + 5*q^3 + 22*q^4 + 6*q^5 + O(q^6)
        sage: g0*g2
        q^2 - 2*q^3 + 3*q^4 + 24*q^5 + O(q^6)
        sage: g1*g1
        q^2 + 10*q^4 - 4*q^5 + O(q^6)
        sage: g1*g2
        q^3 - 2*q^4 + 8*q^5 + O(q^6)
        sage: g2*g2
        q^4 - 4*q^5 + O(q^6)
    """
    G = M.group()
    N = G.level()
    if G is not SL2Z and G == Gamma1(N):
        raise NotImplementedError
    gen_list = M.gen_forms()
    weights_of_generators = [g.weight() for g in gen_list]
    comb = WeightedIntegerVectors(weight, weights_of_generators).list()
    if len(comb) == 0:
        raise ValueError("there is no modular forms of the given weight (%s) for %s"%(weight, M.group()))

    monomials = []
    for exponents in comb:
        monomial_forms = 1
        iter = 0
        for e, g in zip(exponents, gen_list):
            pow_g = prod(g for i in range(0, e))
            monomial_forms *= pow_g
        monomials.append(monomial_forms)
    return monomials

def initialize_matrix_of_coefficients(monomials):
    r"""
    Given a list of monomials, this function returns a matrix of the coefficient of the given monomials.

    EXAMPLES::

    sage: M = ModularFormsRing(Gamma0(6))
        sage: mon = homogeneous_monomials_of_weight(M, 4); mon
        [q^4 - 4*q^5 + O(q^6),
        q^3 - 2*q^4 + 8*q^5 + O(q^6),
        q^2 - 2*q^3 + 3*q^4 + 24*q^5 + O(q^6),
        q^2 + 10*q^4 - 4*q^5 + O(q^6),
        q + 5*q^3 + 22*q^4 + 6*q^5 + O(q^6),
        1 + 48*q^3 + O(q^6)]
        sage: initialize_matrix_of_coefficients(mon)
        [ 0  0  0  0  0  1]
        [ 0  0  0  0  1  0]
        [ 0  0  1  1  0  0]
        [ 0  1 -2  0  5 48]
        [ 1 -2  3 10 22  0]
        [-4  8 24 -4  6  0]
    """
    if len(monomials) == 0:
        raise ValueError
    parent = monomials[0].parent()
    matrix_data = []
    for f in monomials:
        matrix_data.append(f.coefficients([0..parent.sturm_bound()]))
    return Matrix(matrix_data).transpose()

def relations(M, weight):
    r"""
    Given a ring of modular forms `M`, returns the algebraic relations in weight `weight` in terms
    of a matrix. More precisely, the coefficients of each rows of the matrix gives the coefficients
    of the relations between of the monomials of weight `weight`.

    ..NOTE::

        This function currently works only if the weights the generators are all equal.

    EXAMPLES::

        sage: M = ModularFormsRing(Gamma0(6))
        sage: relations(M, 4)
        [    1  2/11  1/11 -1/11     0     0]
        sage: f1,f2,f3,f4,f5,f6 = homogeneous_monomials_of_weight(M, 4)
        sage: (f1 + (2/11)*f2 + (1/11)*f3 - (1/11)*f4).is_zero()
        True
    """
    monomials = homogeneous_monomials_of_weight(M, weight)
    A = initialize_matrix_of_coefficients(monomials)
    ker = A.right_kernel()
    return ker.matrix()

def check_generators_relations(M, relation_matrix, weight):
    r"""
    Verificator function that check if the relations are indeed true.
    """
    monomials = homogeneous_monomials_of_weight(M, weight)
    for row in range(0, relation_matrix.nrows()):
        rel = 0
        for col in range(0, relation_matrix.ncols()):
            rel += relation_matrix[row, col]*monomials[col]
        if not rel.is_zero():
            return False
    return True