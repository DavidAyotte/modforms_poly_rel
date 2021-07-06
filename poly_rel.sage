print("This SageMath script depends on the trac ticket trac.sagemath.org/ticket/31559")
from itertools import combinations_with_replacement

def weights_of_generators(M):
    return [ M.generators()[f][0] for f in range(0, M.ngens()) ]

def check_generators_weights(M):
    if len(set(weights_of_generators(M))) == 1:
        return True
    else:
        return False

def homogeneous_monomials_of_weight(M, weight, check=False):
    if check:
        if not check_generators_weights(M):
            raise ValueError("the generators have different weights")
        k = M.generators()[0][0]
        if not k.is_power_of(weight):
            raise ValueError("the given weight should be the a power of the weight of the generators")
    gen_list = M.gen_forms()
    comb = combinations_with_replacement(gen_list, weight/2)
    monomials = []
    for tuple in list(comb):
        monomials.append( prod(tuple))
    return monomials

def initialize_matrix(monomials):
    if len(monomials) == 0:
        raise ValueError
    parent = monomials[0].parent()
    matrix_data = []
    for f in monomials:
        matrix_data.append(f.coefficients([0..parent.sturm_bound()]))
    return Matrix(matrix_data).transpose()

def relations(M, weight, check=False):
    monomials = homogeneous_monomials_of_weight(M, weight, check)
    A = initialize_matrix(monomials)
    ker = A.right_kernel()
    return ker.matrix()

def check_generators_relations(M, relation_matrix, weight):
    monomials = homogeneous_monomials_of_weight(M, weight)
    for row in range(0, relation_matrix.nrows()):
        rel = 0
        for col in range(0, relation_matrix.ncols()):
            rel += relation_matrix[row, col]*monomials[col]
        if not rel.is_zero():
            return False
    return True


if __name__=='__main__':
    o = open('result.txt','w')

    text = """
The rows of the column gives the coefficients of the relations.
For example, if G = Gamma0(6), the relation matrix is

    [ 0  0  1 -1  2 11]
    
so a polynomial relation in weight 4 is given by:

    g1*g3 - g2^2 + 2*g2*g3 + 11*g3^2 = 0

where g1, g2 and g3 are the generators of M(Gamma0(6))
\n"""

    print("Writing on file 'result.txt' ")
    o.write(text)
    G = Gamma0
    maxN = 25
    print("Maximum level: %s"%(maxN))
    verbose = False
    for N in range(1,maxN+1):
        M = ModularFormsRing(G(N))
        if verbose:
            print("Group: %s, Weight: %s"%(G(N), 2*k))
        if check_generators_weights(M):
            k = weights_of_generators(M)[0]
            relation_matrix = relations(M, 2*k)
            o.write("Group: %s, Weight: %s \nRelation matrix: \n"%(G(N), 2*k))
            o.write(relation_matrix.str() + "\n\n")
            if not check_generators_relations(M, relation_matrix, 2*k):
                print("Problem detected ---> Group: %s, Weight: %s !"%(M, G(N), 2*k))
    print("End of computations")
    o.close()

