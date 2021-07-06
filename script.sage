
r"""
SageMath script that execute the functions defined in `relation.sage`.

The results are saved in a text file.
"""


load('relations.sage')

if __name__=='__main__':
    name_of_text_file = 'result' #Choose the name of the results text file
    o = open(name_of_text_file + '.txt','w')

    text = """
The rows of the relation matrix gives the coefficients of the relations.
For example, if G = Gamma0(6), the relation matrix is

    [ 0  0  1 -1  2 11]
    
so a polynomial relation in weight 4 is given by:

    g1*g3 - g2^2 + 2*g2*g3 + 11*g3^2 = 0

where g1, g2 and g3 are generators for M(Gamma0(6))
\n"""

    print("Writing on file 'result.txt' ")
    o.write(text)
    G = Gamma0 # The code currently only works for Gamma0
    maxN = 25 # Maximum level
    print("Maximum level: %s"%(maxN))
    verbose = False # Set this to True if you want more console outputs
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

