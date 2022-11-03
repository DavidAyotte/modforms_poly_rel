# Script counting the number of occurrences of a weight appearing in the list of
# generators returned by SageMath

max_level = 15
prec = max_level + 1  # q-expansion precision
group = Gamma0  # choose either Gamma0 or Gamma1
file_name = ""  # change this if you want a costum name for the text file containing the results
                # the default name is "{group}_max_{max_level}"

# The script compute a dictionnary named "result" where
#  * the keys are the level up to max_level
#  * the values are lists of the form [forms_weights, count_weight_occurrences]
#    where forms_weights is a list of pairs (f, k) with f a generator and k its
#    weight
result = {}
for N in range(1, max_level + 1):
    M = ModularFormsRing(group(N))
    gens_list = M.gens()
    weights = [f.weight() for f in gens_list]
    forms_weights = [[e[0].qexp(prec), e[1]] for e in zip(gens_list, weights)]
    unique_weights = list(set(weights))
    count_weight_occurrences = [(k, weights.count(k)) for k in unique_weights]

    result[N] = [forms_weights, count_weight_occurrences]

# output the results in a file:
if group == Gamma0:
    group_str = "Gamma0"
elif group == Gamma1:
    group_str = "Gamma1"

if len(file_name) == 0:
    file_name = group_str + '_max_' + str(max_level) + '.txt'

fi = open(file_name, 'w')

fi.write(f"Chosen group: {group_str}\n")
fi.write(f"Maximum level: {max_level}\n")
fi.write(f"Precision: {prec}\n")

E4 = ModularForms(1, 4).gen(0)
E6 = ModularForms(1, 6).gen(0)
fi.write(f"Weight 4 and 6 Eisenstein series:\n  E4 = {E4.qexp(prec)}\n  E6 = {E6.qexp(prec)}\n")

for N in result.keys():
    forms_list_str = str(result[N][0]).replace("], [", "]\n     [")
    fi.write(f"\n{group_str}({N}):\n    {result[N][1]}\n    {forms_list_str}")

fi.close()
