from deap.gp import Primitive
from sympy import parse_expr

infix_map = {
    'add_2': '+',
    'sub_2': '-',
    'multiply': '*',
    'protect_divide': '/',
}


def gene_to_string(gene):
    string = ""
    stack = []
    for node in gene:
        stack.append((node, []))
        while len(stack[-1][1]) == stack[-1][0].arity:
            prim, args = stack.pop()
            # string = prim.format(*args)
            if type(prim) is Primitive:
                string = '('
                if prim.name == 'analytical_quotient':
                    string += f'{args[0]}/sqrt(1+{args[1]}*{args[1]})'
                elif prim.name not in infix_map:
                    string += prim.format(*args)
                else:
                    string += args[0]
                    for a in args[1:]:
                        string += f'{infix_map[prim.name]}{a}'
                string += ')'
            else:
                string = prim.name
            if len(stack) == 0:
                break  # If stack is empty, all nodes should have been seen
            stack[-1][1].append(string)
    return string


def multigene_gp_to_string(label, regr):
    pipes = regr.pipelines
    cur_gene = None
    for i, g in enumerate(regr.best_pop):
        coef = pipes[label]['Ridge'].coef_[i]
        mean = pipes[label]['Scaler'].mean_[i]
        var = pipes[label]['Scaler'].scale_[i]
        if cur_gene is None:
            cur_gene = coef * (parse_expr(gene_to_string(g)) - mean) / var
        else:
            cur_gene += coef * (parse_expr(gene_to_string(g)) - mean) / var
    cur_gene += pipes[label]['Ridge'].intercept_
    return cur_gene
