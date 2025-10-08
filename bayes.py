from factor import Factor, events, marginalize, multiply_factors

class BayesianNetwork:
    """Represents a Bayesian network by its factors, i.e. the conditional probability tables (CPTs).

    Parameters
    ----------
    factors : list[factor.Factor]
        The factors of the Bayesian network
    domains : dict[str, list[str]]
        A dictionary mapping each variable to its possible values
    """

    def __init__(self, factors, domains):
        self.factors = factors
        self.domains = domains
        self.variables = set()
        for factor in self.factors:
            self.variables = self.variables | set(factor.variables)

    def __str__(self):
        return "\n\n".join([str(factor) for factor in self.factors])


def eliminate(bnet, variable):
    """Eliminates a variable from the Bayesian network.

    By "eliminate", we mean that the factors containing the variable are multiplied,
    and then the variable is marginalized (summed) out of the resulting factor.

    Parameters
    ----------
    variable : str
        the variable to eliminate from the Bayesian network

    Returns
    -------
    BayesianNetwork
        a new BayesianNetwork, equivalent to the current Bayesian network, after
        eliminating the specified variable
    """

    with_var = [f for f in bnet.factors if variable in f.variables]
    without_var = [f for f in bnet.factors if variable not in f.variables]

    product = multiply_factors(with_var, bnet.domains)
    reduced = marginalize(product, variable)
    new_factors = without_var + [reduced]
    return BayesianNetwork(new_factors, bnet.domains)


def compute_marginal(bnet, vars):
    """Computes the marginal probability over the specified variables.

    This method uses variable elimination to compute the marginal distribution.

    Parameters
    ----------
    vars : set[str]
        the variables that we want to compute the marginal over
    """

    vars = set(vars)
    bn = BayesianNetwork(bnet.factors[:], bnet.domains)
    for v in (bn.variables - vars):
        bn = eliminate(bn, v) # eliminate other variables

    result = multiply_factors(bn.factors, bnet.domains)
    return result
    
    
def compute_conditional(bnet, event, evidence):
    """Computes the conditional probability of an event given the evidence event."""

    scope_e = set(event) | set(evidence)
    joint = compute_marginal(bnet, scope_e)
    num = joint[{**event, **evidence}]
    denom_factor = compute_marginal(bnet, set(evidence))
    denom = denom_factor[evidence]

    return num / denom