from bayes import BayesianNetwork
from factor import Factor


class FamilyMember:
    """A single member of a family tree."""

    def __init__(self, name, sex, mother, father):
        """
        Parameters
        ----------
        name : str
            The name of the family member.
        sex : str
            The sex of the family member ("male" or "female")
        mother : FamilyMember
            The mother of the family member (or None if unknown)
        father : FamilyMember
            The father of the family member (or None if unknown)
        """

        self.name = name
        self.sex = sex
        self.mother = mother
        self.father = father

    def get_name(self):
        """Returns the name of the family member."""
        return self.name

    def get_sex(self):
        """Returns the sex of the family member."""
        return self.sex


class Male(FamilyMember):
    """A male family member."""

    def __init__(self, name, mother=None, father=None):
        super().__init__(name, "male", mother, father)


class Female(FamilyMember):
    """A female family member."""

    def __init__(self, name, mother=None, father=None):
        super().__init__(name, "female", mother, father)


def romanoffs():
    """A simple example of a family, using four members of the Russian royal family (the Romanoffs)."""
    alexandra = Female("alexandra")
    nicholas = Male("nicholas")
    alexey = Male("alexey", mother=alexandra, father=nicholas)
    anastasia = Female("anastasia", mother=alexandra, father=nicholas)
    return alexandra, nicholas, alexey, anastasia


def create_variable_domains(family):
    """Creates a dictionary mapping each variable to its domain, for the hemophilia network.

    For each family member, we create either 3 or 4 variables (3 if they’re male, 4 if they’re female).
    If N is the name of the family member, then we create the following variables:
        M_N: N’s maternally inherited gene
        P_N: N’s paternally inherited gene (if N is female)
        G_N: the genotype of N
        H_N: whether N has hemophilia

    The variables should be mapped to the following domains:
        - M_N: ['x', 'X']
        - P_N: ['x', 'X']
        - G_N: ['xx', 'xX', 'XX']
        - H_N: ['-', '+']

    Parameters
    ----------
    family : list[FamilyMember]
        the list of family members

    Returns
    -------
    dict[str, list[str]]
        a dictionary mapping each variable to its domain (i.e. its possible values)
    """
    result = dict()
    for person in family:
        name = person.get_name()

        result["M_" + name] = ['x', 'X']
        result["H_" + name] = ['-', '+']

        if person.get_sex() == "female": # if person is female
            result["P_" + name] = ['x', 'X']

            result["G_" + name] = ['xx', 'xX', 'XX']
        else: # if person is male
            result["G_" + name] = ['xy', 'Xy']

    return result

def create_hemophilia_cpt(person):
    """Creates a conditional probability table (CPT) specifying the probability of hemophilia, given one's genotype.

    Parameters
    ----------
    person : FamilyMember
        the family member whom the CPT pertains to

    Returns
    -------
    Factor
        a Factor specifying the probability of hemophilia, given one's genotype
    """
    name = person.get_name()
    H_N = "H_" + name
    G_N = "G_" + name

    if person.get_sex() == "female":
        genotypes = ['xx', 'xX', 'XX']
        hemo = {'xx': '-', 'xX': '-', 'XX': '+'}
    else:
        genotypes = ['xy', 'Xy']
        hemo = {'xy': '-', 'Xy': '+'}

    values = {}
    for g in genotypes:
        values[(g, '-')] = 1.0 if hemo[g] == '-' else 0.0
        values[(g, '+')] = 1.0 if hemo[g] == '+' else 0.0

    return Factor([G_N, H_N], values)



def create_genotype_cpt(person):
    """Creates a conditional probability table (CPT) specifying the probability of a genotype, given one's inherited genes.

    Parameters
    ----------
    person : FamilyMember
        the family member whom the CPT pertains to

    Returns
    -------
    Factor
        a Factor specifying the probability of a genotype, given one's inherited genes
    """

    name = person.get_name()
    G_N = "G_" + name
    M_N = "M_" + name

    values = {}

    if person.get_sex() == "female":
        P_N = "P_" + name
        for m in ['x', 'X']:
            for p in ['x', 'X']:
                if (m == p == 'x'):
                    g = 'xx'
                elif (m == p == 'X'):
                    g = 'XX'
                else:
                    g = 'xX'
                for g_vars in ['xx', 'xX', 'XX']:
                    values[(m, p, g_vars)] = 1.0 if g_vars == g else 0.0
        return Factor([M_N, P_N, G_N], values)
    else:
        for m in ['x', 'X']:
            g = m + 'y'
            for g_vars in ['xy', 'Xy']:
                values[(m, g_vars)] = 1.0 if g_vars == g else 0.0

        return Factor([M_N, G_N], values)


def create_maternal_inheritance_cpt(person):
    """Creates a conditional probability table (CPT) specifying the probability of the gene inherited from one's mother.

    Parameters
    ----------
    person : FamilyMember
        the family member whom the CPT pertains to

    Returns
    -------
    Factor
        a Factor specifying the probability of the gene inherited from the family member's mother.
    """
    name = person.get_name()
    M_N = "M_" + name
    values = {}

    if person.mother == None:
        values = {('x',): 29999/30000, ('X',): 1/30000}
        return Factor([M_N], values)

    mom = person.mother
    G_A = "G_" + mom.get_name()

    for ga in ['xx', 'xX', 'XX']:
        if ga == 'xx':
            values[(ga, 'x')] = 1.0
            values[(ga, 'X')] = 0.0
        elif ga == 'xX':
            values[(ga, 'x')] = 0.5
            values[(ga, 'X')] = 0.5
        else:  # 'XX'
            values[(ga, 'x')] = 0.0
            values[(ga, 'X')] = 1.0

    return Factor([G_A, M_N], values)


def create_paternal_inheritance_cpt(person):
    """Creates a conditional probability table (CPT) specifying the probability of the gene inherited from one's father.

    Parameters
    ----------
    person : FamilyMember
        the family member whom the CPT pertains to

    Returns
    -------
    Factor
        a Factor specifying the probability of the gene inherited from the family member's father.
    """
    name = person.get_name()
    P_N = "P_" + name
    values = {}

    if person.get_sex == "male":
        values = {('y',): 1.0}
        return Factor([P_N], values)

    if person.father == None:
        values = {('x',): 29999/30000, ('X',): 1/30000}
        return Factor([P_N], values)

    dad = person.father
    G_A = "G_" + dad.get_name()

    for ga in ['xy', 'Xy']:
        if ga == 'xy':
            values[(ga, 'x')] = 1.0
            values[(ga, 'X')] = 0.0
        else:  # 'Xy'
            values[(ga, 'x')] = 0.0
            values[(ga, 'X')] = 1.0

    return Factor([G_A, P_N], values)


def create_family_bayes_net(family):
    """Creates a Bayesian network that models the genetic inheritance of hemophilia within a family.

    Parameters
    ----------
    family : list[FamilyMember]
        the members of the family

    Returns
    -------
    BayesianNetwork
        a Bayesian network that models the genetic inheritance of hemophilia within the specified family
    """
    domains = create_variable_domains(family)
    cpts = []
    for person in family:
        if person.get_sex() == "female":
            cpts.append(create_paternal_inheritance_cpt(person))
        cpts.append(create_maternal_inheritance_cpt(person))
        cpts.append(create_genotype_cpt(person))
        cpts.append(create_hemophilia_cpt(person))
    return BayesianNetwork(cpts, domains)
