import itertools
from collections import defaultdict

class Factor:
    def __init__(self, variables, values):
        self.variables = variables
        self.values = values

    def __getitem__(self, event):
        key = []
        for var in self.variables:
            if var not in event:
                raise KeyError(f"Variable {var} not found in given event.")
            key.append(event[var])
        if tuple(key) in self.values:
            return self.values[tuple(key)]
        else:
            raise KeyError(f"No value assigned to event {event}.")

    def __str__(self):
        result = f"{self.variables}:"
        for event, value in self.values.items():
            result += f"\n  {event}: {value}"
        return result

    __repr__ = __str__


def events(vars, domains):
    dom_lists = [domains[v] for v in vars]
    result = []

    # for every combination of values
    for values in itertools.product(*dom_lists):
        event = {}
        for v, val in zip(vars, values):
            event[v] = val
        result.append(event)
    return result


def marginalize(f, var):
    idx = f.variables.index(var) # index of variable in factor
    new_vars = [v for v in f.variables if v != var]
    new_vals = defaultdict(float)

    for event, value in f.values.items(): # iterating through each event
        reduced = event[:idx] + event[idx+1:] # slice var from event
        new_vals[reduced] += value

    return Factor(new_vars, new_vals)


def multiply_factors(factors, domains):
    new_vars = []
    for f in factors:
        for v in f.variables:
            if v not in new_vars:
                new_vars.append(v)

    e = events(new_vars, domains)
    new_vals = {}
    for ev in events(new_vars, domains):
        prod_val = 1.0
        for f in factors:
            prod_val *= f[ev]  # uses only f.variables internally
        key = tuple(ev[v] for v in new_vars)
        new_vals[key] = prod_val


    return Factor(new_vars, new_vals)

