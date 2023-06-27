from functools import wraps

#Decorators

def make_implicit(func):

    @wraps(func)
    def implicit_func(x,y,*args):
        return func(x,*args)-y

    return implicit_func

# Useful Functions

def within(x, bounds, include_bounds = True):
    if include_bounds: return True if min(bounds) <= x <= max(bounds) else False
    if not include_bounds: return True if min(bounds) < x < max(bounds) else False