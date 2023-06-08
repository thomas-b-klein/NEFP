from functools import wraps

#Decorators

def make_implicit(func):

    @wraps(func)
    def implicit_func(x,y):
        return func(x)-y

    return implicit_func