from functools import wraps

#Decorators

def make_implicit(func):

    @wraps(func)
    def implicit_func(x,y,*args):
        return func(x,*args)-y

    return implicit_func