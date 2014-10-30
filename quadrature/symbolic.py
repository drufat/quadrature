import sympy as sp
import itertools
from quadrature.integral import quadrature
import functools

def terms(variables, order):
    '''
    Compute all the possible combinations of terms up to specific `order`
    that may occur in a quadrature formula.

    >>> terms('xy', 2) == [
    ... '',
    ... 'x', 'y',
    ... 'xx', 'xy', 'yy']
    True

    >>> terms('xy', 3) == [
    ... '',
    ... 'x', 'y',
    ... 'xx', 'xy', 'yy',
    ... 'xxx', 'xxy', 'xyy', 'yyy']
    True

    >>> terms('xyz', 2) == [
    ... '',
    ... 'x', 'y', 'z',
    ... 'xx', 'xy', 'xz', 'yy', 'yz', 'zz']
    True

    >>> terms('xyz', 3) == [
    ... '',
    ... 'x', 'y', 'z',
    ... 'xx', 'xy', 'xz', 'yy', 'yz', 'zz',
    ... 'xxx', 'xxy', 'xxz', 'xyy', 'xyz', 'xzz', 'yyy', 'yyz', 'yzz', 'zzz']
    True
    '''
    total = []
    for o in range(order+1):
        total += [''.join(c) for c in
                    itertools.combinations_with_replacement(variables, o)]
    return total

(s, t, u,
 x0, x1, x2, x3,
 y0, y1, y2, y3,
 z0, z1, z2, z3) = sp.symbols('s t u x0:4 y0:4 z0:4')

def symbolic_line(n,m):
    x = (1-s)*x0 + s*x1
    y = (1-s)*y0 + s*y1
    return sp.integrate(x**n * y**m, (s, 0, 1))

def symbolic_triangle(n,m):
    x = (1-s-t)*x0 + s*x1 + t*x2
    y = (1-s-t)*y0 + s*y1 + t*y2
    return 2*sp.integrate(sp.integrate(x**n * y**m, (t, 0, 1-s)), (s, 0, 1))

def symbolic_tetrahedron(n,m,p):
    x = (1-s-t-u)*x0 + s*x1 + t*x2 + u*x3
    y = (1-s-t-u)*y0 + s*y1 + t*y2 + u*y3
    z = (1-s-t-u)*z0 + s*z1 + t*z2 + u*z3
    return 6*sp.integrate(sp.integrate(sp.integrate(
                x**n * y**m * z**p, (u, 0, 1-s-t)), (t, 0, 1-s)), (s, 0, 1))

def tolatex(expression):
    p = sp.factor(expression.as_poly())
    return sp.latex(p)

def flatten(lst):
    '''
    Flatten nested list or tuples.
    >>> flatten([2, [2, [4, 5, [7], [2, [6, 2, 6, [6], 4]], 6]]])
    [2, 2, 4, 5, 7, 2, 6, 2, 6, 6, 4, 6]
    '''

    isnotflat = lambda x: isinstance(x, (list,tuple))
    def f(x,y):
        if isnotflat(y):
            return x+flatten(y)
        return x + [y]
    return functools.reduce(f, lst, [])

def symbolic_quadrature(term, verts):
    '''
    >>> symbolic_quadrature('x', [0,1,2])
    x0/3 + x1/3 + x2/3
    >>> symbolic_quadrature('xy', [0,1])
    x0*y0/3 + x0*y1/6 + x1*y0/6 + x1*y1/3
    '''

    expression, multiple = quadrature(term, verts)
    total = 0
    for term, coeff in expression:
        V = coeff
        for v in  (sp.symbols(str(v)+str(i)) for v,i in term):
            V *= v
        total += V
    return total*multiple
