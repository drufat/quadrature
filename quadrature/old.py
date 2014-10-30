import numpy as np
import operator
import functools
import itertools
import fractions

def coords_symbolic(coords, order=2):
    """

    Generate all the terms in the polynomial up to order `order`.

    >>> coords_symbolic("x", order=1)
    [(), ('x',)]
    >>> coords_symbolic("xy", order=1)
    [(), ('x',), ('y',)]
    >>> coords_symbolic('xy', order=2)
    [(), ('x',), ('y',), ('x', 'x'), ('x', 'y'), ('y', 'y')]
    """
    C = []
    for i in range(order + 1):
        C += itertools.combinations_with_replacement(coords, i)
    return C

def moments_symbolic(C, verts):
    """
    Compute the moments symbolically.

    Edge moments:

    >>> C = coords_symbolic("xy", order=2)
    >>> M = moments_symbolic(C, verts=(0,1))
    >>> (C[1], len(M[1]), M[1]) == \
        (('x',), 2, {(('x', 0),): 1, (('x', 1),): 1})
    True
    >>> (C[3], len(M[3]), M[3]) == \
        (('x', 'x'), 3, {(('x', 0), ('x', 1)): 1, (('x', 1), ('x', 1)): 1, (('x', 0), ('x', 0)): 1})
    True
    >>> (C[4], len(M[4]), M[4]) == \
        (('x', 'y'), 4, {(('x', 0), ('y', 1)): 1, (('x', 0), ('y', 0)): 2, (('x', 1), ('y', 1)): 2, (('x', 1), ('y', 0)): 1})
    True

    Triangle moments:

    >>> C = coords_symbolic("xy", order=2)
    >>> M = moments_symbolic(C, verts=(0,1,2))
    >>> (C[1], len(M[1]), M[1]) == \
        (('x',), 3, {(('x', 0),): 1, (('x', 1),): 1, (('x', 2),): 1})
    True
    >>> (C[3], len(M[3]), M[3]) == \
        (('x', 'x'), 6, {(('x', 1), ('x', 1)): 1, (('x', 0), ('x', 2)): 1, (('x', 2), ('x', 2)): 1, (('x', 0), ('x', 0)): 1, (('x', 0), ('x', 1)): 1, (('x', 1), ('x', 2)): 1})
    True

    """
    M = []
    for c in C:
        P = itertools.permutations(c)
        V = itertools.combinations_with_replacement(verts, len(c))
        X = itertools.product(P, V)
        X = [tuple(sorted(zip(*x))) for x in X]
        D = {}
        for x in X:
            if x in D: D[x] += 1
            else: D[x] = 1
        gcd = functools.reduce(fractions.gcd, list(D.values()))
        for k in D: D[k] /= gcd
        M.append(D)
    return M

def moments_eval(M, coord):
    """
    Evaluate the moments/averages.
    coord[simplices,vertices,dimensions]
    """
    X = [np.zeros(coord.shape[0]) for i in range(len(M))]
    for i, m in enumerate(M):
        S = sum(m.values())
        for k in m:
            x = np.ones(coord.shape[0])
            if k:
                for d, v in k:
                    x *= coord[:, v, d]
            X[i] += m[k] * x / S
    return X

def measure_eval(coord):
    """
    Evaluate the area.
    coord[simplices,vertices,dimensions]
    """
    def cross(u, v):
        return u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0]
    A = np.zeros(coord.shape[0])
    n = coord.shape[1]
    for i in range(n):
        A += 0.5 * cross(coord[:, i, :], coord[:, (i + 1) % n, :])
    return A

def coord_eval(C, coord):
    """
    coord[xyz]
    """
    R = []
    for c in C:
        p = functools.reduce(operator.mul,
                             coord[:, c].T,
                             np.ones(coord.shape[0]))
        R.append(p)
    return R

class Product(object):
    '''
    Compute the product of two polynomials.

    >>> C = coords_symbolic("xy", order=1)
    >>> prod = Product(C)

    1*1=1
    >>> prod([1,0,0], [1,0,0]) == [1,0,0]
    True

    (1+x)*(1+y)=(1+x+y)
    >>> prod([1,1,0], [1,0,1]) == [1,1,1]
    True

    (1+x+y)*(1+y)=(1+x+2y)
    >>> prod([1,1,1], [1,0,1]) == [1,1,2]
    True

    (x)*(y)=(0)
    >>> prod([0,0,1], [0,0,1]) == [0,0,0]
    True

    '''

    def __init__(self, C):
        index = {}
        for i, c in enumerate(C):
            index[c] = i
        coeff = []
        for l, r in itertools.product(C, C):
            p = tuple(sorted(l + r))
            if p in index:
                coeff.append((index[p], index[l], index[r]))
        coeff = sorted(coeff)
        cauchi = [[] for c in C]
        for p, l, r in coeff:
            cauchi[p].append((l, r))
        self.coeff = cauchi

    def __call__(self, A, B):
        def prod(tup):
            l, r = tup
            return A[l] * B[r]
        P = []
        for c in self.coeff:
            P.append(functools.reduce(operator.add, list(map(prod, c))))
        return P

class Derivative(object):
    '''
    Compute the derivative of a polynomial

    >>> C = coords_symbolic("xy", order=2)
    >>> dx = Derivative(C, 'x')
    >>> dy = Derivative(C, 'y')

    d_x(1+x+y+xx+xy+yy)=(1+2x+y)
    >>> dx([1, 1, 1, 1, 1, 1]) == [1, 2, 1, 0, 0, 0]
    True

    d_y(1+x+y+xx+xy+yy)=(1+x+2y)
    >>> dy([1, 1, 1, 1, 1, 1]) == [1, 1, 2, 0, 0, 0]
    True

    '''

    def __init__(self, C, x):
        index = {}
        for i, c in enumerate(C):
            index[c] = i

        self.derivative = []
        for c in C:
            if x not in c:
                continue
            cc = list(c); cc.remove(x); cc = tuple(cc)
            i = index[cc]
            j = index[c]
            k = c.count(x)
            self.derivative += [(i, j, k)]

    def __call__(self, A):

        D = [0 * x for x in A]
        for (i, j, k) in self.derivative:
            D[i] = k * A[j]
        return D

