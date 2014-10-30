import itertools
import math
import functools
import collections
import fractions

def nCk(n,k):
    r'''

    Compute the binomial coefficients :math:`\binom nk = \frac{n!}{k!(n-k)!}`

    >>> nCk(0,0)
    1
    >>> nCk(1,0)
    1
    >>> nCk(4, 2)
    6
    >>> [[nCk(n,k) for k in range(n+1)] for n in range(7)] == [
    ...          [1],
    ...        [1, 1],
    ...       [1, 2, 1],
    ...      [1, 3, 3, 1],
    ...     [1, 4, 6, 4, 1],
    ...    [1, 5, 10, 10, 5, 1],
    ...  [1, 6, 15, 20, 15, 6, 1]                              ]
    True
    '''

    f = math.factorial
    return f(n) // f(k) // f(n-k)

def tostr(expr):
    '''
    Convert quadrature result to a human readable string.

    >>> q =  ([((('x', 0),), 1), ((('x', 1),), 1)], fractions.Fraction(1, 2))
    >>> tostr(q)
    '1/2 (x0 + x1)'
    '''
    expression, multiple = expr
    terms = []
    for t,c in expression:
        s = ''
        if c > 1:
            s += str(c)
            s += ' '
        s += ''.join(map(str, sum(t,())))
        terms.append(s)
    s = ' + '.join(terms)
    if not multiple == 1:
        s = '{}/{} ({})'.format(multiple.numerator, multiple.denominator, s)
    return s

def quadrature(term, verts):
    '''

    Compute the quadrature of a homogenous polynomial given by term='xxy'
    over the vertices verts=(0,1,2) of a simplex.

    >>> quadrature('', (0,1))
    ([((), 1)], Fraction(1, 1))
    >>> quadrature('x', (0,1))
    ([((('x', 0),), 1), ((('x', 1),), 1)], Fraction(1, 2))
    >>> quadrature('x', (0,1,2))
    ([((('x', 0),), 1), ((('x', 1),), 1), ((('x', 2),), 1)], Fraction(1, 3))
    >>> quadrature('xx', (0,1)) == (
    ...    [((('x', 0), ('x', 0)), 1),
    ...     ((('x', 0), ('x', 1)), 1),
    ...     ((('x', 1), ('x', 1)), 1)], fractions.Fraction(1, 3))
    True
    >>> quadrature('xy', (0,1,2)) == (
    ...    [((('x', 0), ('y', 0)), 2),
    ...      ((('x', 0), ('y', 1)), 1),
    ...      ((('x', 0), ('y', 2)), 1),
    ...      ((('x', 1), ('y', 0)), 1),
    ...      ((('x', 1), ('y', 1)), 2),
    ...      ((('x', 1), ('y', 2)), 1),
    ...      ((('x', 2), ('y', 0)), 1),
    ...      ((('x', 2), ('y', 1)), 1),
    ...      ((('x', 2), ('y', 2)), 2)], fractions.Fraction(1, 12))
    True
    '''

    k = len(verts)-1    # dimension of simplex
    q = len(term)       # order of term
    P = itertools.permutations(term)
    V = itertools.combinations_with_replacement(verts, q)
    Xp = itertools.product(P, V)
    X = (tuple(sorted(zip(*x))) for x in Xp)
    counter = collections.Counter(X)

    gcd = functools.reduce(fractions.gcd, counter.values())

    return (sorted([(c, counter[c] // gcd) for c in counter ]),
            fractions.Fraction(gcd, math.factorial(q)*nCk(k+q, q)))

def q(term, verts):
    return tostr(quadrature(term, verts))

def test_line_segment():
    '''

    >>> def q(t): return tostr(quadrature(t, (0,1)))
    >>> q('x')
    '1/2 (x0 + x1)'
    >>> q('xx')
    '1/3 (x0x0 + x0x1 + x1x1)'
    >>> q('xxx')
    '1/4 (x0x0x0 + x0x0x1 + x0x1x1 + x1x1x1)'
    >>> q('xy')
    '1/6 (2 x0y0 + x0y1 + x1y0 + 2 x1y1)'
    >>> q('xxy')
    '1/12 (3 x0x0y0 + x0x0y1 + 2 x0x1y0 + 2 x0x1y1 + x1x1y0 + 3 x1x1y1)'
    '''
    pass

def test_triangle():
    '''

    >>> def q(t): return tostr(quadrature(t, (0,1,2)))
    >>> q('x')
    '1/3 (x0 + x1 + x2)'
    >>> q('xx')
    '1/6 (x0x0 + x0x1 + x0x2 + x1x1 + x1x2 + x2x2)'
    >>> q('xxx')
    '1/10 (x0x0x0 + x0x0x1 + x0x0x2 + x0x1x1 + x0x1x2 + x0x2x2 + x1x1x1 + x1x1x2 + x1x2x2 + x2x2x2)'
    >>> q('xy')
    '1/12 (2 x0y0 + x0y1 + x0y2 + x1y0 + 2 x1y1 + x1y2 + x2y0 + x2y1 + 2 x2y2)'
    '''
    pass

def test_tetrahedron():
    '''

    >>> def q(t): return tostr(quadrature(t, (0,1,2,3)))
    >>> q('x')
    '1/4 (x0 + x1 + x2 + x3)'
    >>> q('xx')
    '1/10 (x0x0 + x0x1 + x0x2 + x0x3 + x1x1 + x1x2 + x1x3 + x2x2 + x2x3 + x3x3)'
    >>> q('xxx')
    '1/20 (x0x0x0 + x0x0x1 + x0x0x2 + x0x0x3 + x0x1x1 + x0x1x2 + x0x1x3 + x0x2x2 + x0x2x3 + x0x3x3 + x1x1x1 + x1x1x2 + x1x1x3 + x1x2x2 + x1x2x3 + x1x3x3 + x2x2x2 + x2x2x3 + x2x3x3 + x3x3x3)'
    >>> q('xy')
    '1/20 (2 x0y0 + x0y1 + x0y2 + x0y3 + x1y0 + 2 x1y1 + x1y2 + x1y3 + x2y0 + x2y1 + 2 x2y2 + x2y3 + x3y0 + x3y1 + x3y2 + 2 x3y3)'
    '''
    pass
