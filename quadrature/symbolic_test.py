from quadrature import *

def test_symbolic_line():
    for term in terms('xy', 3):
        assert( symbolic_line(term.count('x'),term.count('y')) ==
                symbolic_quadrature(term, (0,1)) )

def test_symbolic_triangle():
    for term in terms('xy', 2):
        assert( symbolic_triangle(term.count('x'),term.count('y')) ==
                symbolic_quadrature(term, (0,1,2)) )

def test_symbolic_tetrahedron():
    for term in terms('xyz', 1):
        assert( symbolic_tetrahedron(term.count('x'),
                                     term.count('y'),
                                     term.count('z')) ==
                symbolic_quadrature(term, (0,1,2,3)) )
