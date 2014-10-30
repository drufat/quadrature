Computing Quadratures of Polynomial Functions over k-Simplices :math:`\triangle^k` in :math:`\mathbb{R}^n`
===========================================================================================================

Contents:

.. toctree::
    :maxdepth: 2

    symbolic

Suppose we are given a k-Simplex :math:`\triangle^k` in :math:`\mathbb{R}^n`

.. math::

    \triangle^k=\left\{ \sum_{i=0}^{k}\lambda_{i}\mathbf{x}_{i}\vert\lambda_{i}\ge0,\,0\leq i\le k,\,\sum_{i=0}^{k}\lambda_{i}=1\right\}

where :math:`\mathbf{x}_i` represents the i\ :sup:`th` vertex of the simplex.

.. math::

    \mathbf{x}_{1}=\left(\begin{array}{c}
    x_{1}\\
    y_{1}\\
    z_{1}\\
    \vdots
    \end{array}\right)\qquad\mathbf{x}_{2}=\left(\begin{array}{c}
    x_{2}\\
    y_{2}\\
    z_{2}\\
    \vdots
    \end{array}\right)\dots\qquad\mathbf{x}_{i}=\left(\begin{array}{c}
    x_{i}\\
    y_{i}\\
    z_{i}\\
    \vdots
    \end{array}\right)\dots\qquad\mathbf{x}_{k}=\left(\begin{array}{c}
    x_{k}\\
    y_{k}\\
    z_{k}\\
    \vdots
    \end{array}\right)

We are interested in finding an explicit expression for the quadrature of a homogeneous polynomial
over such a simplex:

.. math::

    \int_{\triangle_{k}}x^{\alpha}y^{\beta}\dots

and the moments, which are the above expression divided by the volume of the simplex:

.. math::
    \langle x^\alpha y^\beta \dots \rangle = \frac{1}{\text{vol} (\triangle_k)}  \int_{\triangle_{k}}x^{\alpha}y^{\beta}\dots

We have implemented in python a simple algorithm that computes exactly the above expressions:


.. autofunction:: quadrature.nCk

.. autofunction:: quadrature.quadrature


Description of Algorithm
===================================

Our algorithm is general, but for simplicity we will just demonstrate it for the specific case of a triangle
in two dimensions. Assume you have a triangle in the plane defined by its three vertices at
:math:`(x_0,y_0)`, :math:`(x_1,y_1)`, and :math:`(x_2,y_2)` - it is a 2-simplex so :math:`k=2`. We will compute
the moment of the q-homoqenous polynomial :math:`\langle xy \rangle` with :math:`q=2`.

.. plot:: fig/triangle.py

The inputs to the algorithm are two tuples - the first is :math:`\text{term}=(x,y)` describing the polynomial and the
other is :math:`\text{verts}=(0,1,2)` describing the labels for the vertices. With this notation, the orders of the simplex
:math:`k` and the degree of the polynomial :math:`q` are given by

.. math::

    k &= \text{len}(\text{verts})-1 = 2 \\
    q &= \text{len}(\text{term}) = 2\\

Now take all the permutations of the :math:`\text{term}` tuple.

.. math::

    P = \mathcal{P}(\text{term}) = (xy, yx)

Compute all the combinations (with repetitions allowed) of size :math:`q=2` of the coordinates of the vertices

.. math::

    V = \mathcal{C}^q(\text{verts}) = (00, 01, 02, 11, 12, 22)

Take the product of :math:`V` and :math:`P`

.. math::

    P \times V = \left (  (xy,00), (xy,01),(xy,02),(xy,11),(xy,12),(xy,22),
                          (yx,00),(yx,01),(yx,02),(yx,11),(yx,12),(yx,22) \right)


For each element in the tuple zip together the variable with the simplex vertices (e.g. :math:`(xy,00) \mapsto x_0y_0, \quad (xy,01) \mapsto x_0y_1` and so on):

.. math::

    \text{zip}(P \times V) = \left (  x_0 y_0, x_0 y_1, x_0 y_2, x_1 y_1, x_1y_2, x_2y_2,
                          y_0x_0, y_0x_1,y_0x_2,y_1x_1,y_1x_2,y_2x_2 \right)

Treat each element as a multiplication and sum all the multiples together

.. math::

    \sum(\text{zip}(P \times V)) = \left ( 2 x_0 y_0 +   x_0 y_1 +   x_0 y_2 +
                                             x_1 y_0 + 2 x_1 y_1 +   x_1 y_2 +
                                             x_2 y_0 +   x_2 y_1 + 2 x_2 y_2 \right)

Divide the result by :math:`\binom{k+q}{q} q!` to obtain the final expression

.. math::

    \frac{1}{\binom{k+q}{q} q!} \sum(\text{zip}(P \times V)) = \frac{1}{12} \left ( 2 x_0 y_0 +   x_0 y_1 +   x_0 y_2 +
                                             x_1 y_0 + 2 x_1 y_1 +   x_1 y_2 +
                                             x_2 y_0 +   x_2 y_1 + 2 x_2 y_2 \right)

which is indeed what one would obtain for :math:`\langle xy \rangle` by manual calculation.

Pseudocode
===========
In pseudocode the algorithm can be concisely given as:

.. math::

    & \text{def  quadrature}(\text{term},\,\text{verts}): \\
    & \qquad  k = \text{len}(\text{verts})-1 \\
    & \qquad q = \text{len}(\text{term}) \\
    & \qquad P = \mathcal{P}(\text{term}) \\
    & \qquad V = \mathcal{C}^q(\text{verts}) \\
    & \qquad \text{return} \quad \frac{1}{\binom{k+q}{q} q!} \sum(\text{zip}(P \times V))

Validation
============

To validate the algorithm we consider a triangle, a line segment, and a tetrahedron, and we compute a
number of low order moments by integrating explicitly. The Mathematica file where all the integration is done
is  provided :download:`for download here <quadratures.nb>`.


Triangle :math:`\triangle^2` in :math:`\mathbb{R}^2`
----------------------------------------------------

The moments of all the homogeneous polynomials over a triangle area given by

.. math::

   \langle x^{n}y^{m}\rangle=\frac{1}{A}\int_{\triangle}x^{n}y^{m}\,\mathbf{d}x\wedge\mathbf{d}y

where :math:`\mathbf{d}x\wedge\mathbf{d}y` is the area form for :math:`\mathbb{R}^2`
and :math:`A=\int_{\triangle}\mathbf{d}x\wedge\mathbf{d}y` is the total area of the triangle.

For example, :math:`\langle x \rangle` and :math:`\langle y \rangle` are the
coordinates of the center of mass of the triangle.

Change of variables :math:`x,y \mapsto s,t`


.. math::

    x	&=	(1-s-t)x_{0}+sx_{1}+tx_{2}=s(x_{1}-x_{0})+t(x_{2}-x_{0}) \\
    y	&=	(1-s-t)y_{0}+sy_{1}+ty_{2}=s(y_{1}-y_{0})+t(y_{2}-y_{0})


so that :math:`(s,t)=(0,0)` corresponds to vertex :math:`(x_{0},y_{0})`, :math:`(s,t)=(1,0)`
corresponds to vertex :math:`(x_{1},y_{1})`, and :math:`(s,t)=(0,1)` corresponds to vertex
:math:`(x_{2},y_{2})`. The area form transforms as follows:

.. math::

    dx\wedge dy	 &= \left((x_{1}-x_{0})(y_{2}-y_{0})-(x_{2}-x_{0})(y_{1}-y_{0})\right)\, ds\wedge dt \\
                 &= 2A\, ds\wedge dt

Thus we have an explicit expression for the quadrature

.. math::

    \langle x^{n}y^{m}\rangle=2\int_{0}^{1}ds\int_{0}^{1-s}dt\left(s(x_{1}-x_{0})+t(x_{2}-x_{0})\right)^{n}\left(s(y_{1}-y_{0})+t(y_{2}-y_{0})\right)^{m}dt

which can be evaluated by hand to obtain the first few moments:

.. math::

    \frac{1}{6} \left(x_{0}^{2} + x_{0} x_{1} + x_{0} x_{2} + x_{1}^{2} + x_{1} x_{2} + x_{2}^{2}\right) \\
    \frac{1}{10} \left(x_{0}^{3} + x_{0}^{2} x_{1} + x_{0}^{2} x_{2} + x_{0} x_{1}^{2} + x_{0} x_{1} x_{2} + x_{0} x_{2}^{2} + x_{1}^{3} + x_{1}^{2} x_{2} + x_{1} x_{2}^{2} + x_{2}^{3}\right) \\
    \frac{1}{12} \left(2 x_{0} y_{0} + x_{0} y_{1} + x_{0} y_{2} + x_{1} y_{0} + 2 x_{1} y_{1} + x_{1} y_{2} + x_{2} y_{0} + x_{2} y_{1} + 2 x_{2} y_{2}\right) \\
    \frac{1}{21} \left(x_{0}^{5} + x_{0}^{4} x_{1} + x_{0}^{4} x_{2} + x_{0}^{3} x_{1}^{2} + x_{0}^{3} x_{1} x_{2} + x_{0}^{3} x_{2}^{2} + x_{0}^{2} x_{1}^{3} + x_{0}^{2} x_{1}^{2} x_{2} + x_{0}^{2} x_{1} x_{2}^{2} + x_{0}^{2} x_{2}^{3} + x_{0} x_{1}^{4} + x_{0} x_{1}^{3} x_{2} + x_{0} x_{1}^{2} x_{2}^{2} + x_{0} x_{1} x_{2}^{3} + x_{0} x_{2}^{4} + x_{1}^{5} + x_{1}^{4} x_{2} + x_{1}^{3} x_{2}^{2} + x_{1}^{2} x_{2}^{3} + x_{1} x_{2}^{4} + x_{2}^{5}\right) \\

.. math::

    \langle 1 \rangle &= 1 \\
    \langle x \rangle &= \frac{1}{3}( x_0 + x_1 + x_2 ) \\
    \langle x^2 \rangle &= \frac{1}{6}( x_0^2 + x_0 x_1 + x_1^2 + x_0 x_2 + x_1 x_2 + x_2^2 ) \\
    \langle x^3 \rangle &= \frac{1}{10}( x_0^3 + x_0^2 x_1 + x_0 x_1^2 + x_1^3 + x_0^2 x_2 + x_0 x_1 x_2 + x_1^2 x_2 + x_0 x_2^2 + x_1 x_2^2 + x_2^3 ) \\
    \langle x y \rangle  &= \frac{1}{12}( 2 x_0 y_0 + x_1 y_0 + x_2 y_0 + x_0 y_1 + 2 x_1 y_1 + x_2 y_1 + x_0 y_2 + x_1 y_2 + 2 x_2 y_2 ) \\

Those match exactly with the results from our algorithm.

.. autofunction:: quadrature.test_triangle

Line Segment :math:`\triangle^1` in :math:`\mathbb{R}^2`
---------------------------------------------------------------

The same be done for line segments defined by their endpoints
:math:`(x_0,y_0)` and :math:`(x_1,y_1)`. Again the first few moments in this
case are:

.. math::

    \langle 1  \rangle 	&= 1 \\
    \langle x  \rangle 	&= \frac{1}{2}(x_0 + x_1) \\
    \langle x^2 \rangle 	&= \frac{1}{3}(x_0^2 + x_0 x1 + x_1^2) \\
    \langle x^3 \rangle 	&= \frac{1}{4}(x_0^3 + x_0^2 x_1 + x_0 x_1^2 + x_1^3) \\
    \langle x y \rangle   &= \frac{1}{6}(2 x_0 y_0 + x_1 y_0 + x_0 y_1 + 2 x_1 y_1) \\
    \langle x^2 y \rangle &= \frac{1}{12}(3 x_0^2 y_0 + 2 x_0 x_1 y_0 + x_1^2 y_0 + x_0^2 y_1 + 2 x_0 x_1 y_1 + 3 x_1^2 y_1) \\

And our algorithm gives

.. autofunction:: quadrature.test_line_segment

Tetrahedron :math:`\triangle^3` in :math:`\mathbb{R}^3`
--------------------------------------------------------

.. math::
    \langle 1 \rangle 	&= 1 \\
    \langle x \rangle  &= \frac{1}{4}(x_0 + x_1 + x_2 + x_3) \\
    \langle x^2 \rangle &= \frac{1}{20}(2 x_0^2 + 2 x_0 x_1 + 2 x_1^2 + 2 x_0 x_2 + 2 x_1 x_2 + 2 x_2^2 + 2 x_0 x_3 + 2 x_1 x_3 + 2 x_2 x_3 + 2 x_3^2) \\
    \langle x^3 \rangle &= \frac{1}{20}(x_0^3 + x_0^2 x_1 + x_0 x_1^2 + x_1^3 + x_0^2 x_2 + x_0 x_1 x_2 + x_1^2 x_2 + x_0 x_2^2 + x_1 x_2^2 + x_2^3 + x_0^2 x_3 + x_0 x_1 x_3 + x_1^2 x_3 + x_0 x_2 x_3 + \\
                   &  + x_1 x_2 x_3 + x_2^2 x_3 + x_0 x_3^2 + x_1 x_3^2 + x_2 x_3^2 + x_3^3) \\
    \langle xy \rangle  &= \frac{1}{20}(2 x_0 y_0 + x_1 y_0 + x_2 y_0 + x_3 y_0 + x_0 y_1 + 2 x_1 y_1 + x_2 y_1 + x_3 y_1 + x_0 y_2 + x_1 y_2 + 2 x_2 y_2 + x_3 y_2 + x_0 y_3 + x_1 y_3 + x_2 y_3 + 2 x_3 y_3) \\

.. autofunction:: quadrature.test_tetrahedron

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

