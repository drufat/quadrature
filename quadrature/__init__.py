"""
The functions below allow us to integrate arbitrary polynomials
on arbitrary simplices.
TODO: Add more extensive documentation, and describe the theory. See:
http://math.stackexchange.com/questions/3990/computing-the-moments-of-a-triangle
"""

from .integral import *

try:
    from .symbolic import *
except:
    pass