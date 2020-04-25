"""
Plotters
========

#. :class:`.ProgressPlotter`
#. :class:`.SelectionPlotter`


The are multiple different plotters, but they are not related by an
inheritance hierarchy. Each plotter has its own distinct API, according
to its needs.

"""

from .progress import *  # noqa
from .selection import *  # noqa
