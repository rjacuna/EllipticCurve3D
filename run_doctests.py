#!/usr/bin/env sage -python
"""Run the doctests of elliptic_curve_3d.sage (this Sage build lacks the sage-runtests script):
    sage -python run_doctests.py
"""
import sys, os
from sage.all import *
from sage.doctest.control import DocTestDefaults, DocTestController
here = os.path.dirname(os.path.abspath(__file__))
DC = DocTestController(DocTestDefaults(), [os.path.join(here, "elliptic_curve_3d.sage")])
sys.exit(DC.run())
