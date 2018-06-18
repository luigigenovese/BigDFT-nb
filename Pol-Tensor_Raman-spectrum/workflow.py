r"""
File defining the abstract base class for a BigDFT workflow, enabling to
do a given series of computation and post-processing, whatever the
calculator.
"""

from __future__ import print_function
from abc import ABCMeta, abstractmethod
import yaml


COORDS = ["x", "y", "z"]
SIGNS = {"+": 1., "-": -1.}


class Workflow(object):
    __metaclass__ = ABCMeta

    def __init__(self, calculator, input_base):
        if "posinp" not in input_base:
            raise ValueError("The initial positions must be included in the "
                             "input file.")
        self.input_base = yaml.load(input_base)
        self.calculator = calculator

    @abstractmethod
    def _post_processing(self):
        pass
