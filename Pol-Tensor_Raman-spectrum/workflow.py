r"""
File defining the abstract base class for a BigDFT workflow, enabling to
do a given series of computation and post-processing, whatever the
calculator.
"""

from __future__ import print_function
from abc import ABCMeta, abstractmethod
import yaml


# The three space coordinates
COORDS = ["x", "y", "z"]
# A dictionary for the two possible signs
SIGNS = {"+": 1., "-": -1.}


class Workflow(object):
    r"""
    This is the abstract base class used to define all other workflow
    classes.
    """
    __metaclass__ = ABCMeta

    def __init__(self, calculator, input_base):
        r"""
        Any workflow requires a calculator and a set of input parameters
        to be initialized.

        :param calculator: A BigDFT Calculator.
        :type calculator: GIBinding or SystemCalculator
        :param input_base: Input parameters.
        :type input_base: str
        """
        # For the moment, the input must be a string that contains the
        # initial positions.
        if "posinp" not in input_base:
            raise ValueError("The initial positions must be included in the "
                             "input file.")
        self.input_base = yaml.load(input_base)
        self.calculator = calculator

    @abstractmethod
    def run(self):
        r"""
        Method running all the desired calculations.
        """
        pass

    @abstractmethod
    def _post_processing(self):
        r"""
        Method running the post-processing of the calculation(s).
        """
        pass
