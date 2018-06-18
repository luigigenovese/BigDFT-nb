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
        self.input_base = input_base
        self._calculator = calculator

    @property
    def input_base(self):
        r"""
        Getter of the base input parameters.
        """
        return self._input_base

    @input_base.setter
    def input_base(self, new_input):
        # For the moment, the input must be a string that contains the
        # initial positions. This is checked by this setter.
        if "posinp" not in new_input:
            raise ValueError("The initial positions must be included in the "
                             "input file.")
        self._input_base = yaml.load(new_input)

    @property
    def calculator(self):
        r"""
        :returns: Calculator used to perform the BigDFT calculations.
        """
        return self._calculator

    @property
    def inputs(self):
        r"""
        :returns: Input parameters for each calculation to be performed.
        :rtype: dict
        """
        return self._inputs

    @property
    def posinp_base(self):
        r"""
        :returns: Positions provided in the base input file.
        :rtype: list
        """
        return self.input_base["posinp"]["positions"]

    def run(self):
        r"""
        Method initializing and running all the calculations and then
        performing the post-processing.
        """
        self._initialize_calculations()
        self._run_calculations()
        self._post_processing()

    @abstractmethod
    def _initialize_calculations(self):
        r"""
        Method initializing the calculation(s) (mainly defining the
        different set of parameters for each calculation).
        """
        self._inputs = {}

    @abstractmethod
    def _run_calculations(self):
        r"""
        Method running the calculation(s).
        """
        pass

    @abstractmethod
    def _post_processing(self):
        r"""
        Method running the post-processing of the calculation(s).
        """
        pass
