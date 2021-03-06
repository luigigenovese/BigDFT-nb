from __future__ import print_function
from copy import deepcopy
import numpy as np
import yaml
from BigDFT.Calculators import GIBinding  # , SystemCalculator
from BigDFT.Logfiles import Logfile
from workflow import Workflow, COORDS, SIGNS


class PolTensorWorkflow(Workflow):
    r"""
    This class allows to run all the calculations allowing to compute
    the polarizability tensor after some post-processing.

    To compute the polarizability tensor, one must compute the variation
    of the dipole moment of the molecule along one direction with
    respect to an infinitesimal variation of the electric field
    amplitude in one direction. To do that by means of second order
    finite difference, six DFT calculations must be performed, that is
    two per space coordinate: one with a positive electric field
    amplitude, the other with a negative one, this being repeated for an
    electric field along the :math:`x`, :math:`y` and :math:`z` axis.
    """

    RADICAL = "EF_along_"

    def __init__(self, calculator, input_base, ef_amplitudes=[1.E-4]*3):
        r"""
        :param calculator: A BigDFT Calculator.
        :type calculator: GIBinding or SystemCalculator
        :param input_base: Input parameters.
        :type input_base: str
        :param ef_amplitudes: Electric field amplitudes to be applied
            along each of the three space coordinates.
        :type ef_amplitudes: list (of length 3)
        """
        super(PolTensorWorkflow, self).__init__(calculator, input_base)
        self.ef_amplitudes = ef_amplitudes

    @property
    def ef_amplitudes(self):
        r"""
        :returns: Amplitudes of the electric field to be applied along
            each space coordinate.
        :rtype: array of length three
        """
        return self._ef_amplitudes

    @ef_amplitudes.setter
    def ef_amplitudes(self, new_ef_amplitudes):
        # This setter makes sure that there are three elements to define
        # the electric field amplitude along each space coordinate
        try:
            assert len(new_ef_amplitudes) == 3
            self._ef_amplitudes = new_ef_amplitudes
        except AssertionError:
            raise ValueError("There must be three electric field amplitudes, "
                             "one for each space coordinate.")

    @property
    def dipoles(self):
        r"""
        :returns: Output dipole of each calculation.
        :rtype: dict
        """
        return self._dipoles

    @property
    def efields(self):
        r"""
        :returns: All the electric fields the system must undergo.
        :rtype: dict
        """
        return self._efields

    def _initialize_calculations(self):
        r"""
        Method initializing the input parameters for each calculation.
        """
        self._set_efields()
        self._initialize_inputs()

    def _set_efields(self):
        r"""
        Function used to define the six electric fields the system must
        undergo from the input amplitudes of electric field.
        """
        self._efields = {}
        for i, coord in enumerate(COORDS):
            amplitude = self.ef_amplitudes[i]
            for sign in SIGNS:
                # Set the electric field vector
                if amplitude is not None and amplitude != 0.:
                    vector = [0.0] * 3
                    vector[i] = SIGNS[sign] * amplitude
                else:
                    vector = None
                # Add a new element to the dictionary of electric fields
                key = self.RADICAL + coord + sign
                self._efields[key] = vector

    def _initialize_inputs(self):
        r"""
        Initialize the six input files required to compute the
        polarizability tensor. They are stored in a dictionary whose
        keys allow to distinguish each calculation.
        """
        self._inputs = {}
        for key, efield in self.efields.iteritems():
            if efield is None:
                self._inputs[key] = None
            else:
                new_input = deepcopy(self.input_base)
                new_input['dft']['elecfield'] = efield
                self._inputs[key] = new_input

    def _run_calculations(self):
        r"""
        Method running all the calculations. 
        
        It also gathers the output dipole of the system for each
        calculation in the 'dipoles' attribute as a dictionary.
        """
        # The dipoles attribute will contain the output dipole of each
        # calculation. It will be a dictionary with the same keys as the
        # inputs dictionary.
        self._dipoles = {}
        # Loop over the inputs to run the calculations
        for key, inp in self.inputs.iteritems():
            # Unwrap coordinate and directions from the key
            coord = key[-2]
            direction = key[-1]
            # Set the message for verbosity
            if SIGNS[direction] == 1.:
                msg = "positive "
            else:
                msg = "negative "
            msg += "electric field applied along {}!".format(coord)
            # Do action based on the existence of an input file
            if inp is None:
                # No calculation if there is no input file
                msg = "no " + msg
                self._dipoles[key] = np.zeros(3)
            else:
                # Run the calculation if there is an input file
                # TODO: The difference between the calculators should be
                #       minimal: use an abstract base class.
                if isinstance(self.calculator, GIBinding):
                    inp['logfile'] = True
                    inp['radical'] = key
                    self.calculator.set(inp)
                    self.calculator.run()
                else:
                    with open(key+".yaml", "w") as f:
                        yaml.dump(inp, f)
                    self.calculator.run(name=key)
                # Treat the logfile to retrieve information
                self._dipoles[key] = self._get_dipole(key)
            print(msg.capitalize())

    def _get_dipole(self, key):
        r"""
        :param key: Key leading to a particular logfile.
        :type key: str
        :returns: Dipole read from the above-mentioned logfile.
        :rtype: numpy array of length 3
        """
        log = Logfile("log-{}.yaml".format(key))
        dipole = log.log['Electric Dipole Moment (AU)']['P vector']
        return np.array(dipole)

    def _post_processing(self):
        r"""
        Method running the post-processing of the calculations, here:
        * computing the polarizability tensor
        """
        self.pol_tensor = self._build_polarizability_tensor()

    def _build_polarizability_tensor(self):
        r"""
        Function returning the polarizability tensor. It corresponds
        to the response of the system (here, the modification of its
        dipole) when an electric field is applied.

        The dipole and the electric field being vectors, the
        polarizability is represented by a tensor. Its elements
        :math:`alpha_{i, j} = d D_i / d E_j` represent the
        proportionality coefficient between the dipole :math:`D` in the
        direction :math:`i` when an electric field of amplitude
        :math:`E` is applied in the direction :math:`j`. (:math:`i` and
        :math:`j` represent one of the :math:`x`, :math:`y` or :math:`z`
        axis).

        :returns: Polarizability tensor
        :rtype: 2D numpy array of dimension 3*3
        """
        pol_tensor = []
        # Loop over the electric field directions
        for i, coord in enumerate(COORDS):
            ef_amplitude = self.ef_amplitudes[i]
            if ef_amplitude is None or ef_amplitude == 0.:
                new_line = np.zeros(3)
            else:
                # Get the delta of the dipoles when the electric field
                # is applied along a given coordinate
                dipole_p = self.dipoles[self.RADICAL+coord+"+"]
                dipole_m = self.dipoles[self.RADICAL+coord+"-"]
                delta_dipoles = dipole_p - dipole_m
                # Find the delta of the electric field (which is twice
                # the absolute value of the electric field amplitude
                # along the direction considered)
                delta_efield = 2 * abs(ef_amplitude)
                # Append the new line of the polarizability tensor,
                # defined as the the ratio of the delta of the dipoles
                # and the delta of the electric fields
                new_line = delta_dipoles / delta_efield
            pol_tensor.append(new_line)
        # Return the polarizability tensor as a numpy array
        return np.array(pol_tensor)
