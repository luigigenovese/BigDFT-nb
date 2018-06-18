from __future__ import print_function
from copy import deepcopy
import numpy as np
import yaml
from BigDFT.Calculators import GIBinding  # , SystemCalculator
from BigDFT.Logfiles import Logfile
from workflow import Workflow, COORDS, SIGNS


# Mass of the different types of atoms in atomic mass units
# TODO: Add more types of atoms
#       (found in $SRC_DIR/bigdft/src/orbitals/eleconf-inc.f90)
MASS_ATOMS = {"H": 1.00794, "He": 4.002602, "Li": 6.941, "Be": 9.012182,
              "B": 10.811, "C": 12.011, "N": 14.00674, "O": 15.9994,
              "F": 18.9984032, "Ne": 20.1797}
# Conversion from atomic to electronic mass unit
AMU_TO_EMU = 1.660538782e-27 / 9.10938215e-31
# Conversion factor from Hartree to cm^-1
HA_TO_CMM1 = 219474.6313705


class RamanWorkflow(Workflow):
    r"""
    This class allows to run all the calculations allowing to compute
    the energies of the Raman spectrum after some post-processing.

    To get the energies of the Raman spectrum, one needs to find the
    eigenvalues of the dynamical matrix, that is closely related to the
    Hessian. To build these matrices, one must find the derivatives of
    the forces when each coordinate of each atom is translated by a
    small amount around the equilibrium positions. To get a better
    precision on the derivative, each coodinate is translated
    positively and negatively, so that the number of BigDFT calculations
    amounts to :math:`2*3*n_{at} = 6 n_{at}`, where :math:`n_{at}` is
    the number of atoms (3 for the coordinates (:math:`x`, :math:`y` and
    :math:`z`), 2 for the number of calculations per coordinates).
    """

    def RADICAL(self, i_at):
        r"""
        :param i_at: Index of the atom that is being translated.
        :type i_at: int
        :returns: Radical for the name of the calculation.
        :rtype: str
        """
        return "atom{:04d}_".format(i_at)

    def __init__(self, calculator, input_base, amplitudes=[0.45/64]*3):
        r"""
        :param calculator: A BigDFT Calculator.
        :type calculator: GIBinding or SystemCalculator
        :param input_base: Input parameters.
        :type input_base: str
        :param ef_amplitudes: Amplitudes of the translations to be
            applied along each of the three space coordinates.
        :type ef_amplitudes: list (of length 3)
        """
        super(RamanWorkflow, self).__init__(calculator, input_base)
        self.amplitudes = amplitudes

    @property
    def amplitudes(self):
        r"""
        :returns: Amplitudes of the displacements to be applied along
            each space coordinate.
        :rtype: array of length three
        """
        return self._amplitudes

    @amplitudes.setter
    def amplitudes(self, new_amplitudes):
        # This setter makes sure that there are three elements to define
        # the displacement amplitude along each space coordinate
        try:
            assert len(new_amplitudes) == 3
            self._amplitudes = new_amplitudes
        except AssertionError:
            raise ValueError("There must be three translation amplitudes, "
                             "one for each space coordinate.")

    @property
    def forces(self):
        r"""
        :returns: Output forces of each calculation.
        :rtype: dict
        """
        return self._forces

    @property
    def displacements(self):
        r"""
        :returns: Displacements each atom of the system must undergo.
            There are six of them (two per space coordinate) in order
            to be able to compute the Hessian matrix by using the
            central difference scheme.
        :rtype: dict of length 6
        """
        return self._displacements

    def _initialize_calculations(self):
        r"""
        Method initializing the input parameters for each calculation.
        """
        self._set_displacements()
        self._initialize_inputs()

    def _set_displacements(self):
        r"""
        Function used to define the six similar displacements each atom
        must undergo from the amplitudes of displacement in each
        direction.
        """
        self._displacements = {}
        for i, coord in enumerate(COORDS):
            amplitude = self.amplitudes[i]
            for sign in SIGNS:
                # Set the displacement vector
                if amplitude is not None and amplitude != 0.:
                    vector = [0.0] * 3
                    vector[i] = SIGNS[sign] * amplitude
                else:
                    vector = None
                # Add a new element to the dictionary of displacements
                key = coord + sign
                self._displacements[key] = vector

    def _initialize_inputs(self):
        r"""
        Initialize the six input files per atom required to compute the
        Raman spectrum. They are stored in a dictionary whose
        keys allow to distinguish each calculation.
        """
        self._inputs = {}
        posinp = self.posinp_base
        for i_at in range(len(posinp)):
            # Loop over the space coordinates along which the atom has
            # to be displaced
            for disp_key, vector in self.displacements.iteritems():
                key = self.RADICAL(i_at) + disp_key
                # Do not create an input if there is no displacement
                if vector is None:
                    self._inputs[key] = None
                # Else, create a new input file
                else:
                    self._inputs[key] = self._translate_atom(i_at, vector)

    def _translate_atom(self, i_at, vector):
        r"""
        Method returning a new input file where the atom i_at has been
        translated according the vector.

        .. warning::

            You have to make sure that the units of the vector match
            those used by the posinp.

        :param i_at: Index of the atom.
        :type i_at: int
        :param vector: Translation vector to apply.
        :type vector: list or numpy.array of length 3
        :returns: A new input file with the tranlsated atom.
        :rtype: dict
        """
        new_input = deepcopy(self.input_base)
        atom_to_translate = new_input["posinp"]["positions"][i_at]
        at_type = atom_to_translate.keys()[0]
        for i, dpos in enumerate(vector):
            atom_to_translate[at_type][i] += dpos
        return new_input

    def _run_calculations(self):
        r"""
        Method running all the calculations.

        It also gathers the output forces of the system for each
        calculation in the 'forces' attribute as a dictionary.
        """
        # The forces attribute will contain the output forces of each
        # calculation. It will be a dictionary with the same keys as the
        # inputs dictionary.
        self._forces = {}
        # Loop over the inputs to run the calculations
        for key, inp in self.inputs.iteritems():
            # Unwrap atom index, coordinate and direction from the key
            i_at = int(key[4:8])
            coord = key[-2]
            direction = key[-1]
            # Set the message for verbosity
            if SIGNS[direction] == 1.:
                msg = "positive "
            else:
                msg = "negative "
            msg += "translation of atom {} along {}!".format(i_at, coord)
            # Do action based on the existence of an input file
            if inp is None:
                # No calculation if there is no input file
                msg = "no " + msg
                # All the forces are set to 0
                self._forces[key] = [{"X": [0.0]*3}
                                     for _ in range(len(self.posinp_base))]
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
                # Store the forces
                log = Logfile("log-{}.yaml".format(key))
                self._forces[key] = log.forces
            print(msg.capitalize())

    def _post_processing(self):
        r"""
        Method running the post-processing of the calculations, here:
        * computing the dynamical matrix,
        * solving it to get the normal modes and their energies.
        """
        # - Set the dynamical matrix
        self.dyn_mat = self._build_dyn_mat()
        # - Find the energies as eigenvalues of the dynamical matrix
        self.energies = {}
        self.energies['Ha'], self.normal_modes = self._solve_dynamical_matrix()
        self.energies['cm^-1'] = self.energies['Ha'] * HA_TO_CMM1

    def _solve_dynamical_matrix(self):
        r"""
        Method solving the dynamical matrix to get the phonon energies
        (converted in Hartree) and the eigenvectors.

        :returns: Tuple made of the eigenvalues (as an array) and the
            eigenvectors (as a matrix).
        :rtype: tuple
        """
        eigs, vecs = np.linalg.eig(self.dyn_mat)
        eigs = np.array([np.sqrt(-e) if e < 0 else np.sqrt(e) for e in eigs])
        return eigs, vecs

    def _build_dyn_mat(self):
        r"""
        Method computing the dynamical matrix of the system. It is
        very similar to the Hessian matrix: its elements are only
        corrected by a weight :math:`w,` which is the inverse of the
        square-root of the product of the atomic masses of the atoms
        involved in the Hessian matrix element.

        The masses are counted in electronic mass units (which is the
        atomic unit of mass, that is different from the atomic mass
        unit).

        :returns: Dynamical matrix.
        :rtype: 2D square numpy array of dimension :math:`3 n_{at}`
        """
        # Numpy does the ratio of arrays intelligently: by making
        # masses an array of the same size as the Hessian, there is
        # nothing but the ratio of both arrays to perform to get
        # the dynamical matrix.
        h = self._build_hessian()
        masses = self._build_masses()
        return h / masses

    def _build_masses(self):
        r"""
        Method computing the masses array used to define the dynamical
        matrix. The masses are counted in electronic mass units (which
        is the atomic unit of mass, that is different from the atomic
        mass unit).

        :returns: Masses matrix.
        :rtype: 2D square numpy array of dimension :math:`3 n_{at}`
        """
        # Get the atoms of the system from the reference posinp
        posinp = self.posinp_base
        atom_types = [atom.keys()[0] for atom in posinp]
        # Build the masses matrix (the loops over range(3) are here
        # to ensure that masses has the same dimension as the Hessian)
        masses = [[MASS_ATOMS[atom1] * MASS_ATOMS[atom2]
                   for atom2 in atom_types for j in range(3)]
                  for atom1 in atom_types for i in range(3)]
        # Return the masses as a numpy array, converted in electronic
        # mass units
        return np.sqrt(masses) * AMU_TO_EMU

    def _build_hessian(self):
        r"""
        Method computing the Hessian of the system. Its dimension is
        :math:`3 n_{at}`, where :math:`n_{at}` is the number of atoms of
        the system.

        :returns: Hessian matrix.
        :rtype: 2D square numpy array of dimension :math:`3 n_{at}`
        """
        # Initialization of variables
        h = []  # Hessian matrix
        n_at = len(self.posinp_base)  # Number of atoms
        # First loop over all atoms
        for i_at in range(n_at):
            # Loop over the coordinates (x, y and z)
            for i_coord, coord in enumerate(COORDS):
                # The Hessian is made of the delta of the forces
                # with respect to the delta of the move distances.
                new_line = []
                # Make sure there is no div. by 0
                amplitude = self.amplitudes[i_coord]
                if amplitude is None or amplitude == 0.0:
                    amplitude = 1.0
                # The Hessian is then built line by line:
                # 1- Find the delta displacement. It is twice the
                #    distance of the positive move along the direction
                #    of the displacement.
                delta_x = 2 * amplitude
                # 2- Find the delta forces for each atom and update
                #    the new line of the Hessian.
                key_base = self.RADICAL(i_at) + coord
                keys = [key_base+sign for sign in SIGNS]
                for j_at in range(n_at):
                    forces = [np.array(self._forces[key][j_at].values()[0])
                              for key in keys]
                    delta_forces = forces[0] - forces[1]
                    new_line += list(delta_forces/delta_x)
                # The new line of the Hessian is now complete
                h.append(new_line)
        # Return the Hessian matrix as a numpy array
        return np.array(h)
