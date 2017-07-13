"""
File containing the Input and Posinp classes.
"""

from __future__ import print_function
import yaml
import os
from copy import deepcopy

class Input(dict):
    """
    Class allowing to initialize, read, write and interact with the 
    input file of a BigDFT calculation.
    """

    def __init__(self, data):
        """
        An input file is created from a yaml dictionary.
        
        :param data: yaml dictionary of the input file.
        :type data: dict
        """
        for key, value in data.items():
            self[key] = value
        #if isinstance(data, dict):
        #    self.data = data
        #else:
        #    raise TypeError("The argument must be a (or inherit from) dict.")


    @classmethod
    def from_file(cls, filename):
        """
        Initializes the Input instance from a file

        :param filename: Name of the file to read.
        :type filename: str
        :returns: Input instance initialized from the file.
        :rtype: Input
        """
        with open(filename, "r") as f:
            data = yaml.load(f)
        return cls(data)


    @classmethod
    def from_string(cls, string):
        """
        Initializes the Input instance from a file

        :param string: String of the input file.
        :type string: str
        :returns: Input instance initialized from the string.
        :rtype: Input
        """
        data = yaml.load(string)
        return cls(data)


    def write(self, filename):
        """
        Method writing an Input instance in filename.

        :param filename: Name of the file to read to initialize the
                         Input instance.
        :type filename: str
        """

        with open(filename, "w") as f:
            for key in self:
                yaml.dump({key: self[key]}, f)



class Posinp:
    """
    Class allowing to initialize, read, write and interact with the 
    input geometry of a BigDFT calculation.
    """

    def __init__(self, filename):
        """
        A file is read to initialize a Posinp instance.
        
        :param filename: Name of the file to read to initialize the
                         Posinp instance.
        :type filename: str
        """
        with open(filename, "r") as f:
            # Read the first line, containing the number of atoms and 
            # the units of the coordinates of each atom
            content = f.readline().split()
            self.n_at = int(content[0])
            self.units = content[1]

            # Read the second line, containing the boundary conditions
            content = f.readline().split()
            self.BC = content[0]

            # Loop over all the atoms to read the positions.
            self.atoms = []
            for i in range(self.n_at):
                content = f.readline().split()
                self.atoms.append({"Type": content[0], \
                    "Position": [float(c) for c in content[1:4]]})


    @classmethod
    def from_string(cls, string):
        """
        Class method allowing to initialize a posinp file from a string.

        :param string: String of the posinp file.
        :type string: str
        :returns: Posinp instance initialized from the string.
        :rtype: Posinp
        """
        # Write the string in a temporary file
        tmpfile = "tmp-posinp.xyz"
        with open(tmpfile, "w") as f:
            f.write(string)
        # Initialize the posinp from the file
        posinp = cls(tmpfile)
        # Remove the temporary file before returning the posinp
        os.remove(tmpfile)
        return posinp


    def __str__(self):
        """
        Method converting the Posinp instance to a string.

        :returns: The Posinp instance as a string.
        :rtype: str
        """
        # Create the first two lines of the posinp file
        pos_str = "{}   {}\n".format(self.n_at, self.units)
        pos_str += "{}\n".format(self.BC)
        # Loop over the atoms to create a line containing their type 
        # and coordinates
        for atom in self.atoms:
            pos_str += "{t}   {pos[0]: .15}   {pos[1]: .15}   {pos[2]: .15}\n"\
            .format(t=atom['Type'], pos=atom['Position'])
        return pos_str


    def __eq__(self, other):
        """
        Two Posinp instances are the same if their string representation 
        atoms at the same positions, and use the same units and 
        boundary conditions.

        :returns: A boolean telling if self and other are the same 
                  Posinp instances.
        :rtype: bool
        """
        return str(self) == str(other)


    def write(self, filename):
        """
        Method writing the posinp in filename.

        :param filename: Name of the file to read to initialize the
                         Posinp instance.
        :type filename: str
        """

        with open(filename, "w") as f:
            f.write(str(self))


    def translate_atom(self, i_at, vector):
        """
        Method returning a new posinp where the atom i_at has been 
        translated according the vector. 

        ..Warning: You have to make sure that the units of the vector 
                   match those used by the posinp.
        
        :param i_at: Index of the atom.
        :type i_at: int
        :param vector: Translation vector to apply.
        :type vector: list or numpy.array of length 3
        :returns: A new posinp with the tranlsated atom
        :rtype: posinp.Posinp
        """
        # Create a deep copy of the posinp and return it after the 
        # update of the position of the atom i_at.
        new_posinp = deepcopy(self)
        for i in range(3):
            new_posinp.atoms[i_at]['Position'][i] += vector[i]

        return new_posinp
