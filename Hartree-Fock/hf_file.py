from numpy import array


def xyz_reader(file_name):
    # Reads an xyz file and returns the number of atoms, atom types and
    # atom coordinates.

    file = open(file_name, 'r')
    print(type(file))

    number_of_atoms = 0
    atom_type = []
    atom_coordinates = []

    for idx, line in enumerate(file):
        # Get number of atoms
        if idx == 0:
            try:
                number_of_atoms = line.split()[0]
            except LookupError:
                print(""""xyz file not in correct format. Make sure the format follows:
                          https://en.wikipedia.org/wiki/XYZ_file_format"""
                      )

        # Skip the comment/blank line
        if idx == 1:
            continue

        # Get atom types and positions
        if idx != 0:
            split = line.split()
            atom = split[0]
            coordinates = array([float(split[1]),
                                 float(split[2]),
                                 float(split[3])])

            atom_type.append(atom)
            atom_coordinates.append(coordinates)

    file.close()

    return int(number_of_atoms), atom_type, atom_coordinates

