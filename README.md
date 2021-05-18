# extract_trajectory_parameters.exe:
Extract geometry parameters from a multi-frame xyz format trajectory file.

## Usage: 
extract_trajectory_parameters.exe [-h | --help]
extract_trajectory_parameters.exe [-t | --traj TRAJ] [-i | --index INDEX]

TRAJ: name of a xyz format trajectory file with multi frames.
INDEX: name of an index file.

The index file should contains several lines, each line contains 2 to 4 integers.
Those integers are the indices of the atoms, if 2 integers are provided, the 
bond length of these two atoms will be calcualted, with the unit of input, 
3 or 4 integers stands for angle (degree) and dihedral angle (degree).

All file names that are not provided in the command argument, 
will be asked interactively.
If INDEX is "-", the program will read indices from STDIN.
In this mode, you need to input the amount of parameters ro be measured first, 
when the screen hints you to do so.

