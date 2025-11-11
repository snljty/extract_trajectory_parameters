/*************************************************************************************
 * This program extracts geometry parameters from a multi-frame xyz trajectory file. *
 * Note: the amount of atoms and their types are assumed to be unique in any frames. * 
 *************************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include <float.h>
# include <stdbool.h>

# define num_coords 3u /* x y z */

double Radian_to_degree(double r);
void Vector_minus(double const *a, double const *b, double *res);
double Get_vector_dot_product(double const *a, double const *b);
double Get_vector_length(double const *a);
void Normalize_vector(double *a);
double Get_distance(double const *a, double const *b);
double Get_vector_angle(double const *a, double const *b);
double Get_angle(double const *a, double const *b, double const *c);
void Vector_cross_product(double const *a, double const *b, double *res);
double Get_dihedral(double const *a, double const *b, double const *c, double const *d);
void Output_parameters(FILE *out_file_ptr, unsigned int num_atoms, double const *const *atom_coords, \
                       unsigned int num_paras, unsigned int const *const *para_index);

int main(int argc, char const **argv)
{
    int iarg = 0;
    char traj_name[BUFSIZ + 1] = {'\0'};
    char index_name[BUFSIZ + 1] = {'\0'};
    char *traj_name_use = traj_name;
    char *index_name_use = index_name;
    char output_name[BUFSIZ + 1] = {'\0'};
    FILE *traj_fp = NULL;
    FILE *index_fp = NULL;
    FILE *output_fp = stdout;
    unsigned int num_paras = 0u;
    char buf[BUFSIZ + 1] = "";
    char *tok = NULL;
    bool read_indices_from_stdin = false;
    unsigned int const max_num_index = 4u;
    unsigned int **para_index = NULL;
    unsigned int i_para = 0u;
    unsigned int num_atoms = 0u;
    unsigned int i_index = 0u;
    unsigned int i_atom = 0u;
    char **atom_name = NULL;
    double **atom_coords = NULL;
    unsigned int const max_atom_name_space = strlen("Bq") + 1u;
    unsigned int i_coord = 0u;
    unsigned int i_frame = 0u;
    double v_double = 0.0;
    unsigned int v_uint = 0u;
    char v_str[24] = "";

    /* Get file names and open files */
    for (iarg = 1; iarg < argc; ++ iarg)
    {
        if (! strcmp(argv[iarg], "-h") || ! strcmp(argv[iarg], "--help") || ! strcmp(argv[iarg], "/?"))
        {
            printf("Usage: %s [-h | --help]\n", argv[0]);
            printf("       %s [-t | --traj TRAJ] [-i | --index INDEX] [-o | --output OUTPUT]\n", argv[0]);
            printf("%s\n", "");
            printf("%s\n", "TRAJ: name of a xyz format trajectory file with multi frames.");
            printf("%s\n", "INDEX: name of an index file.");
            printf("%s\n", "OUTPUT: name of file for output.");
            printf("%s\n", "");
            printf("%s\n", "The index file should contains several lines, each line contains 2 to 4 integers.");
            printf("%s\n", "Those integers are the indices of the atoms, if 2 integers are provided, the ");
            printf("%s\n", "bond length of these two atoms will be calcualted, with the unit of input, ");
            printf("%s\n", "3 or 4 integers stands for angle (degree) and dihedral angle (degree).");
            printf("%s\n", "");
            printf("%s\n", "If either TRAJ or INDEX is omitted in the command arguments, ");
            printf("%s\n", "it will be asked interactively.");
            printf("%s\n", "If INDEX is \"-\", the program will read indices from STDIN.");
            printf("%s\n", "In this mode, you need to input the amount of parameters to be measured first, ");
            printf("%s\n", "when the screen hints you to do so.");
            printf("%s\n", "If OUTPUT is omitted in the command arguments or it is \"-\", stdout will be used.");
            printf("%s\n", "");
            printf("%s\n", "Exiting normally.");
            exit(EXIT_SUCCESS);
        }
    }
    iarg = 1;
    while (iarg < argc)
    {
        if (! strcmp(argv[iarg], "-t") || ! strcmp(argv[iarg], "--traj"))
        {
            ++ iarg;
            if (iarg >= argc)
            {
                fprintf(stderr, "Error! missing argument for \"%s\".\n", argv[iarg - 1]);
                exit(EXIT_FAILURE);
            }
            strncpy(traj_name, argv[iarg], BUFSIZ + 1);
        }
        else if (! strcmp(argv[iarg], "-i") || ! strcmp(argv[iarg], "--index"))
        {
            ++ iarg;
            if (iarg >= argc)
            {
                fprintf(stderr, "Error! missing argument for \"%s\".\n", argv[iarg - 1]);
                exit(EXIT_FAILURE);
            }
            strncpy(index_name, argv[iarg], BUFSIZ + 1);
        }
        else if (! strcmp(argv[iarg], "-o") || ! strcmp(argv[iarg], "--output"))
        {
            ++ iarg;
            if (iarg >= argc)
            {
                fprintf(stderr, "Error! missing argument for \"%s\".\n", argv[iarg - 1]);
                exit(EXIT_FAILURE);
            }
            strncpy(output_name, argv[iarg], BUFSIZ + 1);
        }
        else
        {
            fprintf(stderr, "Error! Unrecognizable argument: \"%s\"\n", argv[iarg]);
            exit(EXIT_FAILURE);
        }
        ++ iarg;
    }
    if (! * traj_name)
    {
        printf("%s\n", "# Input name of the trajectory file:");
        if (! fgets(traj_name, BUFSIZ, stdin))
        {
            exit(EXIT_FAILURE);
        }
        traj_name[strlen(traj_name) - 1] = '\0';
        if (traj_name[0] == '\"')
        {
            traj_name[strlen(traj_name) - 1] = '\0';
            ++ traj_name_use;
        }
    }
    if (strlen(traj_name_use) <= strlen(".xyz") || \
        strcmp(traj_name_use + strlen(traj_name_use) - strlen(".xyz"), ".xyz"))
    {
        fprintf(stderr, "Error! The suffix of a xyz format trajectory file should be \".xyz\".\n");
        exit(EXIT_FAILURE);
    }
    traj_fp = fopen(traj_name_use, "rt");
    if (! traj_fp)
    {
        fprintf(stderr, "Error!\n");
        perror(traj_name_use);
        exit(EXIT_FAILURE);
    }
    if (! * index_name)
    {
        printf("%s\n", "# Input name of the index file:");
        printf("%s\n", "# If you press <Enter> directly, will read indices from STDIN.");
        if (! fgets(index_name, BUFSIZ, stdin))
        {
            exit(EXIT_FAILURE);
        }
        index_name[strlen(index_name) - 1] = '\0';
        if (index_name[0] == '\"')
        {
            index_name[strlen(index_name) - 1] = '\0';
            ++ index_name_use;
        }
    }
    if (! * index_name_use || ! strncmp(index_name_use, "-", strlen("-") + 1))
    {
        read_indices_from_stdin = true;
        printf("%s\n", "# Reading indices from STDIN ...");
        index_fp = stdin;
    }
    else
    {
        index_fp = fopen(index_name_use, "rt");
        if (! index_fp)
        {
            fprintf(stderr, "Error!\n");
            perror(index_name_use);
            exit(EXIT_FAILURE);
        }
    }
    if (! * output_name || ! strcmp(output_name, "-"))
    {
        output_fp = stdout;
    }
    else
    {
        output_fp = fopen(output_name, "rt");
        if (output_fp)
        {
            fclose(output_fp);
            output_fp = NULL;
            fprintf(stderr, "Error! File \"%s\" already exists.\n", output_name);
            exit(EXIT_FAILURE);
        }
        output_fp = fopen(output_name, "wt");
    }

    /* get amount of parameters */
    if (read_indices_from_stdin)
    {
        printf("# Input amount of geometry parameters to be measured: ");
        fflush(stdout);
        for (;;)
        {
            if (! fgets(buf, BUFSIZ, stdin))
            {
                exit(EXIT_FAILURE);
            }
            if (sscanf(buf, "%u", & num_paras) == 1)
            {
                break;
            }
        }
        if (num_paras == 1)
        {
            printf("%s\n", "# Input 1 line of indices:");
        }
        else
        {
            printf("# Input %u lines of indices, one line for one geometry parameter:\n", num_paras);
        }
    }
    else
    {
        /* get amount of parameters */
        while (fgets(buf, BUFSIZ, index_fp))
        {
            tok = strtok(buf, " \n,");
            if (tok)
            {
                ++ num_paras;
            }
        }
        fseek(index_fp, 0, SEEK_SET);
    }

    /* allocate memory for indices */
    para_index = (unsigned int **)malloc(num_paras * sizeof(unsigned int *));
    * para_index = (unsigned int *)malloc(num_paras * max_num_index * sizeof(unsigned int));
    memset(* para_index, 0, num_paras * max_num_index * sizeof(unsigned int));
    for (i_para = 0; i_para < num_paras; ++ i_para)
    {
        para_index[i_para] = * para_index + i_para * max_num_index;
    }

    /* read indices */
    i_para = 0u;
    while (i_para < num_paras && fgets(buf, BUFSIZ, index_fp))
    {
        tok = strtok(buf, " \n,");
        if (tok)
        {
            /* first index of atom */
            if (sscanf(tok, "%u", para_index[i_para]) != 1)
            {
                fprintf(stderr, "Error! Non-blank line %u: the 1st index is not an integer.\n", i_para + 1);
                exit(EXIT_FAILURE);
            }
            tok = strtok(NULL, " \n,");
            /* second index of atom */
            if (! tok)
            {
                fprintf(stderr, "Error! Non-blank line %u: at least 2 indices are required.\n", i_para + 1);
                exit(EXIT_FAILURE);
            }
            if (sscanf(tok, "%u", para_index[i_para] + 1) != 1)
            {
                fprintf(stderr, "Error! Non-blank line %u: the 2nd index is not an integer.\n", i_para + 1);
                exit(EXIT_FAILURE);
            }
            tok = strtok(NULL, " \n,");
            /* third index of atom, optional */
            if (! tok)
            {
                ++ i_para;
                continue;
            }
            if (sscanf(tok, "%u", para_index[i_para] + 2) != 1)
            {
                fprintf(stderr, "Error! Non-blank line %u: the 3rd index is not an integer.\n", i_para + 1);
                exit(EXIT_FAILURE);
            }
            tok = strtok(NULL, " \n,");
            /* fourth index of atom. optional */
            if (! tok)
            {
                ++ i_para;
                continue;
            }
            if (sscanf(tok, "%u", para_index[i_para] + 3) != 1)
            {
                fprintf(stderr, "Error! Non-blank line %u: the 4th index is not an integer.\n", i_para + 1);
                exit(EXIT_FAILURE);
            }
            tok = strtok(NULL, " \n,");
            /* end of indices */
            if (tok)
            {
                fprintf(stderr, "Error! Non-blank line %u: extra non-blank character found.\n", i_para + 1);
                exit(EXIT_FAILURE);
            }
            /* next line */
            ++ i_para;
        }
    }

    /* close index file */
    if (! read_indices_from_stdin)
    {
        fclose(index_fp);
        index_fp = NULL;
    }

    /* print title line of output parameters */
    fprintf(output_fp, "# bond lengthes in input unit, angles and dihedrals in unit degree.\n");
    fprintf(output_fp, "# ");
    fflush(output_fp);
    for (i_para = 0; i_para < num_paras; ++ i_para)
    {
        if (i_para)
        {
            fprintf(output_fp, "%4s", "");
            fflush(stdout);
        }
        if (para_index[i_para][3])
        {
            /* dihedral */
            sprintf(v_str, "%u-%u-%u-%u", para_index[i_para][0], para_index[i_para][1], \
                para_index[i_para][2], para_index[i_para][3]);
            v_uint = strlen(v_str);
            fprintf(output_fp, "%*s%*s%*s", (15 - v_uint) / 2 >= 0 ? (15 - v_uint) / 2 : 0, "", \
                v_uint, v_str, \
                (16 - v_uint) / 2 >= 0 ? (16 - v_uint) / 2 : 0, "");
            fflush(output_fp);
        }
        else
        {
            if (para_index[i_para][2])
            {
                /* angle */
                sprintf(v_str, "%u-%u-%u", para_index[i_para][0], \
                    para_index[i_para][1], para_index[i_para][2]);
                v_uint = strlen(v_str);
                fprintf(output_fp, "%*s%*s%*s", (15 - v_uint) / 2 >= 0 ? (15 - v_uint) / 2 : 0, "", \
                    v_uint, v_str, \
                    (16 - v_uint) / 2 >= 0 ? (16 - v_uint) / 2 : 0, "");
                fflush(output_fp);
            }
            else
            {
                /* length */
                sprintf(v_str, "%u-%u", para_index[i_para][0], para_index[i_para][1]);
                v_uint = strlen(v_str);
                fprintf(output_fp, "%*s%*s%*s", (15 - v_uint) / 2 >= 0 ? (15 - v_uint) / 2 : 0, "", \
                    v_uint, v_str, \
                    (16 - v_uint) / 2 >= 0 ? (16 - v_uint) / 2 : 0, "");
                fflush(output_fp);
            }
        }
    }
    fprintf(output_fp, "\n");

    /*  get number of atoms */
    /* the amount of atoms and their types are assumed to be unique in any frames. */
    if (! fgets(buf, BUFSIZ, traj_fp))
    {
        fprintf(stderr, "Error! Cannot get the first line of trajectory file.\n");
        exit(EXIT_FAILURE);
    }
    if (sscanf(buf, "%u", & num_atoms) != 1)
    {
        fprintf(stderr, "Error! Cannot read the amount of atoms from the first line of trajectory file.\n");
        exit(EXIT_FAILURE);
    }
    if (! fgets(buf, BUFSIZ, traj_fp))
    {
        fprintf(stderr, "Error! Cannot read the title line in frame 1.\n");
        exit(EXIT_FAILURE);
    }

    /* chech whether indices in index file are out of range. */
    for (i_para = 0; i_para < num_paras; ++ i_para)
    {
        for (i_index = 0; i_index < max_num_index; ++ i_index)
        {
            if (para_index[i_para][i_index] > num_atoms)
            {
                fprintf(stderr, "Error! Index %u of non-blank line %u in index file is larger than ", \
                    i_index + 1, i_para + 1);
                fprintf(stderr, "total amount of atoms.\n");
                exit(EXIT_FAILURE);
            }
        }

    }

    /* allocate memory for atom names and atom coordinates */
    atom_name = (char **)malloc(num_atoms * sizeof(char *));
    * atom_name = (char *)malloc(num_atoms * max_atom_name_space * sizeof(char));
    atom_coords = (double **)malloc(num_atoms * sizeof(double *));
    * atom_coords = (double *)malloc(num_atoms * num_coords * sizeof(double));
    for (i_atom = 0; i_atom < num_atoms; ++ i_atom)
    {
        atom_name[i_atom] = * atom_name + i_atom * max_atom_name_space;
        atom_coords[i_atom] = * atom_coords + i_atom * num_coords;
    }

    /*  get atom names  */
    for (i_atom = 0; i_atom < num_atoms; ++ i_atom)
    {
        if (! fgets(buf, BUFSIZ, traj_fp))
        {
            fprintf(stderr, "Error! Cannot read atom %u in frame 1.\n", i_atom + 1);
            exit(EXIT_FAILURE);
        }
        tok = strtok(buf, " ");
        if (! tok)
        {
            fprintf(stderr, "Error! Cannot get the atom name of atom %u, frame 1.\n", i_atom + 1);
            exit(EXIT_FAILURE);
        }
        if (strlen(tok) >= max_atom_name_space)
        {
            printf("# Warning: the name length of atom %u in frame 1 is larger than %u.\n", \
                i_atom + 1, max_atom_name_space - 1);
            printf("# The first %u characters will be used.\n", max_atom_name_space - 1);
        }        
        if (tok[0] > 'Z' || tok[0] < 'A')
        {
            printf("# Warning: the first character of atom %u in frame 1 is not capital.\n", \
                i_atom + 1);
        }
        if (tok[1] != '\0' && (tok[1] > 'z' || tok[1] < 'a'))
        {
            printf("# Warning: the second character of atom %u in frame 1 is not lower-case.\n", \
                i_atom + 1);
        }
        strncpy(atom_name[i_atom], tok, max_atom_name_space);
    }
    fseek(traj_fp, 0, SEEK_SET);

    /* read each frame */
    for (;;)
    {
        /* i_frame count starts from 1. */
        if (! fgets(buf, BUFSIZ, traj_fp) || ! * buf)
        {
            break;
        }
        ++ i_frame;
        if (sscanf(buf, "%u", & v_uint) != 1)
        {
            fprintf(stderr, "Error! Cannot read the amount of atoms from frame %u.\n", i_frame);
            exit(EXIT_FAILURE);
        }
        if (v_uint != num_atoms)
        {
            fprintf(stderr, "Error! Amount of atoms mismatches in frame %u and frame 1.\n", i_frame);
            exit(EXIT_FAILURE);
        }
        if (! fgets(buf, BUFSIZ, traj_fp))
        {
            fprintf(stderr, "Error! Cannot read the title line in frame %u.\n", i_frame);
            exit(EXIT_FAILURE);
        }
        for (i_atom = 0; i_atom < num_atoms; ++ i_atom)
        {
            if (! fgets(buf, BUFSIZ, traj_fp))
            {
                fprintf(stderr, "Error! Cannot read atom %u in frame %u.\n", i_atom + 1, i_frame);
                exit(EXIT_FAILURE);
            }
            tok = strtok(buf, " ");
            if (! tok || sscanf(tok, "%lf", & v_double))
            {
                fprintf(stderr, "Error! Cannot read atom %u in frame %u.\n", i_atom + 1, i_frame);
                exit(EXIT_FAILURE);
            }
            if (strlen(tok) >= max_atom_name_space && i_frame > 1)
            {
                printf("# Warning: the name length of atom %u in frame %u is larger than %u.\n", \
                    i_atom + 1, i_frame, max_atom_name_space - 1);
                printf("# The first %u characters will be used.\n", max_atom_name_space - 1);
            }
            if (strncmp(tok, atom_name[i_atom], max_atom_name_space))
            {
                fprintf(stderr, "Error! The name of atom %u in frame %u mismatches that in frame 1.\n", \
                    i_atom + 1, i_frame);
                exit(EXIT_FAILURE);
            }
            for (i_coord = 0; i_coord < num_coords; ++ i_coord)
            {
                tok = strtok(NULL, " ");
                if (! tok || sscanf(tok, "%lf", atom_coords[i_atom] + i_coord) != 1)
                {
                    fprintf(stderr, "Error! Cannot read %c coordinate of atom %u in frame %u.\n", \
                        i_coord ? (i_coord > 1 ? 'z' : 'y') : 'x', i_atom + 1, i_frame);
                    exit(EXIT_FAILURE);
                }
            }
        }

        /* output geometry parameters of this frame. */
        Output_parameters(output_fp, num_atoms, (double const *const *)atom_coords, num_paras, \
            (unsigned int const *const *)para_index);
    }

    /* release memory */
    free(* para_index);
    memset(para_index, 0, num_paras * sizeof(unsigned int *));
    free(para_index);
    para_index = NULL;
    free(* atom_name);
    memset(atom_name, 0, num_atoms * sizeof(char *));
    free(atom_name);
    atom_name = NULL;
    free(* atom_coords);
    memset(atom_coords, 0, num_atoms * sizeof(double *));
    free(atom_coords);
    atom_coords = NULL;

    /* close trajectory file */
    fclose(traj_fp);
    traj_fp = NULL;

    /* close output file if needed */
    if (* output_name && strcmp(output_name, "-"))
    {
        fclose(output_fp);
        output_fp = NULL;
    }

    return 0;
}

/* note: below double * is a double[num_coords] */

inline double Radian_to_degree(double r)
{
    return r / M_PI * 180;
}

void Vector_minus(double const *a, double const *b, double *res)
{
    unsigned int i = 0u;

    for (i = 0u; i < num_coords; ++ i)
    {
        res[i] = a[i] - b[i];
    }

    return;
}

double Get_vector_dot_product(double const *a, double const *b)
{
    unsigned int i = 0u;
    double ret = 0.0;

    for (i = 0u; i < num_coords; ++ i)
    {
        ret += a[i] * b[i];
    }

    return ret;
}

double Get_vector_length(double const *a)
{
    return sqrt(Get_vector_dot_product(a, a));
}

void Normalize_vector(double *a)
{
    unsigned int i = 0u;
    double vec_len = Get_vector_length(a);

    for (i = 0u; i < num_coords; ++ i)
    {
        a[i] /= vec_len;
    }

    return;
}

double Get_distance(double const *a, double const *b)
{
    double v[num_coords] = {0.0};

    Vector_minus(b, a, v);
    return Get_vector_length(v);
}

double Get_vector_angle(double const *a, double const *b)
{
    return Radian_to_degree(acos(Get_vector_dot_product(a, b) / (Get_vector_length(a) * Get_vector_length(b))));
}

double Get_angle(double const *a, double const *b, double const *c)
{
    double v_1[num_coords] = {0.0};
    double v_2[num_coords] = {0.0};

    Vector_minus(b, a, v_1);
    Vector_minus(b, c, v_2);

    return Get_vector_angle(v_1, v_2);
}

void Vector_cross_product(double const *a, double const *b, double *res)
{
    unsigned int i = 0u;
    unsigned int j = 0u;
    unsigned int k = 0u;

    for (k = 0u; k < num_coords; ++ k)
    {
        i = (k + 1) % num_coords;
        j = (i + 1) % num_coords;
        res[k] = a[i] * b[j] - a[j] * b[i];
    }

    return;
}

double Get_dihedral(double const *a, double const *b, double const *c, double const *d)
{
    double v_1[num_coords] = {0.0};
    double v_2[num_coords] = {0.0};
    double v_3[num_coords] = {0.0};
    double n_1[num_coords] = {0.0};
    double n_2[num_coords] = {0.0};
    double m[num_coords] = {0.0};
    double x = 0.0, y = 0.0;

    Vector_minus(b, a, v_1);
    Vector_minus(c, b, v_2);
    Vector_minus(d, c, v_3);

    Normalize_vector(v_1);
    Normalize_vector(v_2);
    Normalize_vector(v_3);

    Vector_cross_product(v_1, v_2, n_1);
    Vector_cross_product(v_2, v_3, n_2);

    if (Get_vector_length(n_1) <= DBL_EPSILON * 1.E3 || \
        Get_vector_length(n_2) <= DBL_EPSILON * 1.E3)
    {
        /* a-b-c or b-c-d is collinear */
        return (double)NAN;
    }

    Vector_cross_product(n_1, n_2, m);
    y = Get_vector_dot_product(m, v_2);
    /* y /= Get_vector_length(v_2); */ /* if v_2 was not normalized then this is needed */
    x = Get_vector_dot_product(n_1, n_2);

    if (fabs(x) <= DBL_EPSILON * 1.E3 && fabs(y) <= DBL_EPSILON * 1.E3)
    {
        /* arctan2(0, 0) is not defined */
        return (double)NAN;
    }

    return Radian_to_degree(atan2(y, x));
}

void Output_parameters(FILE *out_file_ptr, unsigned int num_atoms, double const *const *atom_coords, \
                       unsigned int num_paras, unsigned int const *const *para_index)
{
    unsigned int i_para = 0u;

    fprintf(out_file_ptr, "%4s", "");
    fflush(out_file_ptr);
    for (i_para = 0u; i_para < num_paras; ++ i_para)
    {
        if (i_para)
        {
            fprintf(out_file_ptr, "%9s", "");
            fflush(out_file_ptr);
        }
        if (para_index[i_para][3])
        {
            /* dihedral */
            fprintf(out_file_ptr, "%10.5lf", Get_dihedral(atom_coords[para_index[i_para][0] - 1], \
                                                          atom_coords[para_index[i_para][1] - 1], \
                                                          atom_coords[para_index[i_para][2] - 1], \
                                                          atom_coords[para_index[i_para][3] - 1]));
            fflush(out_file_ptr);
        }
        else
        {
            if (para_index[i_para][2])
            {
                /* angle */
                fprintf(out_file_ptr, " %9.5lf", Get_angle(atom_coords[para_index[i_para][0] - 1], \
                                                           atom_coords[para_index[i_para][1] - 1], \
                                                           atom_coords[para_index[i_para][2] - 1]));
                fflush(out_file_ptr);
            }
            else
            {
                /* length */
                fprintf(out_file_ptr, "%10.6lf", Get_distance(atom_coords[para_index[i_para][0] - 1], \
                                                              atom_coords[para_index[i_para][1] - 1]));
                fflush(out_file_ptr);
            }
        }
    }
    fprintf(out_file_ptr, "\n");

    return;
}

