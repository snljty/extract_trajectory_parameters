/*************************************************************************************
 * This program extracts geometry parameters from a multi-frame xyz trajectory file. *
 * Note: the amount of atoms and their types are assumed to be unique in any frames. * 
 *************************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include <stdbool.h>

# define num_coords 3u /* x y z */

void Vector_minus(double const *a, double const *b, double *res);
double Get_vector_dot_product(double const *a, double const *b);
double Get_vector_length(double const *a);
double Get_distance(double const *a, double const *b);
double Get_vector_angle(double const *a, double const *b);
double Get_angle(double const *a, double const *b, double const *c);
void Vector_cross_product(double const *a, double const *b, double *res);
double Get_dihedral(double const *a, double const *b, double const *c, double const *d);
void Output_parameters(unsigned int num_atoms, double const *const *atom_coords, \
                       unsigned int num_paras, unsigned int const *const *para_index);

int main(int argc, char const *argv[])
{
    int iarg = 0;
    char traj_name[BUFSIZ + 1] = "";
    char index_name[BUFSIZ + 1] = "";
    char *traj_name_use = traj_name;
    char *index_name_use = index_name;
    FILE *traj_fp = NULL;
    FILE *index_fp = NULL;
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
    double tmp_double = 0.0;
    unsigned int tmp_uint = 0u;
    char tmp_str[24] = "";

    /* Get file names and open files */
    iarg = 1;
    while (iarg < argc)
    {
        if (! strcmp(argv[iarg], "-h") || ! strcmp(argv[iarg], "--help") || ! strcmp(argv[iarg], "/?"))
        {
            printf("Usage: %s [-h | --help]\n", argv[0]);
            printf("       %s [-t | --traj TRAJ] [-i | --index INDEX]\n", argv[0]);
            puts("");
            puts("TRAJ: name of a xyz format trajectory file with multi frames.");
            puts("INDEX: name of an index file.");
            puts("");
            puts("The index file should contains several lines, each line contains 2 to 4 integers.");
            puts("Those integers are the indices of the atoms, if 2 integers are provided, the ");
            puts("bond length of these two atoms will be calcualted, with the unit of input, ");
            puts("3 or 4 integers stands for angle (degree) and dihedral angle (degree).");
            puts("");
            puts("All file names that are not provided in the command argument, ");
            puts("will be asked interactively.");
            puts("If INDEX is \"-\", the program will read indices from STDIN.");
            puts("In this mode, you need to input the amount of parameters ro be measured first, ");
            puts("when the screen hints you to do so.");
            puts("");
            puts("Exiting normally.");
            exit(EXIT_SUCCESS);
        }
        else if (! strcmp(argv[iarg], "-t") || ! strcmp(argv[iarg], "--traj"))
        {
            ++ iarg;
            if (iarg >= argc)
            {
                fprintf(stderr, "Error! missing argument for \"%s\"\n.", argv[iarg - 1]);
                exit(EXIT_FAILURE);
            }
            strncpy(traj_name, argv[iarg], BUFSIZ + 1);
        }
        else if (! strcmp(argv[iarg], "-i") || ! strcmp(argv[iarg], "--index"))
        {
            ++ iarg;
            if (iarg >= argc)
            {
                fprintf(stderr, "Error! missing argument for \"%s\"\n.", argv[iarg - 1]);
                exit(EXIT_FAILURE);
            }
            strncpy(index_name, argv[iarg], BUFSIZ + 1);
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
        puts("# Input name of the trajectory file:");
        if (! fgets(traj_name, BUFSIZ, stdin))
            exit(EXIT_FAILURE);
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
        puts("# Input name of the index file:");
        puts("# If you press <Enter> directly, will read indices from STDIN.");
        if (! fgets(index_name, BUFSIZ, stdin))
            exit(EXIT_FAILURE);
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
        puts("# Reading indices from STDIN, please input ...");
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

    /* get amount of parameters */
    if (read_indices_from_stdin)
    {
        printf("# Input amount of geometry parameters to be measured: ");
        for (;;)
        {
            if (! fgets(buf, BUFSIZ, stdin))
                exit(EXIT_FAILURE);
            if (sscanf(buf, "%u", & num_paras) == 1)
                break;
        }
        if (num_paras == 1)
            puts("# Input 1 line of indices:");
        else
            printf("# Input %u lines of indices, one line for one geometry parameter:\n", num_paras);
    }
    else
    {
        /* get amount of parameters */
        while (fgets(buf, BUFSIZ, index_fp))
        {
            tok = strtok(buf, " \n,");
            if (tok)
                ++ num_paras;
        }
        fseek(index_fp, 0, SEEK_SET);
    }

    /* allocate memory for indices */
    para_index = (unsigned int **)malloc(num_paras * sizeof(unsigned int *));
    * para_index = (unsigned int *)malloc(num_paras * max_num_index * sizeof(unsigned int));
    memset(* para_index, 0, num_paras * max_num_index * sizeof(unsigned int));
    for (i_para = 0; i_para < num_paras; ++ i_para)
        para_index[i_para] = * para_index + i_para * max_num_index;

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
    puts("# bond lengthes in input unit, angles and dihedrals in unit degree.");
    printf("# ");
    for (i_para = 0; i_para < num_paras; ++ i_para)
    {
        if (i_para)
            printf("%4s", "");
        if (para_index[i_para][3])
        {
            /* dihedral */
            sprintf(tmp_str, "%u-%u-%u-%u", para_index[i_para][0], para_index[i_para][1], \
                para_index[i_para][2], para_index[i_para][3]);
            tmp_uint = strlen(tmp_str);
            printf("%*s%*s%*s", (15 - tmp_uint) / 2 >= 0 ? (15 - tmp_uint) / 2 : 0, "", \
                tmp_uint, tmp_str, \
                (16 - tmp_uint) / 2 >= 0 ? (16 - tmp_uint) / 2 : 0, "");
        }
        else
        {
            if (para_index[i_para][2])
            {
                /* angle */
                sprintf(tmp_str, "%u-%u-%u", para_index[i_para][0], \
                    para_index[i_para][1], para_index[i_para][2]);
                tmp_uint = strlen(tmp_str);
                printf("%*s%*s%*s", (15 - tmp_uint) / 2 >= 0 ? (15 - tmp_uint) / 2 : 0, "", \
                    tmp_uint, tmp_str, \
                    (16 - tmp_uint) / 2 >= 0 ? (16 - tmp_uint) / 2 : 0, "");
            }
            else
            {
                /* length */
                sprintf(tmp_str, "%u-%u", para_index[i_para][0], para_index[i_para][1]);
                tmp_uint = strlen(tmp_str);
                printf("%*s%*s%*s", (15 - tmp_uint) / 2 >= 0 ? (15 - tmp_uint) / 2 : 0, "", \
                        tmp_uint, tmp_str, \
                        (16 - tmp_uint) / 2 >= 0 ? (16 - tmp_uint) / 2 : 0, "");
            }
        }
    }
    puts("");

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
            printf("# Warning: the first character of atom %u in frame 1 is not capital.\n", \
                i_atom + 1);
        if (tok[1] != '\0' && (tok[1] > 'z' || tok[1] < 'a'))
            printf("# Warning: the second character of atom %u in frame 1 is not lower-case.\n", \
                i_atom + 1);
        strncpy(atom_name[i_atom], tok, max_atom_name_space);
    }
    fseek(traj_fp, 0, SEEK_SET);

    /* read each frame */
    for (;;)
    {
        /* i_frame count starts from 1. */
        if (! fgets(buf, BUFSIZ, traj_fp) || ! * buf)
            break;
        ++ i_frame;
        if (sscanf(buf, "%u", & tmp_uint) != 1)
        {
            fprintf(stderr, "Error! Cannot read the amount of atoms from frame %u\n.", i_frame);
            exit(EXIT_FAILURE);
        }
        if (tmp_uint != num_atoms)
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
            if (! tok || sscanf(tok, "%lf", & tmp_double))
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
        Output_parameters(num_atoms, (double const *const *)atom_coords, num_paras, \
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

    return 0;
}

/* note: below double * if a double[num_coords] */

void Vector_minus(double const *a, double const *b, double *res)
{
    unsigned int i = 0u;

    for (i = 0u; i < num_coords; ++ i)
        res[i] = a[i] - b[i];

    return;
}

double Get_vector_dot_product(double const *a, double const *b)
{
    unsigned int i = 0u;
    double ret = 0.0;

    ret = 0.0;
    for (i = 0u; i < num_coords; ++ i)
        ret += a[i] * b[i];

    return ret;
}

double Get_vector_length(double const *a)
{
    return sqrt(Get_vector_dot_product(a, a));
}

double Get_distance(double const *a, double const *b)
{
    double tmp[num_coords] = {0.0};

    Vector_minus(b, a, tmp);
    return Get_vector_length(tmp);
}

double Get_vector_angle(double const *a, double const *b)
{
    return acos(Get_vector_dot_product(a, b) / (Get_vector_length(a) * Get_vector_length(b))) / M_PI * 180;
}

double Get_angle(double const *a, double const *b, double const *c)
{
    double tmp_1[num_coords] = {0.0};
    double tmp_2[num_coords] = {0.0};

    Vector_minus(b, a, tmp_1);
    Vector_minus(b, c, tmp_2);

    return Get_vector_angle(tmp_1, tmp_2);
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
    double tmp_1[num_coords] = {0.0};
    double tmp_2[num_coords] = {0.0};
    double tmp_3[num_coords] = {0.0};
    double tmp_p[num_coords] = {0.0};
    double tmp_q[num_coords] = {0.0};
    double ret = 0.0;

    Vector_minus(b, a, tmp_1);
    Vector_minus(b, c, tmp_2);
    Vector_minus(c, d, tmp_3);
    Vector_cross_product(tmp_1, tmp_2, tmp_p);
    Vector_cross_product(tmp_3, tmp_2, tmp_q);
    ret = Get_vector_angle(tmp_p, tmp_q);
    if (Get_vector_dot_product(tmp_p, tmp_3) < 0.0)
        ret = - ret;

    return ret;
}

void Output_parameters(unsigned int num_atoms, double const *const *atom_coords, \
                       unsigned int num_paras, unsigned int const *const *para_index)
{
    unsigned int i_para = 0u;

    printf("%4s", "");
    for (i_para = 0u; i_para < num_paras; ++ i_para)
    {
        if (i_para)
            printf("%9s", "");
        if (para_index[i_para][3])
        {
            /* dihedral */
            printf("%10.5lf", Get_dihedral(atom_coords[para_index[i_para][0] - 1], \
                                           atom_coords[para_index[i_para][1] - 1], \
                                           atom_coords[para_index[i_para][2] - 1], \
                                           atom_coords[para_index[i_para][3] - 1]));
        }
        else
        {
            if (para_index[i_para][2])
            {
                /* angle */
                printf("%10.5lf", Get_angle(atom_coords[para_index[i_para][0] - 1], \
                                            atom_coords[para_index[i_para][1] - 1], \
                                            atom_coords[para_index[i_para][2] - 1]));
            }
            else
            {
                /* length */
                printf("%10.6lf", Get_distance(atom_coords[para_index[i_para][0] - 1], \
                                               atom_coords[para_index[i_para][1] - 1]));
            }
        }
    }
    puts("");

    return;
}

