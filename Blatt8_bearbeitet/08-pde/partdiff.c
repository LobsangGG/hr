/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                 TU München - Institut für Informatik                   **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff.c                                                  **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauß-Seidel and    **/
/**            Jacobi method.                                              **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#define _POSIX_C_SOURCE 200809L

#include <inttypes.h>
#include <malloc.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "mpi.h"
#include "partdiff.h"
#include <omp.h>

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time; /* time when program started                      */
struct timeval comp_time;  /* time when calculation completed                */

/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static void initVariables(struct calculation_arguments *arguments,
                          struct calculation_results *results,
                          struct options const *options) {
  
  arguments->N_columns = (options->interlines * 8) + 9 - 1;

  if (options->method == METH_JACOBI) {

    arguments->N_rows = ((arguments->N_columns - 1) / arguments->size) + 1;

    if ((uint64_t) arguments->rank < (arguments->N_columns - 1) % arguments->size) {
      arguments->N_rows += 1;
    }  
  }
  else {
    arguments->N_rows = arguments->N_columns;
  }
  arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
  arguments->h = 1.0 / arguments->N_columns;

  results->m = 0;
  results->stat_iteration = 0;
  results->stat_precision = 0;
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static void freeMatrices(struct calculation_arguments *arguments) {
  uint64_t i;

  for (i = 0; i < arguments->num_matrices; i++) {
    free(arguments->Matrix[i]);
  }

  free(arguments->Matrix);
  free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static void *allocateMemory(size_t size) {
  void *p;

  if ((p = malloc(size)) == NULL) {
    printf("Speicherprobleme! (%" PRIu64 " Bytes angefordert)\n", size);
    exit(1);
  }

  return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static void allocateMatrices(struct calculation_arguments *arguments) {
  uint64_t i, j;
  uint64_t const N_colums = arguments->N_columns;
  uint64_t const N_rows = arguments->N_rows;

  arguments->M = allocateMemory(arguments->num_matrices * (N_rows + 1) * (N_colums + 1) *
                              sizeof(double));
  arguments->Matrix =
      allocateMemory(arguments->num_matrices * sizeof(double **));

  for (i = 0; i < arguments->num_matrices; i++) {
    arguments->Matrix[i] = allocateMemory((N_rows + 1) * sizeof(double *));

    for (j = 0; j <= N_rows; j++) {
      arguments->Matrix[i][j] =
          arguments->M + (i * (N_rows + 1) * (N_colums + 1)) + (j * (N_colums + 1));
    }
  }
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static void initMatrices(struct calculation_arguments *arguments,
                         struct options const *options) {
  uint64_t g, i, j; /* local variables for loops */

  uint64_t const N_columns = arguments->N_columns;
  uint64_t const N_rows = arguments->N_rows;
  double const h = arguments->h;
  double ***Matrix = arguments->Matrix;

  /* initialize matrix/matrices with zeros */
  for (g = 0; g < arguments->num_matrices; g++) {
    for (i = 0; i <= N_rows; i++) {
      for (j = 0; j <= N_columns; j++) {
        Matrix[g][i][j] = 0.0;
      }
    }
  }

  /* initialize borders, depending on function (function 2: nothing to do) */
  if (options->inf_func == FUNC_F0) {
    for (g = 0; g < arguments->num_matrices; g++) {
      for (i = 0; i <= N_columns; i++) {
        Matrix[g][i][0] = 3 + (1 - (h * i));                  // Linke Kante
        Matrix[g][N_rows][i] = 3 - (h * i);                   // Untere Kante
        Matrix[g][N_columns - i][N_columns] = 2 + h * i;      // Rechte Kante
        Matrix[g][0][N_rows - i] = 3 + h * i;                 // Obere Kante
      }
    }
  }
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static void calculateGS(struct calculation_arguments const *arguments,
                      struct calculation_results *results,
                      struct options const *options) {
  int i, j;           /* local variables for loops */
  int m1, m2;         /* used as indices for old and new matrices */
  double star;        /* four times center value minus 4 neigh.b values */
  double residuum;    /* residuum of current iteration */
  double maxResiduum; /* maximum residuum value of a slave in iteration */

  int const N_columns = arguments->N_columns;
  int const N_rows = arguments->N_rows;
  double const h = arguments->h;

  double pih = 0.0;
  double fpisin = 0.0;

  int term_iteration = options->term_iteration;

  /* initialize m1 and m2 depending on algorithm */

  m1 = 0;
  m2 = 0;
  

  if (options->inf_func == FUNC_FPISIN) {
    pih = PI * h;
    fpisin = 0.25 * TWO_PI_SQUARE * h * h;
  }

  while (term_iteration > 0) {
    double **Matrix_Out = arguments->Matrix[m1];
    double **Matrix_In = arguments->Matrix[m2];

    maxResiduum = 0;

    /* over all rows */
    for (i = 1; i < N_rows; i++) {
      double fpisin_i = 0.0;

      if (options->inf_func == FUNC_FPISIN) {
        fpisin_i = fpisin * sin(pih * (double)i);
      }

      /* over all columns */
      for (j = 1; j < N_columns; j++) {
        star = 0.25 * (Matrix_In[i - 1][j] + Matrix_In[i][j - 1] +
                       Matrix_In[i][j + 1] + Matrix_In[i + 1][j]);

        if (options->inf_func == FUNC_FPISIN) {
          star += fpisin_i * sin(pih * (double)j);
        }

        if (options->termination == TERM_PREC || term_iteration == 1) {
          residuum = Matrix_In[i][j] - star;
          residuum = (residuum < 0) ? -residuum : residuum;
          maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
        }

        Matrix_Out[i][j] = star;
      }
    }

    results->stat_iteration++;
    results->stat_precision = maxResiduum;

    /* exchange m1 and m2 */
    i = m1;
    m1 = m2;
    m2 = i;

    /* check for stopping calculation depending on termination method */
    if (options->termination == TERM_PREC) {
      if (maxResiduum < options->term_precision) {
        term_iteration = 0;
      }
    } else if (options->termination == TERM_ITER) {
      term_iteration--;
    }
  }

  results->m = m2;
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static void calculateJacobi(struct calculation_arguments const *arguments,
                      struct calculation_results *results,
                      struct options const *options) {
  int i, j;           /* local variables for loops */
  int m1, m2;         /* used as indices for old and new matrices */
  double star;        /* four times center value minus 4 neigh.b values */
  double residuum;    /* residuum of current iteration */
  double maxResiduum; /* maximum residuum value of a slave in iteration */

  int const N_columns = arguments->N_columns;
  int const N_rows = arguments->N_rows;
  double const h = arguments->h;

  double pih = 0.0;
  double fpisin = 0.0;

  int term_iteration = options->term_iteration;

  /* initialize m1 and m2 depending on algorithm */
  m1 = 0;
  m2 = 1;

  if (options->inf_func == FUNC_FPISIN) {
    pih = PI * h;
    fpisin = 0.25 * TWO_PI_SQUARE * h * h;
  }

  omp_set_num_threads(options->number);
  while (term_iteration > 0) {
    double **Matrix_Out = arguments->Matrix[m1];
    double **Matrix_In = arguments->Matrix[m2];

    maxResiduum = 0;


    int from = arguments->rank * (arguments->N_rows - 1) + 1;
    if ((uint64_t)arguments->rank >= (arguments->N_columns - 1) % arguments->size) {
      from += (arguments->N_columns - 1) % arguments->size;
    }

    #pragma omp parallel for if (options->method == METH_JACOBI) default(none) private(i, j, residuum, star) shared(from, Matrix_In, Matrix_Out, N_rows, N_columns, fpisin, pih, options, term_iteration) reduction(max:maxResiduum)

    /* over all rows */
    for (i = 1; i < N_rows; i++) {
      double fpisin_i = 0.0;

      if (options->inf_func == FUNC_FPISIN) {

        fpisin_i = fpisin * sin(pih * (double)(from - 1 + i));
      }

      /* over all columns */
      for (j = 1; j < N_columns; j++) {
        star = 0.25 * (Matrix_In[i - 1][j] + Matrix_In[i][j - 1] +
                       Matrix_In[i][j + 1] + Matrix_In[i + 1][j]);

        if (options->inf_func == FUNC_FPISIN) {
          star += fpisin_i * sin(pih * (double)j);
        }

        if (options->termination == TERM_PREC || term_iteration == 1) {
          residuum = Matrix_In[i][j] - star;
          residuum = (residuum < 0) ? -residuum : residuum;
          maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
        }

        Matrix_Out[i][j] = star;
      }
    }

    if (arguments->size > 1) {
      /* communicate border rows with neighboring processes */
      MPI_Request requests[4];
      int nreqs = 0;
      MPI_Status statuses[4];

      int tag_TopRow = 1;
      int tag_BottomRow = 0;

      /* receive from previous and next */
      if (arguments->rank > 0) {
        MPI_Irecv(Matrix_Out[0], N_columns + 1, MPI_DOUBLE, arguments->rank - 1, tag_BottomRow + term_iteration,
                  arguments->comm, &requests[nreqs++]);
      }
      if (arguments->rank < arguments->size - 1) {
        MPI_Irecv(Matrix_Out[N_rows], N_columns + 1, MPI_DOUBLE, arguments->rank + 1, tag_TopRow + term_iteration,
                  arguments->comm, &requests[nreqs++]);
      }


      /* send to previous and next */
      if (arguments->rank > 0) {
        MPI_Issend(Matrix_Out[1], N_columns + 1, MPI_DOUBLE, arguments->rank - 1, tag_TopRow + term_iteration,
                  arguments->comm, &requests[nreqs++]);
      }
      if (arguments->rank < arguments->size - 1) {
        MPI_Issend(Matrix_Out[N_rows - 1], N_columns + 1, MPI_DOUBLE, arguments->rank + 1, tag_BottomRow + term_iteration,
                  arguments->comm, &requests[nreqs++]);
      }

      /* wait for communication to finish */
      MPI_Waitall(nreqs, requests, statuses);
    }

    results->stat_iteration++;
    results->stat_precision = maxResiduum;

    MPI_Barrier(arguments->comm);

    /* exchange m1 and m2 */
    i = m1;
    m1 = m2;
    m2 = i;

    /* check for stopping calculation depending on termination method */
    if (options->termination == TERM_PREC) {
      if (options->method == METH_JACOBI && arguments->size > 1) {
        double global_maxResiduum;
        MPI_Allreduce(&maxResiduum, &global_maxResiduum, 1, MPI_DOUBLE, MPI_MAX, arguments->comm);
        maxResiduum = global_maxResiduum;
      }
      if (maxResiduum < options->term_precision) {
        term_iteration = 0;
      }
    } else if (options->termination == TERM_ITER) {
      term_iteration--;
    } 
    MPI_Barrier(arguments->comm);
  }
  if (options->termination == TERM_ITER) {
    double global_maxResiduum;
    MPI_Allreduce(&maxResiduum, &global_maxResiduum, 1, MPI_DOUBLE, MPI_MAX, arguments->comm);
    results->stat_precision = global_maxResiduum;
  }
  results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static void displayStatistics(struct calculation_arguments const *arguments,
                              struct calculation_results const *results,
                              struct options const *options) {
  int N_columns = arguments->N_columns;
  int N_rows = arguments->N_rows;
  double time = (comp_time.tv_sec - start_time.tv_sec) +
                (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

  printf("Berechnungszeit:    %f s \n", time);
  printf("Speicherbedarf:     %f MiB\n", arguments->size * (N_rows + 1) * (N_columns + 1) * sizeof(double) *
                                             arguments->num_matrices / 1024.0 /
                                             1024.0);
  printf("Berechnungsmethode: ");

  if (options->method == METH_GAUSS_SEIDEL) {
    printf("Gauß-Seidel");
  } else if (options->method == METH_JACOBI) {
    printf("Jacobi");
  }

  printf("\n");
  printf("Interlines:         %" PRIu64 "\n", options->interlines);
  printf("Stoerfunktion:      ");

  if (options->inf_func == FUNC_F0) {
    printf("f(x,y) = 0");
  } else if (options->inf_func == FUNC_FPISIN) {
    printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
  }

  printf("\n");
  printf("Terminierung:       ");

  if (options->termination == TERM_PREC) {
    printf("Hinreichende Genaugkeit");
  } else if (options->termination == TERM_ITER) {
    printf("Anzahl der Iterationen");
  }

  printf("\n");
  printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
  printf("Norm des Fehlers:   %.11e\n", results->stat_precision);
  printf("\n");
}

/****************************************************************************/
/** Beschreibung der Funktion displayMatrix:                               **/
/**                                                                        **/
/** Die Funktion displayMatrix gibt eine Matrix                            **/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**                                                                        **/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.                                      **/
/****************************************************************************/
static void displayMatrix(struct calculation_arguments *arguments,
                          struct calculation_results *results,
                          struct options *options) {
  int x, y;

  double **Matrix = arguments->Matrix[results->m];

  int const interlines = options->interlines;

  printf("Matrix:\n");

  for (y = 0; y < 9; y++) {
    for (x = 0; x < 9; x++) {
      printf("%11.8f", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
    }

    printf("\n");
  }

  fflush(stdout);
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int main(int argc, char **argv) {

  struct options options;
  struct calculation_arguments arguments;
  struct calculation_results results;
  int is_master = 1;
  MPI_Comm newComm = MPI_COMM_NULL;
  int initial_rank;
  int initial_size;

  //Intialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &initial_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &initial_size);

  if (initial_rank == 0) {
    askParams(&options, argc, argv);
  }

  MPI_Bcast(&options, sizeof(struct options), MPI_BYTE, 0, MPI_COMM_WORLD);

  if (options.method == METH_JACOBI) {

    //Finalize processes that are not needed
    int total_rows = (options.interlines * 8) + 9;
    int color = initial_rank < (total_rows - 2) ? 0 : MPI_UNDEFINED ;

    if (initial_size > (total_rows - 2)) {
      MPI_Comm_split(MPI_COMM_WORLD, color, initial_rank, &newComm);
      if (color == MPI_UNDEFINED) {
        MPI_Finalize();
        return 0;
      }
    } else {
      MPI_Comm_dup(MPI_COMM_WORLD, &newComm);
    }

    MPI_Comm_rank(newComm, &arguments.rank);
    MPI_Comm_size(newComm, &arguments.size); 
    arguments.comm = newComm;
    is_master = arguments.rank == 0;
  } else {
    arguments.rank = 0;
    arguments.size = 1;
    arguments.comm = MPI_COMM_NULL;
  }

  initVariables(&arguments, &results, &options);

  allocateMatrices(&arguments);
  initMatrices(&arguments, &options);

  if (is_master) {
    gettimeofday(&start_time, NULL);
  }
  if (options.method == METH_GAUSS_SEIDEL) {
    calculateGS(&arguments, &results, &options);
  } else {
    calculateJacobi(&arguments, &results, &options);
  }
  MPI_Barrier(arguments.comm);
  if (is_master) {
    gettimeofday(&comp_time, NULL);
    displayStatistics(&arguments, &results, &options);
  }

  if  (options.method == METH_GAUSS_SEIDEL || arguments.size == 1) {
    displayMatrix(&arguments, &results, &options);
  } else {
    int from = arguments.rank * (arguments.N_rows - 1) + 1;
    if ((uint64_t)arguments.rank >= (arguments.N_columns - 1) % arguments.size) {
      from += (arguments.N_columns - 1) % arguments.size;
    }
    displayMatrixMPI(&arguments, &results, &options, arguments.rank, arguments.size,
                    from, from + arguments.N_rows - 2);
  }

  freeMatrices(&arguments);
  if (options.method == METH_JACOBI) {
    MPI_Comm_free(&(arguments.comm));
  }
  MPI_Finalize();

  return 0;
}
