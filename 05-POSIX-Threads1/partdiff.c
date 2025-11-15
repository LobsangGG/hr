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
/**            Jacobi method (parallelized with POSIX threads).            **/
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
#include <pthread.h>

#include "partdiff.h"

struct calculation_arguments {
  uint64_t N;            /* number of spaces between lines (lines=N+1)     */
  uint64_t num_matrices; /* number of matrices                             */
  double h;              /* length of a space between two lines            */
  double ***Matrix;      /* index matrix used for addressing M             */
  double *M;             /* two matrices with real values                  */
};

struct calculation_results {
  uint64_t m;
  uint64_t stat_iteration; /* number of current iteration                    */
  double stat_precision;   /* actual precision of all slaves in iteration    */
};

/* ************************************************************************ */
/* Thread-specific structures                                               */
/* ************************************************************************ */

struct thread_data {
  int thread_id;
  int num_threads;
  const struct calculation_arguments *arguments;
  struct calculation_results *results;
  const struct options *options;
  int *m1;
  int *m2;
  int *term_iteration;
  double *maxResiduum;
  int *continue_iteration;
};

/* Global mutex and barrier for thread synchronization */
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_barrier_t barrier;

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
  arguments->N = (options->interlines * 8) + 9 - 1;
  arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
  arguments->h = 1.0 / arguments->N;

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
    printf("Speicherprobleme! (%" PRIu64 " Bytes angefordert)\n", (uint64_t)size);
    exit(1);
  }

  return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static void allocateMatrices(struct calculation_arguments *arguments) {
  uint64_t i, j;

  uint64_t const N = arguments->N;

  arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) *
                                sizeof(double));
  arguments->Matrix =
      allocateMemory(arguments->num_matrices * sizeof(double **));

  for (i = 0; i < arguments->num_matrices; i++) {
    arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double *));

    for (j = 0; j <= N; j++) {
      arguments->Matrix[i][j] =
          arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1));
    }
  }
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static void initMatrices(struct calculation_arguments *arguments,
                         struct options const *options) {
  uint64_t g, i, j; /* local variables for loops */

  uint64_t const N = arguments->N;
  double const h = arguments->h;
  double ***Matrix = arguments->Matrix;

  /* initialize matrix/matrices with zeros */
  for (g = 0; g < arguments->num_matrices; g++) {
    for (i = 0; i <= N; i++) {
      for (j = 0; j <= N; j++) {
        Matrix[g][i][j] = 0.0;
      }
    }
  }

  /* initialize borders, depending on function (function 2: nothing to do) */
  if (options->inf_func == FUNC_F0) {
    for (g = 0; g < arguments->num_matrices; g++) {
      for (i = 0; i <= N; i++) {
        Matrix[g][i][0] = 3 + (1 - (h * i)); // Linke Kante
        Matrix[g][N][i] = 3 - (h * i);       // Untere Kante
        Matrix[g][N - i][N] = 2 + h * i;     // Rechte Kante
        Matrix[g][0][N - i] = 3 + h * i;     // Obere Kante
      }
    }
  }
}

/* ************************************************************************ */
/* thread_calculate: Thread function for parallel Jacobi calculation        */
/* ************************************************************************ */
static void *thread_calculate(void *arg) {
  struct thread_data *data = (struct thread_data *)arg;
  
  int const N = data->arguments->N;
  double const h = data->arguments->h;
  
  double pih = 0.0;
  double fpisin = 0.0;
  
  if (data->options->inf_func == FUNC_FPISIN) {
    pih = PI * h;
    fpisin = 0.25 * TWO_PI_SQUARE * h * h;
  }
  
  /* berechnet die Reihen pro Thread */
  int total_rows = N - 1;  /* äußeren Zeilen nicht zum rechnen */
  int rows_per_thread = total_rows / data->num_threads;
  int rest = total_rows % data->num_threads;
  
  /* die ersten Threads gekommen den Rest */
  int start_row;
  int end_row;
  if (data->thread_id < rest) {
    start_row = 1 + data->thread_id * (rows_per_thread + 1);
    end_row = start_row + (rows_per_thread + 1);
  } else {
    start_row = 1 + rest * (rows_per_thread + 1) + 
                (data->thread_id - rest) * rows_per_thread;
    end_row = start_row + rows_per_thread;
  }
  
  while (*data->continue_iteration) {
    double **Matrix_Out = data->arguments->Matrix[*data->m1];
    double **Matrix_In = data->arguments->Matrix[*data->m2];
    
    double local_maxResiduum = 0.0;
    
    /* nur über Threads spezifische Reihen */
    for (int i = start_row; i < end_row; i++) {
      double fpisin_i = 0.0;
      
      if (data->options->inf_func == FUNC_FPISIN) {
        fpisin_i = fpisin * sin(pih * (double)i);
      }
      
      /* over all columns */
      for (int j = 1; j < N; j++) {
        double star = 0.25 * (Matrix_In[i - 1][j] + Matrix_In[i][j - 1] +
                              Matrix_In[i][j + 1] + Matrix_In[i + 1][j]);
        
        if (data->options->inf_func == FUNC_FPISIN) {
          star += fpisin_i * sin(pih * (double)j);
        }
        
        if (data->options->termination == TERM_PREC || *data->term_iteration == 1) {
          double residuum = Matrix_In[i][j] - star;
          residuum = (residuum < 0) ? -residuum : residuum;
          local_maxResiduum = (residuum < local_maxResiduum) ? local_maxResiduum : residuum;
        }
        
        Matrix_Out[i][j] = star;
      }
    }
    
    /* Update global maxResiduum */
    pthread_mutex_lock(&mutex);
    if (local_maxResiduum > *data->maxResiduum) {
      *data->maxResiduum = local_maxResiduum;
    }
    pthread_mutex_unlock(&mutex);
    
    /* Auf Threads warten */
    pthread_barrier_wait(&barrier);
    
    /* Nur ein Thread!  */
    if (data->thread_id == 0) {
      data->results->stat_iteration++;
      data->results->stat_precision = *data->maxResiduum;
      
      /* tausche m1 und m2 */
      int temp = *data->m1;
      *data->m1 = *data->m2;
      *data->m2 = temp;
      
      /* checked ob nach interation oder Genauigkeit abgebrochen wird */
      if (data->options->termination == TERM_PREC) {
        if (*data->maxResiduum < data->options->term_precision) {
          *data->continue_iteration = 0;
        }
      } else if (data->options->termination == TERM_ITER) {
        (*data->term_iteration)--;
        if (*data->term_iteration <= 0) {
          *data->continue_iteration = 0;
        }
      }
      
      /* Reset maxResiduum für nächste Iteration */
      *data->maxResiduum = 0.0;
    }
    
    /* alle Thread warten auf Thread 0 */
    pthread_barrier_wait(&barrier);
  }
  
  return NULL;
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static void calculate(struct calculation_arguments const *arguments,
                      struct calculation_results *results,
                      struct options const *options) {
  int m1, m2;         /* used as indices for old and new matrices */
  double maxResiduum; /* maximum residuum value of a slave in iteration */
  int term_iteration = options->term_iteration;
  int continue_iteration = 1;
  
  /* initialize m1 and m2 depending on algorithm */
  if (options->method == METH_JACOBI) {
    m1 = 0;
    m2 = 1;
  } else {
    m1 = 0;
    m2 = 0;
  }
  
  /* für Gauß-Seidel */
  if (options->method == METH_GAUSS_SEIDEL) {
    int i, j;
    int const N = arguments->N;
    double const h = arguments->h;
    double pih = 0.0;
    double fpisin = 0.0;
    
    if (options->inf_func == FUNC_FPISIN) {
      pih = PI * h;
      fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }
    
    while (term_iteration > 0) {
      double **Matrix_Out = arguments->Matrix[m1];
      double **Matrix_In = arguments->Matrix[m2];
      
      maxResiduum = 0;
      
      for (i = 1; i < N; i++) {
        double fpisin_i = 0.0;
        
        if (options->inf_func == FUNC_FPISIN) {
          fpisin_i = fpisin * sin(pih * (double)i);
        }
        
        for (j = 1; j < N; j++) {
          double star = 0.25 * (Matrix_In[i - 1][j] + Matrix_In[i][j - 1] +
                                Matrix_In[i][j + 1] + Matrix_In[i + 1][j]);
          
          if (options->inf_func == FUNC_FPISIN) {
            star += fpisin_i * sin(pih * (double)j);
          }
          
          if (options->termination == TERM_PREC || term_iteration == 1) {
            double residuum = Matrix_In[i][j] - star;
            residuum = (residuum < 0) ? -residuum : residuum;
            maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
          }
          
          Matrix_Out[i][j] = star;
        }
      }
      
      results->stat_iteration++;
      results->stat_precision = maxResiduum;
      
      i = m1;
      m1 = m2;
      m2 = i;
      
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
  /* für Jacobi */
  else {
    int num_threads = options->number;
    pthread_t *threads = allocateMemory(num_threads * sizeof(pthread_t));
    struct thread_data *thread_args = allocateMemory(num_threads * sizeof(struct thread_data));
    
    maxResiduum = 0.0;
    
    /* Barriere Initianisierung */
    pthread_barrier_init(&barrier, NULL, num_threads);
    
    /* Thread-Erstellung */
    for (int i = 0; i < num_threads; i++) {
      thread_args[i].thread_id = i;
      thread_args[i].num_threads = num_threads;
      thread_args[i].arguments = arguments;
      thread_args[i].results = results;
      thread_args[i].options = options;
      thread_args[i].m1 = &m1;
      thread_args[i].m2 = &m2;
      thread_args[i].term_iteration = &term_iteration;
      thread_args[i].maxResiduum = &maxResiduum;
      thread_args[i].continue_iteration = &continue_iteration;
      
      /* alle Threads führen thread_calculate mit &thread_args durch(Argumente für die Funktion) */
      pthread_create(&threads[i], NULL, thread_calculate, &thread_args[i]);
    }
    
    /* alle Threads joinen */
    for (int i = 0; i < num_threads; i++) {
      pthread_join(threads[i], NULL);
    }
    
    /* Speicher und Barrier freigeben */
    pthread_barrier_destroy(&barrier);
    free(threads);
    free(thread_args);
    
    results->m = m2;
  }
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static void displayStatistics(struct calculation_arguments const *arguments,
                              struct calculation_results const *results,
                              struct options const *options) {
  int N = arguments->N;
  double time = (comp_time.tv_sec - start_time.tv_sec) +
                (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

  printf("Berechnungszeit:    %f s \n", time);
  printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) *
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

  askParams(&options, argc, argv);

  initVariables(&arguments, &results, &options);

  allocateMatrices(&arguments);
  initMatrices(&arguments, &options);

  gettimeofday(&start_time, NULL);
  calculate(&arguments, &results, &options);
  gettimeofday(&comp_time, NULL);

  displayStatistics(&arguments, &results, &options);
  displayMatrix(&arguments, &results, &options);

  freeMatrices(&arguments);

  return 0;
}
