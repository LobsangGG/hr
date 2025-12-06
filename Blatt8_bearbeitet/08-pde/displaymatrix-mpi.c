#include <inttypes.h>
#include <malloc.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "mpi.h"
#include "partdiff.h"


void displayMatrixMPI(struct calculation_arguments *arguments,
                             struct calculation_results *results,
                             struct options *options, int rank, int size,
                             int from, int to) {
  int const elements = 8 * options->interlines + 9;

  int x, y;
  double **Matrix = arguments->Matrix[results->m];
  MPI_Status status;

  /* first line belongs to rank 0 */
  if (rank == 0)
    from--;
                              
  /* last line belongs to rank size - 1 */
  if (rank + 1 == size)
    to++;

  if (rank == 0)
    printf("Matrix:\n");

  for (y = 0; y < 9; y++) {
    int line = y * (options->interlines + 1);

    if (rank == 0) {
      /* check whether this line belongs to rank 0 */
      if (line < from || line > to) {
        /* use the tag to receive the lines in the correct order
         * the line is stored in Matrix[0], because we do not need it anymore */
        MPI_Recv(Matrix[0], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y,
                 arguments->comm, &status);      
      }
    } else {
      if (line >= from && line <= to) {
        /* if the line belongs to this process, send it to rank 0
         * (line - from + 1) is used to calculate the correct local address */
        MPI_Send(Matrix[line - from + 1], elements, MPI_DOUBLE, 0, 42 + y,
                 arguments->comm);
                      }
    }

    if (rank == 0) {
      for (x = 0; x < 9; x++) {
        int col = x * (options->interlines + 1);

        if (line >= from && line <= to) {
          /* this line belongs to rank 0 */
          printf("%11.8f", Matrix[line][col]);
        } else {
          /* this line belongs to another rank and was received above */
          printf("%11.8f", Matrix[0][col]);
        }
      }

      printf("\n");
    }
  }

  fflush(stdout);
}
