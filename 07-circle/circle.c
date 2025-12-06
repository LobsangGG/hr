#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"

int *init(int N, int is_full_arr, int rank, int size) {
  // TODO
  int *buf = (int *)malloc(sizeof(int) * N);
  if (rank >= size) {
    for (int i = 0; i < N; i++) {
      buf[i] = -1;
    }
  } else if (is_full_arr == 0) {
    N -= 1;
    buf[N] = -1;
  }

  srand(time(NULL) + rank);

  for (int i = 0; i < N; i++) {
    // Do not modify "% 10"
    buf[i] = rand() % 10;
  }

  return buf;
}

int circle(int *buf,int size, int rank, int N) {
  // TODO
  int is_last = rank == size - 1;
  int prev = (rank - 1 + size) % size;
  int next = (rank + 1) % size;
  MPI_Request requests[2*N];
  int first;
  int circle_not_over = 1;
  int iterations = 0;
  if (is_last) {
    MPI_Recv(&first, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else if (rank == 0) {
    first = buf[0];
    MPI_Ssend(&first, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD);
  }
  while (circle_not_over) {
    for (int i = 0; i < N; i++) {
      MPI_Issend(&buf[i], 1, MPI_INT, next, 0, MPI_COMM_WORLD, &requests[i]);
    }
    for (int i = 0; i < N; i++) {
      MPI_Irecv(&buf[i], 1, MPI_INT, prev, 0, MPI_COMM_WORLD, &requests[N + i]);
    }
    MPI_Waitall(2 * N, requests, MPI_STATUSES_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);
    if (is_last) {
      if (first == buf[0]) {
        circle_not_over = 0;
      }
    }
    MPI_Bcast(&circle_not_over, 1, MPI_INT, size - 1, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    iterations++;
  }
  return iterations;
}

void printArray(int *buf, int N, int size) {

        for (int i = 0; i < N; i++) {
            printf("rank %d: %d\n", 0, buf[i]);
        }
        if (size == 1) {
            return;
        }
        int *recvbuf = malloc(sizeof(int) * N);
        for (int src = 1; src < size; src++) {
            MPI_Recv(recvbuf, N, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < N; i++) {
                printf("rank %d: %d\n", src, recvbuf[i]);
            }
        }
        free(recvbuf);
} 


int main(int argc, char **argv) {
  int N;
  int rank;
  int size;
  int *buf;
  int ret;

  if (argc < 2) {
    printf("Arguments error!\nPlease specify a buffer size.\n");
    return EXIT_FAILURE;
  }

  MPI_Init( &argc, &argv ); 
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  int is_master = rank == 0;
  int is_full_arr = 1;

  // Array length
  N = atoi(argv[1]);
  if (N % size > 0) {
    is_full_arr = rank < N % size;
    N = N / size + 1;
  } else {
    N = N / size;
  }
  buf = init(N, is_full_arr, rank, size);

  MPI_Barrier(MPI_COMM_WORLD);

  // TODO
  if (is_master) {
    printf("\nBEFORE\n");

    printArray(buf, N, size);
  } else {
      MPI_Ssend(buf, N, MPI_INT, 0, 0, MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (size > 1) {
    ret = circle(buf, size, rank, N);
  } else {
    ret = 0;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (is_master) {
    printf("\nAFTER\n");

    printArray(buf, N, size);
    printf("\nNumber of iterations: %d\n", ret);
    int termination_value;
    MPI_Recv(&termination_value, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("\nTermination Value: %d\n", termination_value);

  } else {
      MPI_Ssend(buf, N, MPI_INT, 0, 0, MPI_COMM_WORLD);
      if (rank == size - 1) {
        MPI_Ssend(&buf[0], 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  free(buf);
  MPI_Finalize(); 
  return EXIT_SUCCESS;
}
