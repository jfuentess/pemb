/******************************************************************************
 * main.c
 *
 * Parallel construction of succinct plane graphs
 *
 ******************************************************************************
 * Copyright (C) 2016 Leo Ferres, José Fuentes <jfuentess@dcc.uchile.cl>, Travis
 * Gagie, Meng He and Gonzalo Navarro
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 *****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

#ifdef MALLOC_COUNT
#include "malloc_count.h"
#endif

#include "planar_graph.h"

#ifdef POINTER_BASED
#include "pointer_based.h"
#endif

void print_partial_bit_array2(bitRankW32Int *, int32_t);

int main(int argc, char** argv) {

  if(argc < 2) {
    fprintf(stderr, "Usage: %s <input graph>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  struct timespec stime, etime;
  double time;

  Graph *g = read_graph_from_file(argv[1]);

  int *parent = malloc(g->n*sizeof(int));
  uint *count_edges = calloc(2*(g->n-1), sizeof(uint));
  uint *references = calloc(2*(g->n-1), sizeof(uint));
  uint init = 0;

#ifdef MALLOC_COUNT
  malloc_reset_peak();
#else
  if (clock_gettime(CLOCK_THREAD_CPUTIME_ID , &stime)) {
    fprintf(stderr, "clock_gettime failed");
    exit(-1);
  }
#endif

  /*
    t: spanning tree of the planar graph g
    parent: array of parents. parent[a] is the index of the edge (in the
    adjacency list) of node a to a's parent
    count_edges: Number of edges of g\t between two consecutive edges of
    t. count_edges[i] is the number of edges of g\t between i-th edge of t and
    its neighbor in clockwise order
  */

#ifdef NOPARALLEL
  Tree *t = seq_dfs_spanning_tree(g, init, parent, count_edges, references);
  SuccGraph *sg = seq_succ_planar_graph(g, t, init, parent, count_edges, references);
#else
  Tree *t = par_spanning_tree2(g, init, parent, count_edges, references);
  SuccGraph *sg = par_succ_planar_graph(g, t, init, parent, count_edges, references);
#endif

#ifdef MALLOC_COUNT
   printf("%s,%d,%d,%lu,%lu,%lu,%zu\n", argv[1], g->n, g->m, size_graph(g),
	  size_tree(t),  size_planar_graph(sg), malloc_count_peak());
#else
  if (clock_gettime(CLOCK_THREAD_CPUTIME_ID , &etime)) {
    fprintf(stderr, "clock_gettime failed");
    exit(-1);
  }
  
  time = (etime.tv_sec - stime.tv_sec) + (etime.tv_nsec - stime.tv_nsec) /
  1000000000.0;
  printf("%d,%s,%u,%u,%lf", threads, argv[1], g->n, g->m, time);
#endif

  return EXIT_SUCCESS;
}
