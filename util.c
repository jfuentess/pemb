/******************************************************************************
 * util.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util.h"

/*
  Format of the expected input/output file:
  
  <number of nodes>
  <number of edges>
  <source vertex> <target vertex>
  ....

  The list of edges (<source vertex> <target vertex>) must be ordered:
  - First, all edges where vertex 0 is the source vertex
  - Second, all edges where vertex 1 is the source vertex
  - and so on

  For a source vertex v, the corresponding target vertices must be in
  counterclockwise order

  Assumption: To support multiple edges, the indices of the multiple edges must
  be always increasing (in other words, the adjacency list of a node cannot
  start in the middle of a list of multiple edges
*/
Graph* read_graph_from_file(const char* fn) {
  Graph *g = malloc(sizeof(Graph));

  FILE *fp = fopen(fn, "r");
  char line_buffer[BUFSIZ]; /* BUFSIZ is defined if you include stdio.h */

  if (!fp) {
    fprintf(stderr, "Error opening file \"%s\".\n", fn);
    exit(EXIT_FAILURE);
  }

  g->n = atoi(fgets(line_buffer, sizeof(line_buffer), fp));
  g->m = atoi(fgets(line_buffer, sizeof(line_buffer), fp));

  g->V = calloc(g->n,sizeof(Vertex));
  g->E = calloc(2*(g->m),sizeof(Edge));

  for(uint i = 0; i < 2*g->m; i++)
    g->E[i].cmp = -1;
  
  uint source = 0, target = 0, m = 0;

  while (fgets(line_buffer, sizeof(line_buffer), fp)) {
    source = atoi(strtok(line_buffer, " "));
    target = atoi(strtok(NULL, " "));
    g->E[m].src = source;
    g->E[m].tgt = target;

    if(m==0)
      g->V[source].first = m;
    else if(source != g->E[m-1].src) {
      g->V[g->E[m-1].src].last = m-1;
      g->V[source].first = m;
    }
    m++;
  }
  g->V[g->E[m-1].src].last = m-1;

  fclose(fp);

  for(uint i = 0; i < 2*g->m; i++) {
    Vertex target = g->V[g->E[i].tgt];
    int cmp = -1;

    if(g->E[i].cmp != -1)
      continue;
    
    for(uint j = target.first; j <= target.last; j++) {
      /* Condition i != j allows to support self-loops */
      if((g->E[j].cmp == -1) && (g->E[j].tgt == g->E[i].src) && (i != j))
	cmp = j; // Choose the last unvisited (e.cmp=-1) edge, not the first one
    }

    if(cmp != -1) {
      g->E[i].cmp = cmp;
      g->E[cmp].cmp = i;
    }
  }

  return g;
}


void write_graph_to_file(const char* fn, Graph* g) {

  FILE* fp = fopen(fn, "w");

  if (!fp) {
    fprintf(stderr, "Error opening file \"%s\".\n", fn);
    exit(EXIT_FAILURE);
  }

  uint i = 0, j = 0;

  fprintf(fp, "%u\n", g->n);
  fprintf(fp, "%u\n", g->m);
  
  for(i = 0; i < g->n; i++)
    for(j = g->V[i].first; j <= g->V[i].last; j++)
      fprintf(fp, "%u %u\n", g->E[j].src, g->E[j].tgt);

  fclose(fp);
}


Tree* read_tree_from_file(const char* fn) {
  Tree *t = malloc(sizeof(Tree));

  FILE *fp = fopen(fn, "r");
  char line_buffer[BUFSIZ]; /* BUFSIZ is defined if you include stdio.h */

  if (!fp) {
    fprintf(stderr, "Error opening file \"%s\".\n", fn);
    exit(EXIT_FAILURE);
  }

  t->n = atoi(fgets(line_buffer, sizeof(line_buffer), fp));

  t->N = calloc(t->n,sizeof(Node));
  t->E = calloc(2*(t->n-1),sizeof(Edge));

  uint source = 0, target = 0, m = 0;

  while (fgets(line_buffer, sizeof(line_buffer), fp)) {
    source = atoi(strtok(line_buffer, " "));
    target = atoi(strtok(NULL, " "));
    t->E[m].src = source;
    t->E[m].tgt = target;

    if(m==0)
      t->N[source].first = m;
    else if(source != t->E[m-1].src) {
      t->N[t->E[m-1].src].last = m-1;
      t->N[source].first = m;
    }
    m++;
  }
  t->N[t->E[m-1].src].last = m-1;

  fclose(fp);

  uint i = 0, j = 0;

  for(i = 0; i < 2*(t->n-1); i++) {
    Node target = t->N[t->E[i].tgt];
    
    for(j = target.first; j <= target.last; j++)
      if(t->E[j].tgt == t->E[i].src) {
  	t->E[i].cmp = j;
  	break;
      }
  }

  return t;
}


void write_tree_to_file(const char* fn, Tree* t) {

  FILE* fp = fopen(fn, "w");

  if (!fp) {
    fprintf(stderr, "Error opening file \"%s\".\n", fn);
    exit(EXIT_FAILURE);
  }

  uint i = 0, j = 0;

  fprintf(fp, "%u\n", t->n);
  
  for(i = 0; i < t->n; i++)
    for(j = t->N[i].first; j <= t->N[i].last; j++)
      fprintf(fp, "%u %u\n", t->E[j].src, t->E[j].tgt);

  fclose(fp);
}


void print_graph(Graph* g) {

  if(g->m > 50) {
    fprintf(stderr, "Warning: The number of edges is greater than 50 edges. We\
  do not recommend to use this function for graphs with more than 50 edges.\n");
    return;
  }

  uint i, m = 0;

  printf("first / last\n");
  for(i = 0; i < g->n; i++) {
    printf("%u / %u\n", g->V[i].first, g->V[i].last);
  }
  printf("\n");

  printf("edge id: (src, tgt, cmp)\n");
  for(i = 0; i < 2*g->m; i++, m++) {
    printf("%u: (%u, %u, %u)\n", m, g->E[i].src, g->E[i].tgt, g->E[i].cmp);
  }
}

void print_tree(Tree* t) {

  if(t->n > 50) {
    fprintf(stderr, "Warning: The number of nodes is greater than 50. We\
  do not recommend to use this function for graphs with more than 50 nodes.\n");
    return;
  }

  uint i, m = 0;

  printf("first / last\n");
  for(i = 0; i < t->n; i++) {
    printf("%u / %u\n", t->N[i].first, t->N[i].last);
  }
  printf("\n");

  printf("edge id: (src, tgt, cmp)\n");
  for(i = 0; i < 2*(t->n-1); i++, m++) {
    printf("%u: (%u, %u, %u)\n", m, t->E[i].src, t->E[i].tgt, t->E[i].cmp);
  }
}

// Size of a graph in bytes
ulong size_graph(Graph* g) {
  ulong sizeV = g->n*sizeof(Vertex);
  ulong sizeE = g->m*sizeof(Edge);
  ulong others = sizeof(g->V)+sizeof(g->E)+sizeof(g->n)+sizeof(g->m);

  return sizeV + sizeE + others;
}

// Size of a tree in bytes
ulong size_tree(Tree* t) {
  ulong sizeN = t->n*sizeof(Node);
  ulong sizeE = 2*(t->n-1)*sizeof(Edge);
  ulong others = sizeof(t->N)+sizeof(t->E)+sizeof(t->n);

  return sizeN + sizeE + others;
}


void write_bits_to_parentheses(const char* fn, BIT_ARRAY* ba) {
  FILE* fp = fopen(fn, "w");

  if (!fp) {
    fprintf(stderr, "Error opening file \"%s\".\n", fn);
    exit(EXIT_FAILURE);
  }

  for(uint i=0; i<ba->num_of_bits; i++) {
    if(bit_array_get_bit(ba, i ) != 0)
      fprintf(fp, "(");
    else
      fprintf(fp, ")");
  }

  fclose(fp);
}
