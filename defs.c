/******************************************************************************
 * defs.c
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
#include "defs.h"
#include <math.h>

Graph* create_graph(uint n, uint m) {

  Graph *g = malloc(sizeof(Graph));
  g->n = n;
  g->m = m;
  g->V = malloc(g->n*sizeof(Vertex));
  g->E = malloc(2*(g->m)*sizeof(Edge));

  return g;
}

void free_graph(Graph* g) {
  free(g->V);
  free(g->E);
  free(g);
}

Tree* create_tree(uint n) {

  Tree *t = malloc(sizeof(Tree));
  t->n = n;
  t->N = malloc(t->n*sizeof(Node));
  t->E = malloc(2*(t->n-1)*sizeof(Edge));

  return t;
}

void free_tree(Tree* t) {
  free(t->N);
  free(t->E);
  free(t);
}

/*
Compute in parallel the prefix sum of an array of uints
Input: An array A of uints, the size the array and the number of available threads.
Output: None. The prefix sums will be saved in the array A

Note: This algorithm assumes that size >= threads
*/
void parallel_prefix_sum(uint* A, uint size) {
  uint chk = ceil((float)size/threads);
  cilk_for(uint i = 0; i < threads; i++) {
    uint  ll = i*chk, ul = ll + chk;
    if(ul > size)
      ul = size;

    uint acc = 0;
    for(uint j = ll; j < ul; j++) {
      A[j] += acc;
      acc = A[j];
    }
  }
  
  for(uint i = 1; i < threads-1; i++)
    A[((i+1)*chk)-1] += A[i*chk-1];
  
  if(threads > 1)
    A[size-1] += A[(threads-1)*chk-1];
  
  cilk_for(uint i = 1; i < threads; i++) {
    uint ll = i*chk, ul = ll + chk - 1;
    if(ul >= size)
      ul = size - 1;

    uint acc = A[ll-1];
    for(uint j = ll; j < ul; j++) {
      A[j] += acc;
    }
  }
}


/*
 * Stack implementation. The implementation is based on the code in Chapter 1.1
 * of "Advanced Data Structures" by Peter Brass.
 *
 * This implementation is good in terms of performance, but not in terms of
 * memory usage. However, in this project we prefer a good performance than a
 * memory usage.
 */

Stack *stack_create(uint size) {
  Stack *s;
  s = (Stack *) malloc(sizeof(Stack));
  s->base = (item_t *) malloc(size * sizeof(item_t));
  s->size = size;
  s->top = 0;

  return s;
}

int stack_empty(Stack *s) {
  return (s->top == 0);
}

int stack_full(Stack *s) {
  return (s->top == s->size);
}

int stack_push(Stack *s, item_t x) {
  if(!stack_full(s)) {
    s->base[s->top] = x;
    s->top += 1;

    if(s->top > s->size)
      printf("push - top: %d\n", s->top);

    return 1;
  }
  else {
    fprintf(stderr, "Error in function stack_push: The stack is full\n");
    exit(EXIT_FAILURE);
  }
}

item_t stack_pop(Stack *s) {
  if(!stack_empty(s)) {
    s->top -= 1;
    if(s->top > s->size)
      printf("pop - top: %d\n", s->top);
    return s->base[s->top];
  }
  else {
    fprintf(stderr, "Error in function stack_pop: The stack is empty\n");
    exit(EXIT_FAILURE);
  }
}

item_t stack_top(Stack *s) {
  if(!stack_empty(s))
    return s->base[s->top];
  else {
    fprintf(stderr, "Error in function stack_top: The stack is empty\n");
    exit(EXIT_FAILURE);
  }
}

void stack_free(Stack *s) {
  free(s->base);
  free(s);
}

Stack * halving_stack(Stack *s, uint size) {
    
  Stack *s2 = stack_create(size);
  uint half = size/2;

  for(int i=0; i<half; i++)
    stack_push(s2, stack_pop(s));

  return s2;
}

/*
 * Queue implementation. The implementation is based on the code in Chapter 1.2
 * of "Advanced Data Structures" by Peter Brass.
 *
 * This implementation is good in terms of performance, but not in terms of
 * memory usage. However, in this project we prefer a good performance than a
 * memory usage.
 * 
 * Drawback of this implementation: If the input size is s, then the maximum
 * number of elements is s-1. It can be solved be using int variables for rear
 * and front, but in that case the maximum number of elements is reduced
 * (2^31-1)
 */
Queue *queue_create(uint size) {
  Queue *q;
  q = (Queue *) malloc(sizeof(Queue));
  q->base = (item_t *) malloc(size * sizeof(item_t));
  q->size = size;
  q->front = 0;
  q->rear = 0;

  return q;
}

int queue_empty(Queue *q) {
  return (q->front == q->rear);
}

int queue_full(Queue *q) {
  if(q->front == ((q->rear+1)%q->size))
    return 1;

  return 0;
}

int queue_enqueue(Queue *q, item_t x) {
  if(!queue_full(q)) {
    q->base[q->rear] = x;
    q->rear = (q->rear+1)%q->size;
    return 1;
  }
  else {
    fprintf(stderr, "Error in function queue_enqueue: The queue is full\n");
    exit(EXIT_FAILURE);
  }
}

item_t queue_dequeue(Queue *q) {
  if(!queue_empty(q)) {
    uint tmp = q->front;
    q->front = (q->front+1)%q->size;
    return q->base[tmp];
  }
  else {
    fprintf(stderr, "Error in function queue_dequeue: The queue is empty\n");
    exit(EXIT_FAILURE);
  }
}

item_t queue_front(Queue *q) {
  if(!queue_empty(q))
    return q->base[q->front];
  else {
    fprintf(stderr, "Error in function queue_front: The queue is empty\n");
    exit(EXIT_FAILURE);
  }
}

void queue_free(Queue *q) {
  free(q->base);
  free(q);
}


/*
 * basic algorithms. In this section we will include the implementation of
 * some algorithms needed for our parallel planar graph algorithm. The
 * implementations incluude:
 *   - Tree cycle: Verify if the input tree is actually a tree, by finding
 *   cycles.
 *   - Parallel prefix sum
 *   - Parallel list ranking
 */

int tree_cycle(Tree* t) {
  Stack *s = stack_create(t->n); 
  char* visited = calloc(t->n, sizeof(char));
  char* edges = calloc(2*(t->n-1), sizeof(char));
  int par, curr = -1;
  uint edge;
  int first = 1;

  while(!stack_empty(s) || first) {
    if(first) { // Root
      par = -1;
      curr = 0;
      first = 0;
    }
    else {
      edge = stack_pop(s);
      par = t->E[edge].src; // parent
      curr = t->E[edge].tgt; // current
    }
    visited[curr] = 1;
    for(uint i = t->N[curr].first; i <= t->N[curr].last; i++) {
      edges[i] = 1;
      if(!visited[t->E[i].tgt])
	stack_push(s, i);      
      else if(par != t->E[i].tgt && par != -1) {
	printf("\tThe tree has cycles\n");
	return 1;
      }
    }
  }

  for(uint i = 0; i < t->n; i++)
    if(visited[i] == 0) { // There are unvisited nodes
      printf("\tNo connected tree (unvisited nodes)\n");
      return 1;
    }

  for(uint i = 0; i < 2*(t->n-1); i++)
    if(edges[i] == 0) { // There are unvisited edges
      printf("\tNo connected graph (unvisited edges)\n");
      return 1;
    }
  
  return 0;
}

/* 
 * Assuming that the node in position 0 is the head of the list
 * This implementation assumes that the list A has two different values (valueA
 * and valueB) that need to be added.
 */
void parallel_list_ranking(ENode* A, uint size) {

  struct sublist_node {
    int head;
    int next;
    int scratch;
    int valueA;
    int valueB;
  };
  
  uint s = ceil(log2(size)*threads);
  uint chk = size/s;
  if(s > size) {
    s = size;
    chk = 1;     
  }
  /* printf("size: %u, s: %u, chk: %u\n", size, s, chk); */

  struct sublist_node* sublist = malloc(s*sizeof(struct sublist_node));
  
  // Compute the splitters
  cilk_for(uint i = 0; i < s; i++) {
    uint x = i*chk;
    sublist[i].head = x;
    sublist[i].valueA = A[x].rankA;
    sublist[i].valueB = A[x].rankB;
    sublist[i].next = -1;
    sublist[i].scratch = A[x].succ;
    A[x].succ = -(i)-1;
    /* printf("\t%u: %u\n", x, -(A[x].succ)); */
    /* printf("\t[*]%u: %u\n", x, sublist[i].scratch); */
  }
  
  cilk_for(uint i = 0; i < s; i++) {
    int curr = sublist[i].scratch;
    uint tmpA = 0, tmpB = 0, tmp2 = 0;
    
    while(curr > 0) {
      tmp2 = A[curr].rankA;
      A[curr].rankA = tmpA;
      tmpA += tmp2;
      
      tmp2 = A[curr].rankB;
      A[curr].rankB = tmpB;
      tmpB += tmp2;
      
      int aux = A[curr].succ;
      A[curr].succ = -(i)-1;
      curr = aux;
    }
    sublist[i].next = -(curr)-1;
    
    // Special case
    if(curr != 0) {
      sublist[-(curr)-1].valueA = tmpA;
      sublist[-(curr)-1].valueB = tmpB;
    }
  }
  
  int curr = 0;
  int tmpA = 0, tmpB = 0, tmp2 = 0;
  
  while(1) {
    tmp2 = sublist[curr].valueA;
    sublist[curr].valueA += tmpA;
    tmpA += tmp2;

    tmp2 = sublist[curr].valueB;
    sublist[curr].valueB += tmpB;
    tmpB += tmp2;

    curr = sublist[curr].next;
    if(curr < 0)
      break;
  }
  
  cilk_for(uint i = 0; i < s; i++) {
    uint  ll = i*chk, ul = ll + chk;
    if(i == s-1)
      ul = size;
    if(i == 0)
      ll++;

    /* fprintf(stderr, "ll: %u, ul: %u\n", ll, ul); */
    for(uint j = ll; j < ul; j++) {
       int idx = -(A[j].succ)-1;
       A[j].rankA += sublist[idx].valueA;
       A[j].rankB += sublist[idx].valueB;

       /* if(A[j].rankA < 0) */
       /* 	 printf("j:%u, idx: %d, valueA:%d, valueB:%d\n", j, idx, sublist[idx].valueA, */
       /* 		sublist[idx].valueB); */
    }
    // Assuming that the field "succ" will not be longer used
    //    A[i].succ = sublist[i].scratch;
  }
  
  A[0].rankA=0;
  A[0].rankB=0;
  free(sublist);
}

// Verify is the input graph is a connected graph
void connected_graph(Graph* g) {
  Stack *s = stack_create(2*g->n); 
  char* visited = calloc(g->n, sizeof(char));
  int par, curr = -1;
  uint edge;
  int first = 1;
  int components = 1;
  int all_visited = -1;
  int num_vertices = 0;

  while(!stack_empty(s) || first) {
    if(first) { // Root
      par = -1;
      curr = 0;
      first = 0;
    }
    else {
      edge = stack_pop(s);
      par = g->E[edge].src; // parent
      curr = g->E[edge].tgt; // current
    }
    visited[curr] = 1;
    for(uint i = g->V[curr].first; i <= g->V[curr].last; i++) {
      if(!visited[g->E[i].tgt])
	stack_push(s, i);      
    }
  }
  
  printf("num_vertices: %d\n", num_vertices);
  for(uint i = 0; i < g->n; i++)
    if(visited[i] == 0) { // There are unvisited vertices
      num_vertices++;
    }


  printf("unvisited vertices: %d, visited vertices: %d\n", num_vertices, g->n-num_vertices);
}

