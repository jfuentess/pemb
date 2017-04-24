/******************************************************************************
 * defs.h
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

#include <math.h>
#include "bitrank/bitrankw32int.h"
#include "succinct_tree.h"

#ifdef NOPARALLEL
#define cilk_for for
#define cilk_spawn
#define cilk_sync
#define __cilkrts_get_nworkers() 1
#else
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/common.h>
#endif

#define threads  __cilkrts_get_nworkers()

/* #define min(a,b)	      \ */
/*   ({ __typeof__ (a) _a = (a); \ */
/*       __typeof__ (b) _b = (b); \ */
/*     _a < _b ? _a : _b; }) */

/* #define max(a,b) \ */
/*   ({ __typeof__ (a) _a = (a); \ */
/*       __typeof__ (b) _b = (b); \ */
/*     _a > _b ? _a : _b; }) */

typedef struct _vertex_t Vertex;
typedef struct _edge_t Edge;
typedef struct _graph_t Graph;
typedef struct _node_t Node;
typedef struct _euler_node_t ENode;
typedef struct _tree_t Tree;
typedef struct _succinct_planar_graph_t SuccGraph;
typedef unsigned int item_t;
typedef struct _stack_t Stack;
typedef struct _queue_t Queue;

struct _stack_t {
  item_t *base;
  uint top;
  uint size;
};

struct _queue_t {
  item_t *base;
  uint front;
  uint rear;
  uint size;
};

// Used to represent graphs
struct _vertex_t {
  uint first; // Position of the first incident edge of a vertex in E
  uint last; // Position of the last incident edge of a vertex in E
};

// Used to represent graphs and trees
struct _edge_t {
  uint src; // Index of the source vertex of the edge
  uint tgt; // Index of the target vertex of the edge
  int cmp; // Position of the complementary edge (in the adjacency list of
	   // tgt). cmp < 0 means that the field cmp is undefined
};

// Used to represent trees
struct _node_t {
  uint first; // Position of the first incident edge of a node in E
  uint last; // Position of the last incident edge of a node in E
};

struct _graph_t {
  Vertex* V; // Array of vertices of the graph
  Edge* E; // Array of edges of the graph. It is the concatenation of the adjacency lists of all vertices
  uint n; // Number of vertices in the graph
  uint m; // Number of non-repeated edges in the graph
};

struct _tree_t {
  Node* N; // Array of nodes of the tree
  Edge* E; // Array of edges of the tree. It is the concatenation of the adjacency lists of all nodes
  uint n; // Number of nodes in the tree
  // The number of edges is n-1
};

// struct for the Euler tour code
struct _euler_node_t {
  int succ; // stores the index of the succesor in the array. Since
	    // the parallel_list_ranking algorithm uses this fields to
	    // store some negative values, it must be int instead of uint.
  char value;
  int rankA; // Used to build the bitmap A
  int rankB; // Used to build the bitmap B
};

struct _succinct_planar_graph_t {
  bitRankW32Int* A;
  rmMt* B;
  rmMt* B_star;
  uint n;
  uint m;
  bitRankW32Int* B_rs; // rank/select structures for B
  bitRankW32Int* B_star_rs; // rank/select structures for B*
};

Graph* create_graph(uint, uint);

void free_graph(Graph*);

Tree* create_tree(uint);

void free_tree(Tree*);

void parallel_prefix_sum(uint*, uint);
void parallel_list_ranking(ENode*, uint);


Stack *stack_create(uint size);
int stack_empty(Stack *);
int stack_full(Stack *);
int stack_push(Stack *, item_t);
item_t stack_pop(Stack *);
item_t stack_top(Stack *);
void stack_free(Stack *);
Stack * halving_stack(Stack *, uint);

Queue *queue_create(uint);
int queue_empty(Queue *);
int queue_full(Queue *);
int queue_enqueue(Queue *, item_t);
item_t queue_dequeue(Queue *);
item_t queue_front(Queue *);
void queue_free(Queue *);

int tree_cycle(Tree *);
void parallel_list_ranking(ENode *, uint);
void connected_graph(Graph *);
