/******************************************************************************
 * planar_graph.c
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
#include "planar_graph.h"

void _mark_edges(Graph *, uint, uint, uint, uint *);

// init: root of the spanning tree
SuccGraph* par_succ_planar_graph(Graph *g, Tree *t, uint init, int *parent, uint
			     *count_edges, uint *references) {
  SuccGraph *sg = malloc(sizeof(SuccGraph));
  ENode* ET = malloc(2*(g->n-1)*sizeof(ENode));

  BIT_ARRAY* A = bit_array_create(2*g->m);
  BIT_ARRAY* B = bit_array_create(2*g->n);
  BIT_ARRAY* B_star = bit_array_create(2*g->m-2*g->n+4);

  uint chk = 2*(t->n-1)/threads;

  cilk_for(uint h = 0; h < threads; h++) {
    uint ll = h*chk;
    uint ul = ll+chk;
    if(h == threads-1)
      ul = 2*(t->n-1);

    for(uint i = ll; i < ul; i++) {
      Edge e = t->E[i];
      Node n_src = t->N[e.src];
      Node n_tgt = t->N[e.tgt];

      ET[i].rankA = 1 + count_edges[t->E[i].cmp];
      ET[i].rankB = 1;

      /* Forward edges. In the root, all the edges are forward edges */
      /* e.src == init means the root
  	 t->E[parent[e.tgt]].cmp != i means a forward edge */
      if((e.src == init) | (parent[e.src] != i)) {
  	ET[i].value = 1;

      /* Leaf */
  	if(n_tgt.first == n_tgt.last)
  	  ET[i].succ = e.cmp;
  	/* Connect to the first child of the tgt node */
  	else {
  	  uint succ = parent[e.tgt]+1;
  	  if(succ > t->N[e.tgt].last)
  	    succ = t->N[e.tgt].first;

  	  ET[i].succ = succ;
  	}
      }
      /* Backward edges. The root has not backward edges */
      else {
  	ET[i].value = 0;
  	/* Especial case for the last child of the root */
  	/* Assumption: N[init].first in an external edge and it is part of the
  	   spanning tree */
  	if(e.tgt == init && n_tgt.last == e.cmp)
  	  ET[i].succ = 0;
  	else {
  	  uint last = parent[e.tgt] - 1;
  	  if(last < t->N[e.tgt].first)
  	    last = t->N[e.tgt].last;
	  
  	  /* Last child of t->N[i] */
  	  if(e.cmp == last)
  	    ET[i].succ = parent[e.tgt];
  	  /* Intermediate child */
  	  else {
  	    uint next = e.cmp + 1;
  	    if(next > t->N[e.tgt].last)
  	      next = t->N[e.tgt].first;
  	    ET[i].succ = next;

  	  }
  	}
      }
    }
  }

  parallel_list_ranking(ET, 2*(t->n-1));
  
  cilk_for(uint h = 0; h < threads; h++) {
    uint ll = h*chk;
    uint ul = ll+chk;
    if(h == threads-1)
      ul = 2*(t->n-1);
    
    for(uint i = ll; i < ul; i++) {
      parallel_or_bit_array_set_bit(A, ET[i].rankA);
      if(ET[i].value) {
	// TODO: Eliminate the addition (+1). Include that work into the
	// computation of the parallel_list_ranking
  	parallel_or_bit_array_set_bit(B, ET[i].rankB+1);
      }
      
    }
  }
  parallel_or_bit_array_set_bit(B, 0);

  /*
   * To construct the bitvector B*, we need to order the edges of G\T according
   * to the traversal of T. Using the implicit order of B* in ET (rankA-rankB), we
   * write in marked_edges the indices in G of the edges that are part of
   * B*. The order of the indices in marked_edges will be the final order in
   * B*. Additionally, we reuse the field 'src' of the edges in G\T to store its
   * position in marked_edges. It will be used later.
   */
  uint *marked_edges = malloc((2*g->m-2*g->n+2)*sizeof(uint));
      
  cilk_for(uint h = 0; h < threads; h++) {
    uint ll = h*chk;
    uint ul = ll+chk;
    if(h == threads-1)
      ul = 2*(t->n-1);
    
    for(uint i = ll; i < ul; i++)
      if(count_edges[t->E[i].cmp])
  	_mark_edges(g, references[t->E[i].cmp], count_edges[t->E[i].cmp],
  		    ET[i].rankA-ET[i].rankB, marked_edges);
  }

  free(ET);

  /*
   * Finally, to decide if a bit in B* must be 0 o 1, we need to check the field
   * 'src' of the corresponding edge in G\T and its complement. Let e be the
   * current edge that we are checking and let e' be the complement edge of
   * e. If e.src (the position in marked_edges and the final position in B*) is
   * less than e'.src, then e is a forward edge. Therefore, we write a
   * 0. Otherwise, we write a 1.
   */
  chk = (2*g->m-2*g->n+2)/threads;
  cilk_for(uint h = 0; h < threads; h++) {
    uint ll = h*chk;
    uint ul = ll+chk;
    if(h == threads-1)
      ul = (2*g->m-2*g->n+2);
    
    for(uint i = ll; i < ul; i++) {
      Edge e = g->E[marked_edges[i]];
      if(e.src <= g->E[e.cmp].src)
	// TODO: Eliminate the addition (+1). Include that work into the
	// computation of the parallel_list_ranking
  	parallel_or_bit_array_set_bit(B_star, i+1);
    }
  }
  parallel_or_bit_array_set_bit(B_star, 0);

  free(marked_edges);

  // bitarray A is not a balanced parentheses sequence, so, we cannot use the
  // rmMt
  sg->A = createBitRankW32Int(A->words, A->num_of_bits, 1, 20); // Factor 20
  sg->B = st_create(B, B->num_of_bits);
  sg->B_star = st_create(B_star, B_star->num_of_bits);
  sg->B_rs = createBitRankW32Int(B->words, B->num_of_bits, 1, 20); // Factor 20

  sg->n = g->n;
  sg->m = g->m;
  return sg;
}


/* 
   This function store the indices of the edges in G\T in the final order of
   B*. Besides, it reuses the field 'src' of the edges to store the final
   position of that edge in B*. In each call, the function visits as many nodes
   as edges of G\T between two consecutive edges of T (starting from the edge of
   G in ref).
   
   To Do: The number of consecutives edges can be greater than O(n/p). Break the
   range in parallel tasks of at most O(n/p) edges.
*/
void _mark_edges(Graph *g, uint ref, uint cnt, uint pos, uint *marked_edges) {
  uint limit = ref + cnt;
  Edge e = g->E[ref];
  if(limit <= g->V[e.src].last)
    for(uint i = ref+1; i <= limit; i++) {
      g->E[i].src = pos;
      marked_edges[pos++] = i;
    }
  else {
    for(uint i = ref+1; i <= g->V[e.src].last; i++) {
      g->E[i].src = pos;
      marked_edges[pos++] = i;
    }

    limit = g->V[e.src].first - g->V[e.src].last + limit;

    for(uint i = g->V[e.src].first; i < limit; i++) {
      g->E[i].src = pos;
      marked_edges[pos++] = i;
    }
  }
}


// init: root of the spanning tree
SuccGraph* seq_succ_planar_graph(Graph *g, Tree *t, uint init, int *parent, uint
			     *count_edges, uint *references) {
  SuccGraph *sg = malloc(sizeof(SuccGraph));

  Stack *s = stack_create(2*t->n);
  stack_push(s, t->N[init].first);
  char* visited = calloc(2*t->n, sizeof(char));//edges
  int par, curr = -1;
  
  BIT_ARRAY* A = bit_array_create(2*g->m);
  BIT_ARRAY* B = bit_array_create(2*g->n);
  BIT_ARRAY* B_star = bit_array_create(2*g->m-2*g->n+4);

  uint *marked_edges = malloc((2*g->m-2*g->n+2)*sizeof(uint));

  uint idx = 0;
  uint ii = 0;
  uint pos = 0;
  for(uint i = 1; i <= 2*(t->n-1); i++) { 
    Edge e = t->E[idx];
    Node n_src = t->N[e.src];
    Node n_tgt = t->N[e.tgt];


    // precompute forward and backward edges for B*
    uint ref = references[t->E[idx].cmp];
    uint limit = ref + count_edges[t->E[idx].cmp];
    Edge ee = g->E[ref];
    if(limit <= g->V[ee.src].last)
      for(uint j = ref+1; j <= limit; j++) {
	g->E[j].src = pos;
	marked_edges[pos++] = j;
      }
    else {
      for(uint j = ref+1; j <= g->V[ee.src].last; j++) {
	g->E[j].src = pos;
	marked_edges[pos++] = j;
      }
      
      limit = g->V[ee.src].first - g->V[ee.src].last + limit;
      
      for(uint j = g->V[ee.src].first; j < limit; j++) {
	g->E[j].src = pos;
	marked_edges[pos++] = j;
      }
    }
        
    bit_array_set_bit(A, ii);
    ii += count_edges[t->E[idx].cmp]+1;

    /* Forward edges. In the root, all the edges are forward edges */
    /* e.src == init means the root
       t->E[parent[e.tgt]].cmp != i means a forward edge */
    if((e.src == init) | (parent[e.src] != idx)) {
      bit_array_set_bit(B, i);

      /* Leaf */
      if(n_tgt.first == n_tgt.last)
	idx = e.cmp;
      /* Connect to the first child of the tgt node */
      else {
	uint succ = parent[e.tgt]+1;
	if(succ > t->N[e.tgt].last)
	  succ = t->N[e.tgt].first;

	idx = succ;
      }
    }
    else {
      /* Especial case for the last child of the root */
      /* Assumption: N[init].first in an external edge and it is part of the
	 spanning tree */
      if(e.tgt == init && n_tgt.last == e.cmp)
	break;
      else {
	uint last = parent[e.tgt] - 1;
	if(last < t->N[e.tgt].first)
	  last = t->N[e.tgt].last;
	
	/* Last child of t->N[i] */
	if(e.cmp == last)
	  idx = parent[e.tgt];
	/* Intermediate child */
	else {
	  uint next = e.cmp + 1;
	  if(next > t->N[e.tgt].last)
	    next = t->N[e.tgt].first;
	  idx = next;	  
	}
      }
    }
  }
  bit_array_set_bit(B, 0);


  uint ul = 2*g->m-2*g->n+2;
  for(uint i = 0; i < ul; i++) {
    Edge e = g->E[marked_edges[i]];
    if(e.src <= g->E[e.cmp].src)
      bit_array_set_bit(B_star, i+1);
  }
  
  bit_array_set_bit(B_star, 0);

  free(marked_edges);

  sg->A = createBitRankW32Int(A->words, A->num_of_bits, 1, 20); // Factor 20
  sg->B = st_create(B, B->num_of_bits);
  sg->B_star = st_create(B_star, B_star->num_of_bits);
  sg->B_rs = createBitRankW32Int(B->words, B->num_of_bits, 1, 20); // Factor 20
  sg->B_star_rs = createBitRankW32Int(B_star->words, B_star->num_of_bits, 1, 20); // Factor 20

  sg->n = g->n;
  sg->m = g->m;
  
  return sg;
}


/* Assuming indices start with 0 */
uint first(SuccGraph *sg, uint v) {
  if(v >= 0) {
    uint pos = select1(sg->B_rs, v+1);
    uint edge = select1(sg->A, pos);
    if(v == 0) // The root of the spanning tree
      return edge;
    else
      return edge+1;
  } else
    return -1;
}

/* Assuming indices start with 0 */
uint mate(SuccGraph *sg, uint i) {
  if(isBitSet(sg->A, i) == 1) {
    uint pos_in_B = rank(sg->A, i); // rank1
    uint match_in_B = match(sg->B, pos_in_B);
    return select1(sg->A, match_in_B);
  }
  else
										 {
    uint pos_in_B_star = i - rank(sg->A, i) + 1; // rank0
    uint match_in_B_star = match(sg->B_star, pos_in_B_star);
    return select0(sg->A, match_in_B_star);
  }  
  return -1;
}

/* Assuming indices start with 0 */
uint next(SuccGraph *sg, uint i) {
  if(i > sg->A->n)
    return -1;
  
  if(isBitSet(sg->A, i) == 0) {
    return i+1;
  }
  else {
    uint pos_in_B = rank(sg->A, i); // rank1
    if(bit_array_get_bit(sg->B->bit_array, pos_in_B) == 1) {
      return mate(sg, i) + 1;
    }
  }
  return -1;
}

// Size of a succinct planar graph in bytes
ulong size_planar_graph(SuccGraph *sg) {
  ulong sizeSG= sizeof(SuccGraph);
  ulong sizeA = length_in_bits(sg->A)/8;
  ulong sizeB = size_rmMt(sg->B);
  ulong sizeB_star = size_rmMt(sg->B_star);
  
  return sizeSG + sizeA + sizeB + sizeB_star;
}

uint degree(SuccGraph *sg, uint v) {
  if(v >= sg->n)
    return 0;
  
  uint dg = 0;
  uint nxt = first(sg, v);
  
  while(nxt < 2*sg->m) {
    nxt = next(sg, nxt);
    dg++;
  }
  
  return dg;
}

uint vertex(SuccGraph *sg, uint v) {
  uint pos_in_A = rank(sg->A, v); // rank1
  if(isBitSet(sg->A, v)==1) {
    if(bit_array_get_bit(sg->B->bit_array, pos_in_A+1) == 0) {
      uint match_pos = match(sg->B, pos_in_A+1);
      return (rank(sg->B_rs, match_pos) - 1);
    }
    else {
      uint rank1_B = rank(sg->B_rs, pos_in_A+1);
      return parent_t(sg->B, rank1_B-1);
    }
  }
  else {
    if(bit_array_get_bit(sg->B->bit_array, pos_in_A+1) == 0) {
      return rank(sg->B_rs, pos_in_A)-1;
    }
    else {
      uint match_pos = match(sg->B, pos_in_A+1);
      uint rank1_B = rank(sg->B_rs, match_pos);
      return parent_t(sg->B, rank1_B-1);
    }
  }
}

void list_neighbors(SuccGraph *sg, uint v) {
  if(v >= sg->n)
    return;
  
  uint nxt = first(sg, v);
  while(nxt < 2*sg->m) {
    if(nxt < 2*sg->m) {
      uint mt = mate(sg, nxt);
      uint x = vertex(sg, mt);
    }
    nxt = next(sg, nxt);	
  }
}

int balanced(BIT_ARRAY* ba) {
  Stack* s = stack_create(ba->num_of_bits);

  for(uint i=0; i<ba->num_of_bits; i++) {
    if(bit_array_get_bit(ba, i ) != 0) {
      stack_push(s, 1);
    }
    else
      stack_pop(s);
  }

  if(stack_empty(s))
    return 1;
    
  return 0;
}
