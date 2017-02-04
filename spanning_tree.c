/******************************************************************************
 * spanning_tree.h
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

#include "spanning_tree.h"

uint global_cnt = 0;
void _stub_spanning_tree(Graph *, uint, uint, BIT_ARRAY *, int *, Stack**,
			 uint); 
void _local_spanning_tree(Graph *, BIT_ARRAY *, int *, Stack*);
void _local_spanning_tree2(Graph *, BIT_ARRAY *, int *, Stack*, uint);

Tree *seq_dfs_spanning_tree(Graph *g, uint init, int *parent, uint *count_edges,
  uint * references) {
  Tree *t = create_tree(g->n);
  BIT_ARRAY *visited = bit_array_create(g->n);
  uint *edges = (uint *)calloc(2*g->m,sizeof(uint));
  uint num_tree_edges = 2*(t->n-1);
  
  Stack *s = stack_create(g->n);
  bit_array_set_bit(visited, init);
  stack_push(s, init);
  parent[init] = -1;
  
  while(!stack_empty(s)) {
    uint curr = (uint)stack_pop(s);

    for(uint i = g->V[curr].first; i <= g->V[curr].last; i++) {
      uint tgt = g->E[i].tgt;
      
      if(!bit_array_get_bit(visited, tgt)) { // Not visited	
	bit_array_set_bit(visited, tgt);
	stack_push(s, tgt);
	parent[tgt] = g->E[i].cmp; // Edge child-to-parent
      }
    }
  }

  // Marking the edges of G that are in T
  for(uint i = 0; i < init; i++) {
    edges[(uint)parent[i]] = 1;
    edges[(uint)g->E[parent[i]].cmp] = 1;
  }
  for(uint i = init+1; i < g->n; i++) {
    edges[(uint)parent[i]] = 1;
    edges[(uint)g->E[parent[i]].cmp] = 1;
  }

  uint m = 0;
  /* Counting edges */
  for(uint i = 0; i < g->n; i++) {
    uint first = g->V[i].first;
    uint last = g->V[i].last;
    uint zeros = 0;
    int last_one = -1;
    uint carry_last = 0;
    
    for(uint j = first; j <= last; j++) {   
      if(edges[j] == 1) {
	t->E[m].src = g->E[j].src;
	t->E[m].tgt = g->E[j].tgt;

	if(last_one == -1)
	  carry_last = zeros;	
	else
	  count_edges[last_one] += zeros;

	references[m] = j;
	last_one = m;
	zeros = 0;

	if(m == 0)
	  t->N[t->E[m].src].first = m;
	else if(t->E[m].src != t->E[m-1].src) {	
	  t->N[t->E[m-1].src].last = m-1;
	  t->N[t->E[m].src].first = m;
	  
	}
	m++;
      }
      else
	zeros++;
      edges[j] = m;
    }

    count_edges[last_one] += zeros + carry_last;
    t->N[t->E[m-1].src].last = m-1;
  }

  for(uint i = 0; i < init; i++)
      parent[i] = edges[parent[i]]-1;
  for(uint i = init+1; i < g->n; i++)
      parent[i] = edges[parent[i]]-1;
  
  for(uint i = 0; i < num_tree_edges; i++) {
    Node tgt = t->N[t->E[i].tgt];
    
    for(uint j = tgt.first; j <= tgt.last; j++)
      if(t->E[j].tgt == t->E[i].src) {
  	t->E[i].cmp = j;
  	break;
      }
  }

  free(visited);
  free(edges);
  stack_free(s);
    
  return t;
}

Tree* par_spanning_tree2(Graph* g, uint init, int *parent, uint *count_edges,
			uint * references) {

  Tree *t = create_tree(g->n);

  BIT_ARRAY *visited = bit_array_create(g->n);
  uint *edges = (uint *)calloc(2*g->m,sizeof(uint));
  uint size = 10*threads; /* Maximum number of elements in the stack of the stub
			  spanning tree */ 
  
  // One stack per thread
  Stack **s_threads = (Stack **) malloc(threads * sizeof(Stack *));

  // The spawn depth (number of tasks/frames that a thread can manage) is
  // 1024. This number is fixed and defined in
  // libcilkrts/runtime/global_state.cpp. Use s_size to control the spawn depth
  uint s_size = g->n/(100*threads); // Size of the stack of each thread

  parent[init] = -1;
  
  _stub_spanning_tree(g, init, size, visited, parent, s_threads, s_size);

  cilk_for(uint i = 0; i < threads; i++) {
    _local_spanning_tree2(g, visited, parent, s_threads[i], s_size);
  }
       
  uint chk = ceil((double)g->n/threads);
  cilk_for(uint i = 0; i < threads; i++) {
    uint ll = i*chk;
    uint ul = (i+1)*chk;
    if(ul > g->n)
      ul = g->n;

    /* Make it faster: Thread-safe only at the beginning and the end*/
    for(uint j = ll; j < ul; j++) {
      if(init != j) {
  	edges[(uint)parent[j]] = 1;
  	edges[(uint)g->E[parent[j]].cmp] = 1;
      }
    }
  }
    
  // Prefix Sum
  uint *offsets = calloc(threads, sizeof(uint));
  chk = ceil((double)2*g->m/threads);
  cilk_for(uint i = 0; i < threads; i++) {
    uint  ll = i*chk, ul = ll + chk;
    if(ul > 2*g->m)
      ul = 2*g->m;

    uint acc = 0;
    for(uint j = ll; j < ul; j++)
      acc += edges[j];
    offsets[i] = acc;
  }

  uint tmp = 0;
  uint acc = 0;
  for(uint i = 0; i < threads; i++) {
    acc += offsets[i];
    offsets[i] = tmp;
    tmp = acc;
  }

  // Especial case for the complement of the first edge
  uint first_edge = 0;
  for(uint i = g->V[init].first; i < g->V[init].last; i++)
    if(edges[i]) {
      first_edge = i;
      break;
    }
	
  cilk_for(uint i = 0; i < threads; i++) {
    uint ll = i*chk, ul = ll + chk;
    if(ul >= 2*g->m)
      ul = 2*g->m;

    uint local_offset = offsets[i];
    for(uint j = ll; j < ul; j++) {
      if(edges[j]) {
  	t->E[local_offset].src = g->E[j].src;
  	t->E[local_offset].tgt = g->E[j].tgt;
  	edges[j] = local_offset;
  	local_offset++;
      }
    }
  }
 
  cilk_for(uint i = 0; i < g->n; i++) {
    uint first = g->V[i].first;
    uint last = g->V[i].last;
    if(i == 0)
      first++;
    uint zeros = 0;
    int last_one = -1;
    uint carry_last = 0;

    // NOTE: Assuming that, for each vertex v of G, degree(v) < O(n/p)
    // TODO: Implement the case degree(v) > O(n/p)
    for(uint j = first; j <= last; j++) {
      if(edges[j]) {
  	if(last_one == -1)
  	  carry_last = zeros;
  	else
  	  count_edges[last_one] += zeros;

  	references[edges[j]] = j;
  	last_one = edges[j];
  	zeros = 0;
      }
      else
  	zeros++;
    }

    count_edges[last_one] += zeros + carry_last;
  }
    
  uint cmp_first_edge = 0;
  cilk_for(uint i = 0; i < threads; i++) {
    uint ll = i*chk, ul = ll + chk;
    if(ul >= 2*g->m)
      ul = 2*g->m;

    uint local_offset = offsets[i];
    for(uint j = ll; j < ul; j++) {
      if(edges[j]) {
  	// Especial case for the complement of the first edge
  	if(g->E[j].cmp == first_edge)
  	  cmp_first_edge = j;
  	t->E[edges[j]].cmp = edges[g->E[j].cmp];
      }
    }
  }

  // Especial case for the complement of the first edge
  t->E[0].cmp = edges[cmp_first_edge];
  references[0] = first_edge;
    
  chk = ceil((double)g->n/threads);
  cilk_for(uint i = 0; i < threads; i++) {
    uint ll = i*chk, ul = ll + chk;
    if(ul >= g->n)
      ul = g->n;
 
    for(uint j = ll; j < ul; j++) {
      uint new_first = 0;
      uint new_last = 0;
      uint idx = 0;
      
      // NOTE: Assuming that, for each vertex v of G, degree(v) < O(n/p)
      // TODO: Implement the case degree(v) > O(n/p)
      idx = g->V[j].first;
      while(!edges[idx])
  	idx++;
      t->N[j].first = edges[idx];
      
      idx = g->V[j].last;
      while(!edges[idx])
  	idx--;
      t->N[j].last = edges[idx];
    }
  }

  t->N[0].first = 0;  

  chk = ceil((double)g->n/threads);
  cilk_for(uint i = 0; i < threads; i++) {
    uint ll = i*chk, ul = ll + chk;
    if(ul > g->n)
      ul = g->n;
    for(uint j = ll; j < ul; j++) {
      if(parent[j] >= 0)
  	parent[j] = edges[parent[j]];
    }
  }

  return t;
}


void _stub_spanning_tree(Graph* g, uint init, uint size, BIT_ARRAY *visited,
			 int *parent, Stack **s_threads, uint s_size) {

  Stack *s = stack_create(size);

  bit_array_set_bit(visited, init);
  stack_push(s, init);
  
  uint cnt = 1;
  uint ele_stack = ceil(size/threads); // Elements per stack
  
  while((cnt < size) && !stack_empty(s)) {
    uint curr = (uint)stack_pop(s);
    cnt--;
    for(uint i = g->V[curr].first; cnt < size && i <= g->V[curr].last; i++) {
      uint tgt = g->E[i].tgt;
      
      if(!bit_array_get_bit(visited, tgt)) { // Not visited
	bit_array_set_bit(visited, tgt);
	stack_push(s, tgt);
	cnt++;
	parent[tgt] = g->E[i].cmp; // Edge child-to-parent
      }
    }
  }
    
  for(uint i = 0; i < threads; i++) {
    s_threads[i] = stack_create(s_size);
    for(uint j = 0; j < ele_stack && cnt > 0; j++, cnt--)
      stack_push(s_threads[i], stack_pop(s));
  }

  return;
}

void _local_spanning_tree2(Graph* g, BIT_ARRAY *visited, int *parent, Stack
			   *s, uint s_size) {
  while(!stack_empty(s)) {
    uint curr = (uint)stack_pop(s);

    for(uint i = g->V[curr].first; i <= g->V[curr].last; i++) {
      uint src = g->E[i].src;
      uint tgt = g->E[i].tgt;
      
      if(!bit_array_get_bit(visited, tgt)) { // Not visited
	if(stack_full(s)) {
	  Stack *s2 = halving_stack(s, s_size);
	  cilk_spawn _local_spanning_tree2(g, visited, parent, s2, s_size);
	}
	  
	stack_push(s, tgt);

	parent[tgt] = g->E[i].cmp; // Edge child-to-parent
	if(tgt >= visited->num_of_bits || tgt < 0)
	  printf("ERROR = src: %d, tgt:%d, num_bits:%d, node:%d, edge:%d\n",
		 src, tgt, visited->num_of_bits, curr, i);
	parallel_or_bit_array_set_bit(visited, tgt);
      }
    }
  }
  return;
}
