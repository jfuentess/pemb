/******************************************************************************
 * succinct_tree.c
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

#include "lookup_tables.h"
#include "binary_trees.h"
#include "succinct_tree.h"
#include "bit_array.h"
#include "util.h"
#include "bitrank/basic.h"

/* ASSUMPTIONS:
 * - s = 256 (8 bits) (Following the sdsl/libcds implementations)
 * - k = 2 (Min-max tree will be a binary tree)
 * - Each thread has to process at least one chunk with parentheses (Problem with n <= s)
 */

rmMt* init_rmMt(unsigned long n) {
  rmMt* st = (rmMt*)malloc(sizeof(rmMt));
  st->s = 256;
  st->k = 2;
  st->n = n;
  st->num_chunks = ceil((double)n/st->s);
  st->height = ceil(log(st->num_chunks)/log(st->k)); // heigh = logk(num_chunks), Heigh of the min-max tree
  st->internal_nodes = (pow(st->k,st->height)-1)/(st->k-1); // Number of internal nodes;

  return st;
}

void print_rmMt(rmMt* st) {
  fprintf(stderr, "Chunk size: %u\n", st->s);
  fprintf(stderr, "Arity: %u\n", st->k);
  fprintf(stderr, "Number of parentheses: %lu\n", st->n);
  fprintf(stderr, "Number of chunks (leaves): %u\n", st->num_chunks);
  fprintf(stderr, "Height: %u\n", st->height);
  fprintf(stderr, "Number of internal nodes: %u\n", st->internal_nodes);
}

rmMt* st_create(BIT_ARRAY* bit_array, unsigned long n) {
  rmMt* st = init_rmMt(n);
  /* print_rmMt(st); */

  st->e_prime = (depth_t*)calloc(st->num_chunks,sizeof(depth_t));
  // num_chunks leaves plus internal nodes
  st->m_prime = (depth_t*)calloc(st->num_chunks + st->internal_nodes,sizeof(depth_t));
  // num_chunks leaves plus internal nodes
  st->M_prime = (depth_t*)calloc(st->num_chunks + st->internal_nodes,sizeof(depth_t));  
  st->n_prime = (int16_t*)calloc(st->num_chunks + st->internal_nodes,sizeof(int16_t));
  st->bit_array = bit_array;
  
  if(st->s >= n){
    fprintf(stderr, "Error: Input size is smaller or equal than the chunk size (input size: %lu, chunk size: %u)\n", n, st->s);
    exit(0);
  }
  
  /*
   * STEP 2: Computation of arrays e', m', M' and n'
   */
  unsigned int num_threads;
  if(st->num_chunks < threads)
    num_threads = st->num_chunks;
  else
    num_threads = threads;

  // Each thread works on 'chunks_per_thread' consecutive chunks of the bit_array 
  unsigned int chunks_per_thread = ceil((double)st->num_chunks/num_threads);

  /*
   * STEP 2.1: Each thread computes the prefix computation in a range of the bit array
   */

  cilk_for(unsigned int thread = 0; thread < num_threads; thread++) {
    unsigned int chunk = 0;
    unsigned chunk_limit; // It is possible that the last thread process less chunks
    
    if((thread == num_threads - 1) && (st->num_chunks%chunks_per_thread != 0))
      chunk_limit = st->num_chunks%chunks_per_thread;
    else
      chunk_limit = chunks_per_thread;

    depth_t min = 0, max = 0, partial_excess = 0;

    // Each thread traverses their chunks
    for(chunk = 0; chunk < chunk_limit; chunk++) {
      int16_t num_mins = 1; // Number of occurrences of the minimum value in the chunk
      unsigned int llimit = 0, ulimit = 0;
      unsigned int global_chunk = thread*chunks_per_thread+chunk;
      
      // Compute the limits of the current chunk
      if(st->num_chunks-1 < global_chunk) {
	llimit = 0;
	ulimit = 0;
      }
      else if(global_chunk == st->num_chunks-1 && chunk == (chunk_limit-1) && n % (st->num_chunks * st->s) != 0){
	llimit = thread*chunks_per_thread*st->s+(st->s*chunk);
	ulimit = n;
      }
      else {
	llimit = thread*chunks_per_thread*st->s + (st->s*chunk);
	ulimit = llimit + st->s;
	if(st->n < st->s)
	  ulimit = n;
	if(ulimit > st->n)
	  ulimit = st->n;
      }
      
      unsigned int symbol=0;
	
      for(symbol=llimit; symbol<ulimit; symbol++) {

	// Excess computation
	if(bit_array_get_bit(bit_array, symbol) == 0)
	  --partial_excess;
	else
	  ++partial_excess;

	// Minimum computation
	if(symbol==llimit) {
	  min = partial_excess; // By default the minimum value is the first excess value
	  max = partial_excess; // By default the maximum value is the first excess value
	  num_mins = 1;
	}
	else {
	  if(partial_excess < min) {
	    min = partial_excess;
	    num_mins = 1;
	  } else if(partial_excess == min)
	    num_mins++;

	  if(partial_excess > max)
	    max = partial_excess;	  
	}
      }

      if(global_chunk < st->num_chunks) {
	st->e_prime[thread*chunks_per_thread+chunk] = partial_excess;
	st->m_prime[st->internal_nodes + thread*chunks_per_thread+chunk] = min;
	st->M_prime[st->internal_nodes + thread*chunks_per_thread+chunk] = max;
	st->n_prime[st->internal_nodes + thread*chunks_per_thread+chunk] = num_mins;
      }
    }
  }

  /*
   * STEP 2.2: Computation of the final prefix computations (desired values)
   */
  for(unsigned int thread=1; thread < num_threads-1; thread++) {
    unsigned int global_chunk = thread*chunks_per_thread+chunks_per_thread-1;
    
    if(global_chunk < st->num_chunks) {  
      st->e_prime[global_chunk] +=
	st->e_prime[(thread-1)*chunks_per_thread+chunks_per_thread-1];
    }
  }  

  cilk_for(unsigned int thread=1; thread < num_threads; thread++) {
    unsigned int chunk = 0;
    unsigned int ul = chunks_per_thread;
    unsigned int global_chunk = thread*chunks_per_thread+chunks_per_thread-1;

    if(thread == num_threads-1) {
      ul = st->num_chunks - (num_threads-1)*chunks_per_thread;
      global_chunk = st->num_chunks-1;
    }

   /*
     * Note 1: Thread 0 does not need to update their excess values
     * Note 2:Thread 0 does not need to update the minimum value of its first chunk
     */
    if(global_chunk < st->num_chunks) {
      for(chunk=0; chunk < ul; chunk++) {
	if((thread == num_threads-1) || (chunk < chunks_per_thread -1))
	  st->e_prime[thread*chunks_per_thread+chunk] +=
	    st->e_prime[(thread-1)*chunks_per_thread+chunks_per_thread-1]; 
	st->m_prime[st->internal_nodes + thread*chunks_per_thread+chunk] +=
	  st->e_prime[(thread-1)*chunks_per_thread+chunks_per_thread-1]; 
	st->M_prime[st->internal_nodes + thread*chunks_per_thread+chunk] +=
	  st->e_prime[(thread-1)*chunks_per_thread+chunks_per_thread-1]; 
      }
    }
  }
    
  /*
   * STEP 2.3: Completing the internal nodes of the min-max tree
   */
      
  int p_level = ceil(log(num_threads)/log(st->k)); /* p_level = logk(num_threads), level at which each thread has at least one 
						  subtree to process in parallel */
  unsigned int num_subtrees = pow(st->k,p_level); /* num_subtrees = k^p_level, number of subtrees of the min-max tree 
						 that will be computed in parallel at level p_level.
						 num_subtrees is O(num_threads) */
  
  //unsigned int subtree = 0;

  uint total_chunks = st->internal_nodes + st->num_chunks;
  cilk_for(unsigned int subtree = 0; subtree < num_subtrees; subtree++) {
      for(int lvl = st->height-1; lvl >= p_level; lvl--){ //The current level that is being constructed.
	//Note: The last level (leaves) is already constructed
	unsigned int num_curr_nodes = pow(st->k, lvl-p_level); //Number of nodes at curr_level level that belong to the subtree
      
      for(unsigned int node = 0; node < num_curr_nodes; node++) {
  	unsigned int pos = pow(st->k,lvl)-1 + node + subtree*num_curr_nodes;// Position in the final array of 'node'.
  									    //Note: It should be less than the offset
  	unsigned int lchild = pos*st->k+1, rchild = (pos+1)*st->k; //Range of children of 'node' in the final array
	
  	/* for(unsigned int child = lchild; (child <= rchild) && (child < st->num_chunks); child++) { */
  	for(unsigned int child = lchild; (child <= rchild) && (child <
  	total_chunks); child++) {	  
  	  if(child == lchild){// first time
  	    st->m_prime[pos] = st->m_prime[child];
  	    st->M_prime[pos] = st->M_prime[child];
  	    st->n_prime[pos] = st->n_prime[child];
  	  }
  	  else {
  	    if(st->m_prime[child] < st->m_prime[pos]) {
  	      st->m_prime[pos] = st->m_prime[child];
  	      st->n_prime[pos] = 1;
	    }
	    else if(st->m_prime[child] == st->m_prime[pos])
	      st->n_prime[pos]++;
	    
  	    if(st->M_prime[child] > st->M_prime[pos])
  	      st->M_prime[pos] = st->M_prime[child];
  	  }
  	}
      }
    }
  }
   
  for(int lvl=p_level-1; lvl >= 0 ; lvl--){ // O(num_threads)
    
    unsigned int num_curr_nodes = pow(st->k, lvl); // Number of nodes at curr_level level that belong to the subtree
    unsigned int node = 0, child = 0;
    
    for(node = 0; node < num_curr_nodes; node++) {
      unsigned int pos = (pow(st->k,lvl)-1)/(st->k-1) + node; // Position in the final array of 'node'
      unsigned int lchild = pos*st->k+1, rchild = (pos+1)*st->k; // Range of children of 'node' in the final array
      for(child = lchild; child <= rchild; child++){
	if(st->m_prime[child] == st->M_prime[child])
	  continue;
	
	if(child == lchild) { // first time
	  st->m_prime[pos] = st->m_prime[child];
	  st->M_prime[pos] = st->M_prime[child];
	  st->n_prime[pos] = st->n_prime[child];
	}
	else {
	  if(st->m_prime[child] < st->m_prime[pos]) {
	    st->m_prime[pos] = st->m_prime[child];
  	    st->n_prime[pos] = 1;
	  }
	  else if(st->m_prime[child] == st->m_prime[pos])
	    st->n_prime[pos]++;

	  if(st->M_prime[child] > st->M_prime[pos])
	    st->M_prime[pos] = st->M_prime[child];
	}
      }
    }
  }
  
  /*
   * STEP 3: Computation of all universal tables
   */

  T = create_lookup_tables();

  return st;
}

int32_t sum(rmMt* st, int32_t idx){

  if(idx >= st->n)
    return -1;
  
  int32_t chk = idx/st->s;
  int32_t excess = 0;

  // Previous chunk
  if(chk)
    excess += st->e_prime[chk-1];
  
  int llimit = chk*st->s;
  int rlimit = (idx/8)*8;

  word_t j=0;
  for(j=llimit; j<rlimit; j+=8) {
#ifdef ARCH64
    int32_t sum_idx = (((st->bit_array)->words[j>>logW]) & (0xFFL<<(j&(word_size-1)))) >> (j&(word_size-1));
#else
    int32_t sum_idx = (((st->bit_array)->words[j>>logW]) & (0xFF<<(j&(word_size-1)))) >> (j&(word_size-1));
#endif

    excess += T->word_sum[sum_idx];
  }

  for(uint i=j; i<=idx; i++)
    excess += 2*bit_array_get_bit(st->bit_array,i)-1;

  return excess;
}

int32_t check_leaf(rmMt* st, int32_t i, int32_t d) {
  int end = (i/st->s+1)*st->s;
  int llimit = (((i)+8)/8)*8;
  int rlimit = (end/8)*8;
  int32_t excess = d;
  int32_t output;
  int32_t j = 0;

  for(j=i+1; j< min(end, llimit); j++){
    excess += 2*bit_array_get_bit(st->bit_array,j)-1;
    if(excess == d-1)
      return j;
  }

  for(j=llimit; j<rlimit; j+=8) {
    int32_t desired = d - 1 - excess; // desired value must belongs to the range [-8,8]
    
#ifdef ARCH64
    int32_t sum_idx = (((st->bit_array)->words[j>>logW]) & (0xFFL<<(j&(word_size-1)))) >> (j&(word_size-1));
#else
    int32_t sum_idx = (((st->bit_array)->words[j>>logW]) & (0xFF<<(j&(word_size-1)))) >> (j&(word_size-1));
#endif
    
    if (desired >= -8 && desired <= 8) {
    uint16_t ii = (desired+8<<8) + sum_idx;
        
    int8_t x = T->near_fwd_pos[ii];
    if(x < 8)
      return j+x;
  }
    excess += T->word_sum[sum_idx];
  }
  
  for (j=max(llimit,rlimit); j < end; ++j) {
    excess += 2*bit_array_get_bit(st->bit_array,j)-1;
    if (excess == d-1) {
      return j;
    }
  }
  
  return i-1;
}


int32_t check_sibling(rmMt* st, int32_t i, int32_t d) {
  int llimit = i;
  int rlimit = i+st->s;
  int32_t output;
  int32_t excess = st->e_prime[(i-1)/st->s];
  int32_t j = 0;

  for(j=llimit; j<rlimit; j+=8) {
    int32_t desired = d - 1 - excess; // desired value must belongs to the range [-8,8]  
    
#ifdef ARCH64
    int32_t sum_idx = (((st->bit_array)->words[j>>logW]) & (0xFFL<<(j&(word_size-1)))) >> (j&(word_size-1));
#else
    int32_t sum_idx = (((st->bit_array)->words[j>>logW]) & (0xFF<<(j&(word_size-1)))) >> (j&(word_size-1));
#endif
    
    if (desired >= -8 && desired <= 8) {
    uint16_t ii = (desired+8<<8) + sum_idx;
    
    int8_t x = T->near_fwd_pos[ii];
    if(x < 8)
      return j+x;
  }
    excess += T->word_sum[sum_idx];
  }
    
  return i-1;
}

int32_t fwd_search2(rmMt* st, int32_t i) {
    // Excess value up to the ith position 
    int32_t d = sum(st, i);
    
    int chunk = i / st->s;
    int32_t output;
    long j;
    
    // Case 1: Check if the chunk of i contains fwd_search(B, i, d)
    
    output = check_leaf(st, i, d);
    if(output > i)
      return output;
    
    // Case 2: The answer is not in the chunk of i, but it is in its sibling
    // (assuming a binary tree, if i%2==0, then its right sibling is at position
    // i+1)
    
    if(chunk%2 == 0) { // The current chunk has a right sibling
      // The answer is in the right sibling of the current node
      if(st->m_prime[chunk+1] <= d && d <= st->M_prime[chunk+1]) { 
	output = check_sibling(st, st->s*(chunk+1), d);
	if(output >= st->s*(chunk+1))
	  return output;
      }
    }
  
    // Case 3: It is necessary up and then down in the min-max tree
    
    long node = parent(st->internal_nodes + chunk); // Initial node
    // Go up the tree
    while (!is_root(node)) {
      if (is_left_child(node)) { // if the node is a left child
	node = right_sibling(node); // choose right sibling
	
	if (st->m_prime[node] <= d-1 && d-1 <= st->M_prime[node])
	  break;
      }
      node = parent(node); // choose parent
    }

    // Go down the tree
    if (!is_root(node)) { // found solution for the query
      while (!is_leaf(node, st->internal_nodes)) {
	node = left_child(node); // choose left child
	if (!(st->m_prime[node] <= d-1 && d-1 <= st->M_prime[node])) {
	  node = right_sibling(node); // choose right child == right sibling of the left child
	  if(st->m_prime[node] > d-1 || d-1 > st->M_prime[node]) {
	    return -1;
	  }
	}
      }

      chunk = node - st->internal_nodes;
      return check_sibling(st, st->s*chunk, d);
    }
    return -1;
}

// Check a leaf from left to right
int32_t check_leaf_r(rmMt* st, int32_t i, int32_t d) {
  int end = (i/st->s+1)*st->s;
  int llimit = (((i)+8)/8)*8;
  int rlimit = (end/8)*8;
  int32_t excess = d;
  int32_t output;
  int32_t j = 0;
  
  for(j=i+1; j< min(end, llimit); j++){
    excess += 2*bit_array_get_bit(st->bit_array,j)-1;
    if(excess == d-1)
      return j;
  }

  for(j=llimit; j<rlimit; j+=8) {
    int32_t desired = d - 1 - excess; // desired value must belongs to the range [-8,8]
    
#ifdef ARCH64
    int32_t sum_idx = (((st->bit_array)->words[j>>logW]) & (0xFFL<<(j&(word_size-1)))) >> (j&(word_size-1));
#else
    int32_t sum_idx = (((st->bit_array)->words[j>>logW]) & (0xFF<<(j&(word_size-1)))) >> (j&(word_size-1));
#endif
    
    if (desired >= -8 && desired <= 8) {
    uint16_t ii = (desired+8<<8) + sum_idx;
        
    int8_t x = T->near_fwd_pos[ii];
    if(x < 8)
      return j+x;
  }
    excess += T->word_sum[sum_idx];
  }
  
  for (j=max(llimit,rlimit); j < end; ++j) {
    excess += 2*bit_array_get_bit(st->bit_array,j)-1;
    if (excess == d-1) {
      return j;
    }
  }
  
  return i-1;
}

// Check siblings from left to right
int32_t check_sibling_r(rmMt* st, int32_t i, int32_t d) {
  int llimit = i;
  int rlimit = i+st->s;
  int32_t output;
  int32_t excess = st->e_prime[(i-1)/st->s];
  int32_t j = 0;

  for(j=llimit; j<rlimit; j+=8) {
    int32_t desired = d - excess; // desired value must belongs to the range [-8,8]  
    
#ifdef ARCH64
    int32_t sum_idx = (((st->bit_array)->words[j>>logW]) & (0xFFL<<(j&(word_size-1)))) >> (j&(word_size-1));
#else
    int32_t sum_idx = (((st->bit_array)->words[j>>logW]) & (0xFF<<(j&(word_size-1)))) >> (j&(word_size-1));
#endif
    
    if (desired >= -8 && desired <= 8) {
      uint16_t ii = (desired+8<<8) + sum_idx;
      
      int8_t x = T->near_fwd_pos[ii];
      
      if(x < 8)
	return j+x;
    }
    excess += T->word_sum[sum_idx];
  }
    
  return i-1;
}

int32_t fwd_search(rmMt* st, int32_t i, int32_t d) {
    // Excess value up to the ith position 
    int32_t target = sum(st, i) + d - 1;
    
    int chunk = i / st->s;
    int32_t output;
    long j;
    
    // Case 1: Check if the chunk of i contains fwd_search(bit_array, i, target)
    output = check_leaf_r(st, i, target);
    if(output > i)
      return output;
    
    // Case 2: The answer is not in the chunk of i, but it is in its sibling
    // (assuming a binary tree, if i%2==0, then its right sibling is at position i+1)
    if(chunk%2 == 0) { // The current chunk has a right sibling
      // The answer is in the right sibling of the current node
      if(st->m_prime[st->internal_nodes + chunk+1] <= target && target <=
	 st->M_prime[st->internal_nodes+ chunk+1]) {

	output = check_sibling_r(st, st->s*(chunk+1), target);
	if(output >= st->s*(chunk+1))
	  return output;
      }
    }
  
    // Case 3: It is necessary up and then down in the min-max tree
    long node = parent(chunk + st->internal_nodes); // Initial node
    // Go up the tree
    while (!is_root(node)) {
      if (is_left_child(node)) { // if the node is a left child
	node = right_sibling(node); // choose right sibling
	
	if (st->m_prime[node] <= target && target <= st->M_prime[node])
	  break;
      }
      node = parent(node); // choose parent
    }

    // Go down the tree
    if (!is_root(node)) { // found solution for the query
      while (!is_leaf(node, st->internal_nodes)) {
	node = left_child(node); // choose left child
	if (!(st->m_prime[node] <= target && target <= st->M_prime[node])) {
	  node = right_sibling(node); // choose right child == right sibling of the left child
	  if(st->m_prime[node] > target || target > st->M_prime[node]) {
	    return i;
	  }
	}
      }
      
      chunk = node - st->internal_nodes;

      return check_sibling_r(st, st->s*chunk, target);
    }
    return i;
}

int32_t find_close(rmMt* st, int32_t i){
  if(bit_array_get_bit(st->bit_array,i) == 0)
    return i;

  return fwd_search(st, i, 0);
}


// Naive implementation of fwd_search
int32_t naive_fwd_search(rmMt* st, int32_t i, int32_t d) {
  int begin = i+1;
  int end = st->n;
  int32_t excess = sum(st, i);
  int32_t target = excess + d - 1;
  int32_t j = 0;

  for(j=begin; j < end; j++) {
    excess += 2*bit_array_get_bit(st->bit_array,j)-1;
    if(excess == target)
      return j;
  }

  return i;
}

// Semi naive implementation of fwd_search
int32_t semi_fwd_search(rmMt* st, int32_t i, int32_t d) {
  int32_t excess = sum(st, i);
  int32_t target = excess + d - 1;
  int32_t j = 0;
  int chunk = i/st->s;

  int begin = i+1;
  int end = (chunk+1)*st->s;
  for(j=begin; j < end; j++) {
    excess += 2*bit_array_get_bit(st->bit_array,j)-1;
    if(excess == target)
      return j;
  }

  begin = chunk+1;
  end = st->num_chunks;
  for(j=begin; j<end;j++) {
    uint idx = st->internal_nodes + j;
    if(st->m_prime[idx] <= target && st->M_prime[idx] >= target) {
      chunk = j;
      break;
    }
  }

  begin = chunk*st->s;
  end = (chunk+1)*st->s;
  excess = st->e_prime[chunk-1];

  for(j=begin; j < end; j++) {
    excess += 2*bit_array_get_bit(st->bit_array,j)-1;
    if(excess == target)
      return j;
  }
  return i;
}

int32_t find_close_naive(rmMt* st, int32_t i){
  if(bit_array_get_bit(st->bit_array,i) == 0)
    return i;

  return naive_fwd_search(st, i, 0);
}

int32_t find_close_semi(rmMt* st, int32_t i){
  if(bit_array_get_bit(st->bit_array,i) == 0)
    return i;

  return semi_fwd_search(st, i, 0);
}

int32_t rank_0(rmMt* st, int32_t i) {
  // Excess value up to the ith position
  if(i >= st->n)
    i = st->n-1;
    int32_t d = sum(st, i);
  
  return (i+1-d)/2;
}


// Naive implementation of bwd_search
int32_t naive_bwd_search(rmMt* st, int32_t i, int32_t d) {
  int begin = 0;
  int32_t excess = sum(st, i);
  int32_t target = excess + d;
  int32_t j = 0;

  for(j=i; j >= begin; j--) {
    excess += 2*bit_array_get_bit(st->bit_array,j)-1;
    if(excess == target) {
      return j;
    }
  }

  return i;
}

// Semi implementation of bwd_search
int32_t semi_bwd_search(rmMt* st, int32_t i, int32_t d) {
  int32_t excess = sum(st, i);
  int32_t target = excess + d;
  int32_t j = 0;

  if(target == 0 && i == st->n-1)
    return 0;
  
  int chunk = i/st->s;
  int begin = i;
  int end = chunk*st->s;

  for(j=begin; j >= end; j--) {
    excess += 1 - 2*bit_array_get_bit(st->bit_array,j);
    if(excess == target)
      return j;
  }

  begin = chunk-1;
  end = 0;

  for(j=begin; j >= end; j--) {
    int idx = st->internal_nodes + j;

    if(st->m_prime[idx] <= target && st->M_prime[idx] >= target) {
      chunk = j;
      break;
    }
  }

  begin = (chunk+1)*st->s-1;
  end = chunk*st->s;
  excess = st->e_prime[chunk];

  // Special case 2
  if(excess == target)
    return begin+1;

  for(j=begin; j >= end; j--) {
    excess += 1 - 2*bit_array_get_bit(st->bit_array,j);
    if(excess == target)
      return j;
  }

  return i;
}

// Check a leaf from right to left
int32_t check_leaf_l(rmMt* st, int32_t i, int32_t target, int32_t excess) {
  int rlimit = (i/8)*8;
  int begin = (i/st->s)*st->s;
  int llimit = ((begin+8)/8)*8;
  if(llimit > rlimit)
    llimit = rlimit;
  int32_t output;
  int32_t j = 0;

  for(j=i; j >= max(rlimit, llimit); j--){
    excess += 2*bit_array_get_bit(st->bit_array,j)-1;
    if(excess == target) {
      return j;
    }
  }
  for(j = rlimit-8; j >= llimit; j-=8) {
    int32_t desired = excess - target; // desired value must belongs to the range [-8,8]
    
#ifdef ARCH64
    int32_t sum_idx = (((st->bit_array)->words[j>>logW]) & (0xFFL<<(j&(word_size-1)))) >> (j&(word_size-1));
#else
    int32_t sum_idx = (((st->bit_array)->words[j>>logW]) & (0xFF<<(j&(word_size-1)))) >> (j&(word_size-1));
#endif
    if (desired >= -8 && desired <= 8) {
      uint16_t ii = (desired+8<<8) + sum_idx;
      
      int8_t x = T->near_bwd_pos[ii];
      if(x < 8)
	return j+x;
    }
    excess += T->word_sum[sum_idx];
  }

  for (j=min(llimit,rlimit)-1; j >= begin; j--) {
    excess += 2*bit_array_get_bit(st->bit_array,j)-1;
    if (excess == target) {
      return j;
    }
  }
  
  return i;
}

// Check a left sibling
int32_t check_sibling_l(rmMt* st, int32_t i, int32_t excess, int32_t d) {
  int llimit = i;
  int rlimit = i+st->s;

  int32_t e = st->e_prime[i/st->s];
  int32_t output;
  int32_t j = 0;

  for(j = rlimit-8; j >= llimit; j-=8) {
    int32_t desired =  excess - d - e; // desired value must belongs to the range [-8,8]
    
#ifdef ARCH64
    int32_t sum_idx = (((st->bit_array)->words[j>>logW]) & (0xFFL<<(j&(word_size-1)))) >> (j&(word_size-1));
#else
    int32_t sum_idx = (((st->bit_array)->words[j>>logW]) & (0xFF<<(j&(word_size-1)))) >> (j&(word_size-1));
#endif
    if (desired >= -8 && desired <= 8) {
      uint16_t ii = (desired+8<<8) + sum_idx;
      
      int8_t x = T->near_bwd_pos[ii];
      if(x < 8)
	return j+x;
    }
    e -= T->word_sum[sum_idx];
  }
  
  return i-1;
}

int32_t bwd_search(rmMt* st, int32_t i, int32_t d) {
  int32_t excess = sum(st, i);
  int32_t target = excess + d;

  int chunk = i / st->s;
  int32_t output = i;
  long j;

  // Case 1: Check if the chunk of i contains bwd_search(bit_array, i, target)
  output = check_leaf_l(st, i, target, excess);
  if(output < i)
    return output;
  
  // Case 2: The answer is not in the chunk of i, but it is in its sibling
  // (assuming a binary tree, if i%2==1, then its left sibling is at position i-1)
  if(chunk%2 == 1) { // The current chunk has a left sibling
    // The answer is in the left sibling of the current node
    if(st->m_prime[st->internal_nodes + chunk - 1] <= excess-d+1 && excess-d+1 <=
       st->M_prime[st->internal_nodes + chunk - 1]) {
      
      output = check_sibling_l(st, st->s*(chunk-1), excess, d);
      if(output >= st->s*(chunk-1))
  	return output;
    }
  }

  // Case 3: It is necessary up and then down in the min-max tree
  long node = parent(chunk + st->internal_nodes); // Initial node
  // Go up the tree
  while (!is_root(node)) {
    if (is_right_child(node)) { // if the node is a left child
      node = left_sibling(node); // choose right sibling
      
      //      if (st->m_prime[node] <= target && target <= st->M_prime[node])
      if (st->m_prime[node] <= excess-d && excess-d <= st->M_prime[node])
	break;      
    }
    node = parent(node); // choose parent
  }

  // Go down the tree
  if (!is_root(node)) { // found solution for the query
    while (!is_leaf(node, st->internal_nodes)) {
      node = right_child(node); // choose right child

      if (!(st->m_prime[node] <= excess-d && excess-d <= st->M_prime[node])) {
	node = left_sibling(node); // choose left child == left sibling of the	right child
	
	if(st->m_prime[node] > excess-d || excess-d > st->M_prime[node]) {
	  return i;
	}
      }
    }

    chunk = node - st->internal_nodes;

    // Special case: if the result is at the beginning of chunk i, then,
    // the previous condition will select the chunk i-1
    //    if(st->e_prime[chunk] == target) { // If the last value (e') of chunk is equal
    if(st->e_prime[chunk] == excess-d) { // If the last value (e') of chunk is equal
    // to the target, then the answer is in the first position of the next chunk

      return (chunk+1)*st->s;
    }
    
    return check_sibling_l(st, st->s*chunk, excess, d);
  }
  else {// Special case: Pair of parentheses wrapping the parentheses sequence
        //(at positions 0 and n-1)
    if(i == st->n-1 && excess==0)
      output = 0;
  }

  return output;
}

int32_t find_open_naive(rmMt* st, int32_t i){
  if(bit_array_get_bit(st->bit_array,i) == 1)
    return i;

  return naive_bwd_search(st, i, 0);  
}

int32_t find_open(rmMt* st, int32_t i){
  if(bit_array_get_bit(st->bit_array,i) == 1)
    return i;

  return bwd_search(st, i, 0);  
}

int32_t find_open_semi(rmMt* st, int32_t i){
  if(bit_array_get_bit(st->bit_array,i) == 1)
    return i;

  return semi_bwd_search(st, i, 0);  
}

int32_t rank_1(rmMt* st, int32_t i) {
    // Excess value up to the ith position 
  if(i >= st->n)
    i = st->n-1;
  int32_t d = sum(st, i);
  
  return (i+1+d)/2;
}


int32_t check_chunk(rmMt* st, int32_t i, int32_t d) {
  int llimit = i;
  int rlimit = i+st->s;
  int32_t output;
  int32_t excess = st->e_prime[(i-1)/st->s];
  int32_t j = 0;

  for(j=llimit; j<rlimit; j+=8) {
    int32_t desired = d - 1 - excess; // desired value must belongs to the range [-8,8]  
    
#ifdef ARCH64
    int32_t sum_idx = (((st->bit_array)->words[j>>logW]) & (0xFFL<<(j&(word_size-1)))) >> (j&(word_size-1));
#else
    int32_t sum_idx = (((st->bit_array)->words[j>>logW]) & (0xFF<<(j&(word_size-1)))) >> (j&(word_size-1));
#endif
    
    if (desired >= -8 && desired <= 8) {
      uint16_t ii = (desired+8<<8) + sum_idx;
      
      int8_t x = T->near_fwd_pos[ii];
      if(x < 8)
	return j+x;
    }
    excess += T->word_sum[sum_idx];
  } 
    
  return i-1;
}

// ToDo: Implement it more efficiently
int32_t select_0(rmMt* st, int32_t i){

  int32_t j = 0;

  // The answer is after the position 2*i-1
  int32_t excess = sum(st,2*i-1);

  // Note: The answer is not beyond the position 2*i-1+depth_max, where
  // depth_max is the maximal depth (excess) of the input tree
  int32_t llimit = 2*i-1;
  int32_t rlimit = llimit + st->M_prime[0];
  int32_t d = 0;

  for (j=llimit+1; j <=rlimit; ++j,++d) {
    if (excess == d)
      return j-1;
    excess += 2*bit_array_get_bit(st->bit_array,j)-1;
  }

    return -1;
}

// ToDo: Implement it more efficiently
int32_t select_1(rmMt* st, int32_t i){
  int32_t j = 0;

  // Note: The answer is in the range [0,2*i-1] not beyond the position 2*i-1
  int32_t excess = 0;
  int32_t llimit = 0;
  int32_t rlimit = 2*i-1;
  int32_t d = 2*i-1;

  for (j=llimit; j <=rlimit; ++j,--d) {
    excess += 2*bit_array_get_bit(st->bit_array,j)-1;
    if (excess == d)
      return j;
  }

    return -1;
}


int32_t match(rmMt* st, int32_t i) {
  if(bit_array_get_bit(st->bit_array,i))
    return find_close(st, i);
  else
    return find_open(st, i);
}

int32_t match_naive(rmMt* st, int32_t i) {
  if(bit_array_get_bit(st->bit_array,i))
    return find_close_naive(st, i);
  else
    return find_open_naive(st, i);
}

int32_t match_semi(rmMt* st, int32_t i) {
  if(bit_array_get_bit(st->bit_array,i))
    return find_close_semi(st, i);
  else
    return find_open_semi(st, i);
}


int32_t parent_t(rmMt* st, int32_t i) {
  if(!bit_array_get_bit(st->bit_array,i))
    i = find_open(st, i);
  
  return bwd_search(st, i, 2);
}

int32_t depth(rmMt* st, int32_t i) {
  return 2*rank_1(st, i)-i-1;
}

int32_t first_child(rmMt* st, int32_t i) {
  if(i >= st->n-1)
    return -1;

  if(!bit_array_get_bit(st->bit_array,i))
    return -1;
  
  if(bit_array_get_bit(st->bit_array,i+1))
    return i+1;
  else 
    return -1;
}

int32_t next_sibling(rmMt* st, int32_t i) {
  if(i >= st->n-1)
    return -1;
  
  if(!bit_array_get_bit(st->bit_array,i))
    return -1;

  i = find_close(st, i);  
  
  if(bit_array_get_bit(st->bit_array,i+1))
    return i+1;
  else 
    return -1;
}

int32_t is_leaf_t(rmMt* st, int32_t i) {
  if(i >= st->n-1)
    return 0;

  if(!bit_array_get_bit(st->bit_array,i))
    return 0;

  if(!bit_array_get_bit(st->bit_array,i+1))
    return 1;
  
  return 0;
}

ulong size_rmMt(rmMt *st) {
  ulong sizeRmMt = sizeof(rmMt);
  ulong sizeBitArray = st->bit_array->num_of_bits/8;
  ulong sizePrimes = 3*((st->num_chunks + st->internal_nodes)*sizeof(int16_t)) +
    st->num_chunks*sizeof(int16_t);

  return sizeRmMt + sizeBitArray + sizePrimes;
}


void print_partial_bit_array(rmMt *st, int32_t i) {
  for(uint j=0; j<=i; j++)
    printf("%d", bit_array_get_bit(st->bit_array,j));
  printf("\n");
}

void print_partial_bit_array2(bitRankW32Int *A, int32_t i) {
  for(uint j=0; j<=i; j++)
    printf("%d", isBitSet(A,j));
  printf("\n");
}
