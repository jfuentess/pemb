/*
 * Adaptation of the file "bp_support_algorithm.hpp" of sdsl-lite
 * sdsl-lite: https://github.com/simongog/sdsl-lite
 */

#include "lookup_tables.h"
#include "util.h"
#include <stdio.h>
#include <stdlib.h>

// Given an excess value x in [-8,8] and a 8-bit
// word w interpreted as parentheses sequence.
// near_fwd_pos[(x+8)<<8 | w] contains the minimal position
// p in [0..7] where the excess value x is reached, or 8
// if x is not reached in w.

lookup_table* create_lookup_tables() {
  
  lookup_table* T = (lookup_table *)malloc(sizeof(lookup_table));
  
  //  int32_t x;
  cilk_for (int32_t x = -8; x < 8; ++x) {
    for (uint16_t w=0; w < 256; ++w) {
      uint16_t i = (x+8)<<8|w;
      T->near_fwd_pos[i] = 8;
      int8_t p=0;
      int8_t excess = 0;
      do {
  	excess += 1-2*((w&(1<<p))==0);
  	if (excess == x) {
  	  T->near_fwd_pos[i] = p;
  	  break;
  	}
  	++p;
      } while (p < 8);
      
      T->near_bwd_pos[i] = 8;
      p = 7;
      excess = 0;
      do {
      	excess += 1-2*((w&(1<<p))>0);
      	if (excess == x) {
      	  T->near_bwd_pos[i] = p;
      	  break;
      	}
      	--p;
      } while (p > -1);
    }
  }
  
  cilk_for (uint16_t w = 0; w < 256; ++w) {
    uint16_t p;
    int8_t excess = 0;
    uint32_t ones = 0;
    for (p=0; p<8; ++p) {
      ones += (w&(1<<p))!=0;
      excess += 1-2*((w&(1<<p))==0);
    }
    T->word_sum[w] = excess;
  }
    
  cilk_for(uint16_t w = 0; w < 256; ++w) {
    uint32_t packed_mins[8];
    uint32_t packed_maxs[8];

    int8_t excess = 0;
    int8_t rev_excess = 0;
    int32_t min_excess_of_open = 17;
    int32_t min_excess_of_open_pos = 0;
    uint32_t ones = 0;
    T->min[w] = 8;
    packed_mins[0] = 0x99999999U;
    packed_maxs[0] = 0x99999999U;
    uint16_t p;
    
    for (p=0; p<8; ++p) {
      ones += (w&(1<<p))!=0;
      excess += 1-2*((w&(1<<p))==0);
      if (excess <= T->min[w]) {
    	T->min[w] = excess;
    	T->min_pos_max[w] = p;
      }
      if (excess < 0 && packed_mins[-excess-1] == 9) {
    	packed_mins[-excess-1] = p;
      }
      if (w&(1<<p) && excess+8 <= min_excess_of_open) {
    	min_excess_of_open     = excess+8;
    	min_excess_of_open_pos = p;
      }
      rev_excess += 1-2*((w&(1<<(7-p)))>0);
      if (rev_excess < 0 && packed_maxs[-rev_excess-1] == 9) {
    	packed_maxs[-rev_excess-1] = 7-p;
      }
    }
    T->word_sum[w] = excess;
    T->min_match_pos_packed[w] = packed_mins[0];
    T->max_match_pos_packed[w] = packed_maxs[0];
    T->min_open_excess_info[w] = (min_excess_of_open) |
      (min_excess_of_open_pos << 8) |
      (ones << 12);
  }
  
  return T;
}
