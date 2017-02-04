/*
 * Adaptation of the file "bp_support_algorithm.hpp" of sdsl-lite
 * sdsl-lite: https://github.com/simongog/sdsl-lite
 */

#ifndef LOOKUP_TABLES_H
#define LOOKUP_TABLES_H

#include <stdint.h>

struct _lookup_table {
  // Given an excess value x in [-8,8] and a 8-bit
  // word w interpreted as parentheses sequence.
  // near_fwd_pos[(x+8)<<8 | w] contains the minimal position
  // p in [0..7] where the excess value x is reached, or 8
  // if x is not reached in w.
  uint8_t near_fwd_pos[(8-(-8))*256];
  
  // Given an excess value of x in [-8,8] and a 8-bit
  // word w interpreted as parentheses sequence.
  // near_bwd_pos[(x+8)<<8 | w] contains the maximal position
  // p in [0..7] where the excess value x is reached, or 8
  // if x is not reached in w.
  uint8_t near_bwd_pos[(8-(-8))*256];
  
  // Given a 8-bit word w. word_sum[w] contains the
  // excess value of w.
  int8_t word_sum[256];
  
  // Given a 8-bit word w. min[w] contains the
  // minimal excess value in w.
  int8_t min[256];
  
  // Given a 8-bit word w. min_pos_max[w] contains
  // the maximal position p in w, where min[w] is
  // reached
  int8_t min_pos_max[256];
  
  // Given an excess value x in [1,8] and a 8-bit
  // word w interpreted as parentheses sequence.
  // min_match_pos_packed[w]:[(x-1)*4,x*4] contains
  // the minimal position, where excess value
  // -x is reached and 9, if there is no such position.
  uint32_t min_match_pos_packed[256];
  
  // Given an excess value x in [1,8] and a 8-bit
  // word w interpreted as parentheses sequence.
  // max_match_pos_packed[w]:[(x-1)*4,x*4] contains
  // the maximal position, where excess value
  // -x is reached and 9, if there is no such position.
  uint32_t max_match_pos_packed[256];
  
  // Given a 8-bit word w. x=min_and_info[w] contains
  // the following information.
  // * [0..7] the minimum excess value in w + 8 of an opening parenthesis
  // * [8..11] the maximal position of the minimal excess value
  // * [12..15] the number of ones in the word
  // if w != 0, and 17 for w=0.
  uint16_t min_open_excess_info[256];
  
  
};

typedef struct _lookup_table lookup_table;

lookup_table * create_lookup_tables();

#endif // LOOKUP_TABLES_H
