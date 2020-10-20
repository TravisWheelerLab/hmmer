#include <stdio.h>
#include "hmmer.h"

#define AA_COUNT 21
#define codon_DNA_length 192

// int sm_get_row(int index, int codon_count_array[]);
//void clear_array(long double array[], int l);

int build_aa_transition_matrix(float DNA_transitions[][4], long double aa_trans_table[AA_COUNT][AA_COUNT])
{
  typedef enum {A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y, ST} AA; 

  int aa_codon_count[AA_COUNT] = {4, 2, 2, 2, 2, 4, 2, 3, 2, 6, 1, 2, 4, 2, 6, 6, 4, 4, 1, 2, 3}; 

  int codon_DNA_table[] = {
    2,1,0, 2,1,1, 2,1,2, 2,1,3,
    3,2,1, 3,2,3,
    2,0,1, 2,0,3,
    2,0,0, 2,0,2,
    3,3,1, 3,3,3,
    2,2,0, 2,2,1, 2,2,2, 2,2,3,
    1,0,1, 1,0,3, 
    0,3,0, 0,3,1, 0,3,3,
    0,0,0, 0,0,2,
    1,3,0, 1,3,1, 1,3,2, 1,3,3, 3,3,0, 3,3,2,
    0,3,2,
    0,0,1, 0,0,3,
    1,1,0, 1,1,1, 1,1,2, 1,1,3,
    1,1,0, 1,0,2,
    0,2,0, 0,2,2, 1,2,0, 1,2,1, 1,2,2, 1,2,3,
    0,2,1, 0,2,3, 3,1,0, 3,1,1, 3,1,2, 3,1,3,
    0,1,0, 0,1,1, 0,1,2, 0,1,3, 
    2,3,0, 2,3,1, 2,3,2, 2,3,3,
    3,2,2,
    3,0,1, 3,0,3,
    3,0,0, 3,0,2, 3,2,0,
  };

  int i, j, k, l, m, row_1_index = 0, row_2_index, *codon_1, *codon_2; 

  int nucleotide_index_1, nucleotide_index_2;

  long double aa_trans_prob, codon_trans_prob, nuc_trans_prob;

  long double aa_tr_pr_row[AA_COUNT], aa_tr[AA_COUNT];

  for (i = 0; i < AA_COUNT; i++) {
    row_1_index = sm_get_row(i, aa_codon_count); 
    
    for (j = 0; j < AA_COUNT; j++) {
      row_2_index = sm_get_row(j, aa_codon_count);
      aa_trans_prob = 0;
      
      for(k = 0; k < aa_codon_count[i]; k++) {
        codon_1 = codon_DNA_table + row_1_index + 3*k;
        codon_trans_prob = 0;

        for(l = 0; l < aa_codon_count[j]; l++) {
          codon_2 = codon_DNA_table + row_2_index + 3*l;
          nuc_trans_prob = 1;
          
          for (m = 0; m < 3; m++) { 
            nuc_trans_prob *= DNA_transitions[*(codon_1+m)][*(codon_2+m)]; 
//printf("Nuc Freqs: %d -> %d = %f\n",  *(codon_1+m), *(codon_2+m), DNA_transitions[*(codon_1+m)][*(codon_2+m)]);
          }
          //printf("%Lf += %Lf\n\n", codon_trans_prob, nuc_trans_prob);
          codon_trans_prob += nuc_trans_prob;
        }
        
        //printf("\n!! %Lf += (%Lf / %d) !!\n\n",  aa_trans_prob, codon_trans_prob, aa_codon_count[j]);
        aa_trans_prob += (codon_trans_prob / aa_codon_count[i]) ;
      }
      
      //printf("!!! %d->%d TRANS PROB: %Lf !!!\n\n", i, j, aa_trans_prob);  
      aa_trans_table[i][j] = aa_trans_prob; 
    }
  }

  return 0;
}

int sm_get_row(int index, int codon_count_array[]) 
{
  if (index == 0)
    return 0;

  int sum = 0;

  for (index -= 1; index >= 0; index--) {
    sum += codon_count_array[index]*3;   
  }

  return sum; 
}

