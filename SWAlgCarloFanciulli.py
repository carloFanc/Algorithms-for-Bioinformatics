#!/usr/bin/env python
#
#Version of Smith Watermann algorithm modified for the Exam of 24/06
#@author: Carlo Fanciulli 198793

import numpy as np 
import argparse
import sys
    
#procedure filling each matrix position with calculated scores 
def scoring(scoring_matrix, seq1, seq2, match_score, mismatch_score, gap_pen):
    # i = row , j = column 
	for i in range(1, len(seq1) + 1):		
		for j in range(1, len(seq2) + 1):	
			result = 0;   
			if(seq1[i-1] == seq2[j-1]):  
				result = match_score + scoring_matrix[i-1, j-1]	 
			else:    
				result = mismatch_score + scoring_matrix[i-1, j-1]	 
		
			gap_r = gap_pen + scoring_matrix[i-1, j] 
			gap_c = gap_pen + scoring_matrix[i, j-1]  

			# assign the max value between result, gaps and zero.
			scoring_matrix[i, j] = max(result, gap_r, gap_c, 0)

	return scoring_matrix

#procedure starts at the score inserted as parameter(the highest) and proceeds until a cell with score zero is encountered
def traceback(scoring_matrix, seq1, seq2, match_score, mismatch_score, gap_pen, score):

	r,c = np.where(scoring_matrix == score) 
    
   # it could be more than one max  
	if (r[0]):
		r = r[0]
	if (c[0]):
		c = c[0]

	score = scoring_matrix[r,c] 
	seq1_aligned = ""
	seq2_aligned = ""
  
	while(score > 0): # because in the local alignment the last value must be 0 

		if ((seq1[r-1] == seq2[c-1] and ((scoring_matrix[r-1,c-1] + match_score) == scoring_matrix[r,c])) or # match
			(seq1[r-1] != seq2[c-1] and ((scoring_matrix[r-1,c-1] + mismatch_score) == scoring_matrix[r,c]))): # mismatch
			seq1_aligned = seq1[r-1] + seq1_aligned
			seq2_aligned = seq2[c-1] + seq2_aligned
			r -= 1
			c -= 1 

		elif (scoring_matrix[r-1,c] + gap_pen == scoring_matrix[r,c]): # gap from the above
			seq1_aligned = seq1[r-1] + seq1_aligned
			seq2_aligned = '-' + seq2_aligned
			r -= 1 
	
		elif (scoring_matrix[r,c-1] + gap_pen == scoring_matrix[r,c]): # gap from the left
			seq1_aligned = '-' + seq1_aligned
			seq2_aligned = seq2[c-1] + seq2_aligned
			c -= 1 
		  
		score = scoring_matrix[r,c] 

	return (seq1_aligned, seq2_aligned)

def get_pipe(s1,s2):
    pipe = ""
    for k in range(0,len(s1)):
        mid = " "
        if s1[k] == s2[k]:
            mid = "|"
        pipe = pipe + mid
    return pipe 
 
if __name__ == '__main__':	 
	
	#inputs 
	parser = argparse.ArgumentParser()
	parser.add_argument("sequence1", help="first sequence")
	parser.add_argument("sequence2", help="second sequence")

	args = parser.parse_args()
	 
	seq1 = args.sequence1
	seq2 = args.sequence2
	#seq1 = 'TGTTACGG'	 
	#seq2 = 'GGTTGAGTA'	 
	print("Sequence_1 = \"" + seq1 + "\"")
	print("Sequence_2 = \"" + seq2 + "\"")

	#match and mismatch scores and gap penalty
	match_score = 3
	mismatch_score = -3
	gap_pen = -2
 
    #inizialization of the scoring matrix 
	scoring_matrix = np.zeros((len(seq1) + 1 , len(seq2) + 1), np.int)
	
	#matrix filling with scores
	scoring_matrix = scoring(scoring_matrix, seq1, seq2, match_score, mismatch_score, gap_pen)
	print("\nScoring matrix: ")
	print(scoring_matrix)
	 
	#find highest score
	score = np.amax(scoring_matrix)
	print("\nScore: " + str(score))
	 
	#traceback
	s1,s2 = traceback(scoring_matrix, seq1, seq2, match_score, mismatch_score, gap_pen, score)
	
	#create the pipes which highlight the alignment 
	pipe = get_pipe(s1,s2)
	
	#result
	print("\nResult: ")
	print(s1)
	print(pipe)
	print(s2)
 
