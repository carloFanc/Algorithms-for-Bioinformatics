#!/usr/bin/env python 
#
#Implementation of Smith Watermann algorithm modified for the exam
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

#procedure starts at the score inserted as parameter and proceeds until a cell with score zero is encountered
def traceback(scoring_matrix, r, c, seq1, seq2, match_score, mismatch_score, gap_pen):
 
	score = scoring_matrix[r,c] 
	score_record = scoring_matrix[r,c]
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

	return (seq1_aligned, seq2_aligned, score_record)

#procedure which finds all the possible alignments with the required parameters
def find_all_alignments(scoring_matrix): 
	all_allign1, all_allign2, scores, lengths, gaps, matches = [],[],[],[],[],[]
	for r in range(scoring_matrix.shape[0]):		 
		for c in range(scoring_matrix.shape[1]):
			countMatches = 0	 
			if scoring_matrix[r,c] > 5: #first verification in which I check if the score is greater than 5 
				s1,s2,sc = traceback(scoring_matrix, r, c, seq1, seq2, match_score, mismatch_score, gap_pen)
				
				for k in range(0,len(s1)):
					if s1[k] == s2[k]:
						countMatches += 1
				#second verification in which I check if the matches are more than 4 and the gaps more than 0 
				if countMatches > 4 and ((s1.count("-")+ s2.count("-")) > 0): 
					all_allign1.append(s1) 
					all_allign2.append(s2)
					scores.append(sc)
					lengths.append(len(s1))
					gaps.append(s1.count("-")+ s2.count("-"))
					matches.append(countMatches)
					
	return all_allign1, all_allign2, scores, lengths, gaps, matches
	
#sorting procedure by decreasing score
def sortbyscore(alignments1, alignments2, scores, lengths, gaps, matches):
	for i in range(len(scores)):  
		for j in range(i+1, len(scores)):   
			 if (scores[i] < scores[j]):   
				 
				 temp1= lengths[i]
				 lengths[i]= lengths[j]
				 lengths[j]=temp1
				 
				 temp2=alignments1[i]  
				 alignments1[i]=alignments1[j]  
				 alignments1[j]=temp2
				
				 temp3=alignments2[i]  
				 alignments2[i]=alignments2[j]  
				 alignments2[j]=temp3 
				
				 temp4=scores[i]  
				 scores[i]=scores[j]  
				 scores[j]=temp4 
				
				 temp5=gaps[i]  
				 gaps[i]=gaps[j]  
				 gaps[j]=temp5  
				 
				 temp6=matches[i]  
				 matches[i]=matches[j]  
				 matches[j]=temp6  

#create the pipes which highlight the matches				 
def get_pipe(s1, s2):  
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
	
	# matrix filling with scores
	scoring_matrix = scoring(scoring_matrix, seq1, seq2, match_score, mismatch_score, gap_pen)
	print("\nScoring matrix: ")
	print(scoring_matrix)
	
	#finding all possible alignments with the required parameters (matches >4, gap>0, score>5)
	alignments1, alignments2, scores, lengths, gaps, matches = find_all_alignments(scoring_matrix)
	
	#sort by decreasing score
	sortbyscore(alignments1, alignments2, scores, lengths, gaps, matches) 
	
	#print all the results
	print("\nList of all", len(alignments1), "possible alignments sorted by decreasing score: ")
	for i in range(len(alignments1)):
		print("\nScore:", scores[i],"Length:", lengths[i], "Gaps:", gaps[i],"Matches:", matches[i] )
		pipe = get_pipe(alignments1[i], alignments2[i])
		print(alignments1[i])
		print(pipe)
		print(alignments2[i]) 
		
	

 
