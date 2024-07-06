// Test Bench fot the implementation of the Smith-Waterman (Gotoh)

#include <stdio.h>
#include <stdlib.h>
#include "smith_waterman.h"

////////////////////////////////////////////////////
void to_string(char * seqA, char * seqB, int ws, int wd, int gap_opening, int enlargement);
void printMatrix(int **P, int **D, int **Q, int lenA, int lenB);

int main(int argc, char const *argv[]) {
	const short int lenA = 4;
	const short int lenB = 4; 				//length of the sequences

	int retval = 0;							//used to compare results
	
	char *seqA = "CGGA"; 					
	char *seqB = "CCGA"; 					//sequences

	const short int wd = 1; 				//match score
	const short int ws = -1;				//mismatch score
	const short int gap_opening = -3;
	const short int enlargement = -1;

	int **P, **D, **Q; 						//Alignement Matrix
	//-----------------------------------------//

	//Setting up the matrix
	P = (int **)malloc((lenA + 1) * sizeof(int*));
	D = (int **)malloc((lenA + 1) * sizeof(int*));
	Q = (int **)malloc((lenA + 1) * sizeof(int*));

	for(int i = 0; i < lenA + 1; i++) {
		P[i] = (int *)malloc((lenB + 1) * sizeof(int));
		D[i] = (int *)malloc((lenB + 1) * sizeof(int));
		Q[i] = (int *)malloc((lenB + 1) * sizeof(int));
	}

	smithWaterman(P, D, Q, wd, ws, gap_opening, enlargement, lenA, seqA, lenB, seqB);
	to_string(seqA, seqB, ws, wd, gap_opening, enlargement);

	//cheging testBench esit
	printMatrix(P, D, Q, lenA, lenB);
	retval = system("diff --brief -w test_out.txt test_out_golden.txt");
	if (retval != 0) {
		printf("\x1b[31m Test failed  !!! \x1b[0m \n"); 
		retval=1;
	} else {
		printf("\x1b[32m Test passed ! \x1b[0m \n");
	}

	return 0;
}


void to_string(char * seqA, char * seqB, int ws, int wd, int gap_opening, int enlargement) {
	printf("\n+++++++++++++++++++++\n");
	printf("+ Sequence A: %s\n", seqA);
	printf("+ Sequence B: %s\n", seqB);
	printf("+ Match Score: %d\n", wd);
	printf("+ Mismatch Score: %d\n", ws);
	printf("+ Gap Opening: %d\n", gap_opening);
	printf("+ Enlargement: %d\n", enlargement);
	printf("+++++++++++++++++++++\n");
}

void printMatrix(int **P, int **D, int **Q, int lenA, int lenB){
	int i;

	FILE *test_out = fopen("test_out.txt", "w");
	if (test_out == NULL) {
		perror("Error while opening the file");
		return;
	}

	fprintf(test_out, "--- P ---\n");
	for (i = 0; i < lenA + 1; ++i) {
		for (int j = 0; j < lenB + 1; ++j) {
			if(i == 0 && j == 0)
				fprintf(test_out, "\t \t");
			else if (i == 0)
				fprintf(test_out, "\t -00 \t");
			else if (j == 0)
				fprintf(test_out, "\t - \t");
			else
				fprintf(test_out, "\t %d \t", P[i][j]);
		}
		fprintf(test_out, "\n");
	}

	fprintf(test_out, "--- D ---\n");
	for (i = 0; i < lenA + 1; ++i) {
		for (int j = 0; j < lenB + 1; ++j) {
			fprintf(test_out, "\t %d \t", D[i][j]);
		}
		fprintf(test_out, "\n");
	}

	fprintf(test_out, "--- Q ---\n");
	for (i = 0; i < lenA + 1; ++i) {
		for (int j = 0; j < lenB + 1; ++j) {
			if(i == 0 && j == 0)
				fprintf(test_out, "\t \t");
			else if (j == 0)
				fprintf(test_out, "\t -00 \t");
			else if (i == 0)
				fprintf(test_out, "\t - \t");
			else
				fprintf(test_out, "\t %d \t", Q[i][j]);
		}
		fprintf(test_out, "\n");
	}

	fclose(test_out);
}