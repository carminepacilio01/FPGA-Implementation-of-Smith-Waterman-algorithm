// Test Bench for the implementation of the Smith-Waterman (Gotoh)

#include <stdio.h>
#include <stdlib.h>

#include "smith_waterman.h"

////////////////////////////////////////////////////
void to_string(char *seqA, char *seqB, int ws, int wd, int gap_opening, int enlargement);
void printMatrix(short int *P, short int *D, short int *Q, int lenA, int lenB);
void fprintMatrix(short int *P, short int *D, short int *Q, int lenA, int lenB);

int main(int argc, char const *argv[]) {
	int lenA = 4, lenB = 4; 												//	length of the sequences

	int retval = 0;															//	used to compare results

	const char *seqA = "CGGA";
	const char *seqB = "CCGA"; 												//	sequences
	const short int wd = 1; 												//	match score
	const short int ws = -1;												//	mismatch score
	const short int gap_opening = -3;
	const short int enlargement = -1;

	short int *P = (short int *) malloc((lenA + 1) * (lenB + 1) * sizeof(short int));
	short int *D = (short int *) malloc((lenA + 1) * (lenB + 1) * sizeof(short int));
	short int *Q = (short int *) malloc((lenA + 1) * (lenB + 1) * sizeof(short int));

	//--------------------------------------------------------------------------------//
	//	Printing current configuration
	to_string(seqA, seqB, ws, wd, gap_opening, enlargement);

	//Calling kernel
	smithWaterman(lenA, seqA, lenB, seqB, P, D, Q, wd, ws, gap_opening, enlargement);

	//checking testBench
	printMatrix(P, D, Q, lenA, lenB);
	int itr, i, j;

	retval = system("diff --brief -w test_out.txt test_out_golden.txt");
	if (retval != 0) {
		printf("-- Test failed  !!! --\n");
		retval=1;
	} else {
		printf("-- Test passed ! --\n");
	}

	free(P);
	free(D);
	free(Q);

	return 0;
}


void to_string(char *seqA, char *seqB, int ws, int wd, int gap_opening, int enlargement) {
	printf("\n+++++++++++++++++++++\n");
	printf("+ Sequence A: %s\n", seqA);
	printf("+ Sequence B: %s\n", seqB);
	printf("+ Match Score: %d\n", wd);
	printf("+ Mismatch Score: %d\n", ws);
	printf("+ Gap Opening: %d\n", gap_opening);
	printf("+ Enlargement: %d\n", enlargement);
	printf("+++++++++++++++++++++\n");
}

void printMatrix(short int *P, short int *D, short int *Q, int lenA, int lenB) {
	int i;

	printf("--- P ---\n");
	for (i = 0; i < lenA + 1; ++i) {
		for (int j = 0; j < lenB + 1; ++j) {
			if(i == 0 && j == 0)
				printf("\t \t");
			else if (i == 0)
				printf("\t -00 \t");
			else if (j == 0)
				printf("\t - \t");
			else
				printf("\t %d \t", P[i * MAX_DIM + j]);
		}
		printf("\n");
	}

	printf("--- D ---\n");
	for (i = 0; i < lenA + 1; ++i) {
		for (int j = 0; j < lenB + 1; ++j) {
			printf("\t %d \t", D[i * MAX_DIM + j]);
		}
		printf("\n");
	}

	printf("--- Q ---\n");
	for (i = 0; i < lenA + 1; ++i) {
		for (int j = 0; j < lenB + 1; ++j) {
			if(i == 0 && j == 0)
				printf("\t \t");
			else if (j == 0)
				printf("\t -00 \t");
			else if (i == 0)
				printf("\t - \t");
			else
				printf("\t %d \t", Q[i * MAX_DIM + j]);
		}
		printf("\n");
	}
}

void fprintMatrix(short int *P, short int *D, short int *Q, int lenA, int lenB) {
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
				fprintf(test_out, "\t %d \t", P[i * MAX_DIM + j]);
		}
		fprintf(test_out, "\n");
	}

	fprintf(test_out, "--- D ---\n");
	for (i = 0; i < lenA + 1; ++i) {
		for (int j = 0; j < lenB + 1; ++j) {
			fprintf(test_out, "\t %d \t", D[i * MAX_DIM + j]);
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
				fprintf(test_out, "\t %d \t", Q[i * MAX_DIM + j]);
		}
		fprintf(test_out, "\n");
	}

	fclose(test_out);
}

