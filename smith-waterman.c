// Implementation of the Smith-Waterman (Gotoh)      
//
//  Author: Carmine Pacilio							 
//	Year: 2024										  
//	Version: 1.0									  

///////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

//////////////////////////////////////////////////////
//matrix type
typedef int ** Matrix;

typedef struct node
{
	char baseA;
	char baseB;
	char match;
	struct node * next;
}t_node;

typedef t_node * alignment;
////////////////////////////////////////////////////

int gap_penalty(int k, int gap_opening, int enlargement);
void traceback(Matrix P, Matrix D, Matrix Q, int i, int j, char * seqA, char * seqB);

//For Debug
void to_string(char * seqA, char * seqB, int ws, int wd, int gap_opening, int enlargement);
void printMatrix(Matrix P, Matrix D, Matrix Q, int lenA, int lenB);

int main(int argc, char const *argv[]) {
	int lenA, lenB; //length of the sequences
	char * seqA; //sequences
	char * seqB; //sequences
	int wd; //match score
	int ws;	//mismatch score
	int gap_opening;
	int enlargement;

	Matrix P, D, Q; //Alignement Matrix

	//helper for traceback
	int max = INT_MIN;
	int occur = 0;

	// -- Input Phase -- //
	//printf("Insert the lenght of sequence A: ");
	scanf("%d", &lenA);
	seqA = malloc(lenA * sizeof(char) + 1);
	//printf("Insert sequence A: ");
	scanf("%s", seqA);

	//printf("Insert the lenght of sequence B: ");
	scanf("%d", &lenB);
	seqB = malloc(lenB * sizeof(char) + 1);
	//printf("Insert sequence B: ");
	scanf("%s", seqB);

	//printf("Insert the match score: ");
	scanf("%d", &wd);
	//printf("Insert the mismatch score: ");
	scanf("%d", &ws);
	//printf("Insert the gap opening: ");
	scanf("%d", &gap_opening);
	//printf("Insert the gap enlargement: ");
	scanf("%d", &enlargement);

	//Print the current settings
	to_string(seqA, seqB, wd, ws, gap_opening, enlargement);
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

	for (int i = 0; i < lenA + 1; ++i) {
		P[i][0] = INT_MIN / 2;
		Q[i][0] = INT_MIN / 2;
		D[i][0] = 0;

	}
	for (int i = 0; i < lenB + 1; ++i) {
		P[0][i] = INT_MIN / 2;
		Q[0][i] = INT_MIN / 2;
		D[0][i] = 0;
	}
	D[0][0] = 0;

	//Filling matrix phase
	for (int i = 1; i < lenA + 1; ++i) {
		for (int j = 1; j < lenB + 1; ++j) {
			//Calculate P = min{Di-1,j + g(1); Pi-1,j}
			P[i][j] = (D[i-1][j] + gap_penalty(1, gap_opening, enlargement) > P[i-1][j] + enlargement) ? 
			D[i-1][j] + gap_penalty(1, gap_opening, enlargement) : P[i-1][j] + enlargement;

			//Calculate Q = min{Di,j-1 + g(1); Qi,j-1}
			Q[i][j] = (D[i][j-1] + gap_penalty(1, gap_opening, enlargement) > Q[i][j-1] + enlargement) ?
			D[i][j-1] + gap_penalty(1, gap_opening, enlargement) : Q[i][j-1] + enlargement;

			//Calculate D = min{Di-1,j-1 + w(ai, bj); Pi,j; Qij}
			int tmp = D[i-1][j-1];
			tmp += (seqA[i-1] == seqB[j-1]) ? wd : ws;

			if (tmp > P[i][j] && tmp > Q[i][j] && tmp > 0){
				D[i][j] = tmp;
			} else if (P[i][j] > Q[i][j] && P[i][j] > 0) {
				D[i][j] = P[i][j];
			} else if (Q[i][j] > P[i][j] && Q[i][j] > 0) {
				D[i][j] = Q[i][j];
			} else {
				D[i][j] = 0;
			}

			if (D[i][j] > max)
			{
				max = D[i][j];
				occur = 1;
			} else if (D[i][j] == max){
				occur++;
			}
		}
	}

	printMatrix(P, D, Q, lenA, lenB);

	//traceback phase
	printf("--- Results ---\n");
	for (int i = 1; i < lenA + 1; ++i) {
		for(int j = 1; j < lenB + 1; ++j){
			if (D[i][j] == max) {
				printf("alignment n. %d\n", occur);
				traceback(P, D, Q, i, j, seqA, seqB);
				printf("\n");
				occur--;
			}
		}

		if(occur == 0) break;
	}
	
	return 0;
}

//calculate g(k) function
int gap_penalty(int k, int gap_opening, int enlargement) {
	return gap_opening + enlargement * k;
}

void traceback(Matrix P, Matrix D, Matrix Q, int i, int j, char * seqA, char * seqB) {
	int score = 0;
	alignment head = NULL;

	while(D[i][j] != 0){
		t_node * tmp = malloc(sizeof(t_node));
		tmp->baseA = seqA[i-1];
		tmp->baseB = seqB[j-1];
		tmp->match = (tmp->baseA != tmp->baseB) ? '|' : '*';
		tmp->next = head;
		head = tmp;

		if (i - 1 == 0) {
			i--;
			j--;
		} else if (((D[i-1][j-1] >= D[i-1][j]) && (D[i-1][j-1] >= D[i][j-1])) || (D[i-1][j-1] == 0)) {
			i--;
			j--;
		} else if (D[i-1][j] > D[i][j-1]) {
			i--;
		} else {
			j--;
		}
	}

	t_node * cursor = head;
	while(cursor != NULL) {
		printf("%c", cursor->baseA);
		cursor = cursor->next;
	}
	printf("\n");
	cursor = head;
	while(cursor != NULL) {
		printf("%c", cursor->match);
		cursor = cursor->next;
	}
	printf("\n");
	cursor = head;
	while(cursor != NULL) {
		printf("%c", cursor->baseB);
		cursor = cursor->next;
	}
	printf("\n");
}

////////////////////////////////////////////////////////////////////////////////////////////

void to_string(char * seqA, char * seqB, int ws, int wd, int gap_opening, int enlargement) {
	printf("\n+++++++++++++++++++++\n");
	printf("+ Sequence A: %s\n", seqA);
	printf("+ Sequence B: %s\n", seqB);
	printf("+ Match Score: %d\n", ws);
	printf("+ Mismatch Score: %d\n", wd);
	printf("+ Gap Opening: %d\n", gap_opening);
	printf("+ Enlargement: %d\n", enlargement);
	printf("+++++++++++++++++++++\n");
}
void printMatrix(Matrix P, Matrix D, Matrix Q, int lenA, int lenB){

	printf("--- P ---\n");
	for (int i = 0; i < lenA + 1; ++i) {
		for (int j = 0; j < lenB + 1; ++j) {
			if(i == 0 && j == 0)
				printf("\t \t");
			else if (i == 0)
				printf("\t -\u221e \t");
			else if (j == 0)
				printf("\t - \t");
			else 
				printf("\t %d \t", P[i][j]);
		}
		printf("\n");
	}

	printf("--- D ---\n");
	for (int i = 0; i < lenA + 1; ++i) {
		for (int j = 0; j < lenB + 1; ++j) {
			printf("\t %d \t", D[i][j]);
		}
		printf("\n");
	}

	printf("--- Q ---\n");
	for (int i = 0; i < lenA + 1; ++i) {
		for (int j = 0; j < lenB + 1; ++j) {
			if(i == 0 && j == 0)
				printf("\t \t");
			else if (j == 0)
				printf("\t -\u221e \t");
			else if (i == 0)
				printf("\t - \t");
			else 
				printf("\t %d \t", Q[i][j]);
		}
		printf("\n");
	}
}