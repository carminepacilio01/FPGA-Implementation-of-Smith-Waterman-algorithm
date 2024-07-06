// Implementation of the Smith-Waterman (Gotoh)      
//
//  Author: Carmine Pacilio							 
//	Year: 2024										  
//	Version: 1.0									  
////////////////////////////////////////////////////
#include <limits.h>

void smithWaterman(int **P, int **D, int **Q, int wd, int ws, 
	int gap_opening, int enlargement, int lenA, char *seqA, int lenB, char *seqB);
int gap_penalty(int k, int gap_opening, int enlargement);
//void traceback(char **P, char **D, char **Q, int i, int j, char *seqA, char *seqB);

void smithWaterman(int **P, int **D, int **Q, int wd, int ws, 
	int gap_opening, int enlargement, int lenA, char *seqA, int lenB, char *seqB){

    int i, j;

    for (i = 0; i < lenA + 1; ++i) {
		P[i][0] = INT_MIN / 2;
		Q[i][0] = INT_MIN / 2;
		D[i][0] = 0;
	}
    for (i = 0; i < lenB + 1; ++i) {
		P[0][i] = INT_MIN / 2;
		Q[0][i] = INT_MIN / 2;
		D[0][i] = 0;
	}
	D[0][0] = 0;

    //Filling matrix phase
    for (i = 1; i < lenA + 1; ++i) {
        for (j = 1; j < lenB + 1; ++j) {
            //Calculate P = min{Di-1,j + g(1); Pi-1,j}
            P[i][j] = (D[i - 1][j] + gap_penalty(1, gap_opening, enlargement) > P[i - 1][j] + enlargement) ?
                             D[i - 1][j] + gap_penalty(1, gap_opening, enlargement) : P[i - 1][j] + enlargement;

            //Calculate Q = min{Di,j-1 + g(1); Qi,j-1}
            Q[i][j] = (D[i][j - 1] + gap_penalty(1, gap_opening, enlargement) > Q[i][j - 1] + enlargement) ?
                             D[i][j - 1] + gap_penalty(1, gap_opening, enlargement) : Q[i][j - 1] + enlargement;

            //Calculate D = min{Di-1,j-1 + w(ai, bj); Pi,j; Qij}
            int tmp = D[i - 1][j - 1];
            tmp += (seqA[i - 1] == seqB[j - 1]) ? wd : ws;

            if (tmp > P[i][j] && tmp > Q[i][j] && tmp > 0) {
                D[i][j] = tmp;
            } else if (P[i][j] > Q[i][j] && P[i][j] > 0) {
                D[i][j] = P[i][j];
            } else if (Q[i][j] > P[i][j] && Q[i][j] > 0) {
                D[i][j] = Q[i][j];
            } else {
                D[i][j] = 0;
            }
        }
    }
}

//calculate g(k) function
int gap_penalty(int k, int gap_opening, int enlargement) {
	return gap_opening + enlargement * k;
}

/*void traceback(char **P, char **D, char **Q, int i, int j, char *seqA, char *seqB) {
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
}*/
