#ifndef SMITH_WATERMAN_H
#define SMITH_WATERMAN_H

void smithWaterman(int **P, int **Q, int **D, int wd, int ws, int gap_opening, int enlargement, int lenA, char *seqA, int lenB, char *seqB);
int gap_penalty(int k, int gap_opening, int enlargement);
//void traceback(char **P, char **D, char **Q, int i, int j, char *seqA, char *seqB);

#endif