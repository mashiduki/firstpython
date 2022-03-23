#include <stdio.h>
#include <stdlib.h>
#include<string.h>

#define nb 3	// number of blocks in one direction
#define nb2 12	// number of small blocks in one direction
#define MC 4000

double fnorm[nb2][nb2];
double fnorm_tmp[nb2][nb2];

double evaluate_off_norm(int *perm, double *max_off_norm, int *max_perm)
{
	int i, j, k, lb2;
	double off_norm = 0.0;
	for (i=0; i<nb2; i++)
		for (j=0; j<nb2; j++)
			fnorm_tmp[perm[i]][perm[j]] = fnorm[i][j];

	lb2 = nb2 / nb;	// number of small blocks in one block
	for (k=0; k<nb; k++)
 		for (i=k*lb2; i<(k+1)*lb2; i++)
			for (j=k*lb2; j<(k+1)*lb2; j++)
				off_norm += fnorm_tmp[i][j];
//	printf("off_norm = %lf\n", off_norm);
	if (off_norm > *max_off_norm)
	{
		*max_off_norm = off_norm;
		printf("current max_off_norm = %lf, ", *max_off_norm);
		printf("perm = ");
		for (i=0; i<nb2; i++)
		{
			max_perm[i] = perm[i];
			printf("%d ", max_perm[i]);
		}
		printf("\n");
	}
	return 0;
}

int generate_permutation(int level, int *perm, double *max_off_norm,
			 int *max_perm)
{
	int i, j, icall;

	for (i=0; i<nb2; i++)
	{
		icall = 1;
		for (j=0; j<level; j++)
			if (perm[j] == i) icall = 0;
		if (icall == 1)
		{
			perm[level] = i;
			if (level == nb2-1)
			{
//				for (j=0; j<nb2; j++)
//					printf("%d ", perm[j]);
//				printf("\n");
				evaluate_off_norm(perm, max_off_norm, max_perm);
			}
			else
			{
				generate_permutation(level+1, perm, 
						     max_off_norm, max_perm);
			}
		}
	}
	return 0;
}

int main(void)
{
	int perm[nb2], max_perm[nb2];
	int level;
	double max_off_norm = 0.0;
	
	char *fmm="f.csv";
    FILE *csv;
    char bff[MC];
    char *c;
    int i,j;

    csv=fopen(fmm,"r");
    i=j=0;
    while(fgets(bff,MC,csv)!=NULL){
        c=bff;
        while((c=strtok(c,",\n"))!=NULL){
             fnorm[i][j++]=strtod(c,NULL);
            c=NULL;
        }
    i++;
    j=0;
    }
    // printf(" fnorm[%d][%d]:\n", N,N);
    // for(i=0;i<N;i++){
    //     for(j=0;j<N;j++){
    //         printf("%9f", fnorm[i][j]);
    //     }
    //     printf("\n");
    // }

	level = 0;
	generate_permutation(level, perm, &max_off_norm, max_perm);

    for (i=0; i<nb2; i++)
        for (j=0; j<nb2; j++)
            fnorm_tmp[max_perm[i]][max_perm[j]] = fnorm[i][j];
    printf(" fnorm[%d][%d]:\n", nb2,nb2);
    for(i=0;i<nb2;i++){
         for(j=0;j<nb2;j++){
             printf("%5.0f ", fnorm_tmp[i][j]);
         }
         printf("\n");
     }
    double sum=0;
    int blocksize=nb2/nb; //対角ブロック
     for(i=0;i<nb;i++){
         for(j=0;j<blocksize;j++){
             for(int k=0;k<blocksize;k++){
                 sum+=fnorm_tmp[j+(i*(blocksize))][k+(i*(blocksize))]*fnorm_tmp[j+(i*(blocksize))][k+(i*(blocksize))];
             }
         }
     }
    printf("対角ブロックの和:%5.0f\n",sum);
    int fclose(FILE *f);
	return 0;
}
