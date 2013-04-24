/*
 * main.c
 *
 *  Created on: 23/03/2013
 *      Author: arran
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "grantham.h"
#include "main.h"

int main(int argc, char *argv[]){
	if(argc!=5){
		fprintf(stderr, "Usage: ./grantham <msa> <deleterious> <neutral> <classify>\n");
		return 1;
	}

	FILE *fp[4];
	msa_t msa;
	variant_t *variants[3];
	int nVariants[3];
	float coeff[5];

	aaProp_t *props = granthamInit();
	granthamCoefficients(&(coeff[0]));

	int f;
	for(f=0; f<4; f++){
		if(!(fp[f]=fopen(argv[f+1], "r"))){
			fprintf(stderr, "Could not open file: %s\n", argv[f+1]);
			closeFiles(fp, f);
			return 1;
		}
	}

	getMSA(fp[0], &msa);

	int v;
	for(v=0; v<3; v++){
		nVariants[v] = getVariants(fp[v+1], &(variants[v]), &msa, v==2, argv[v+2]);
	}

	fprintf(stderr, "%f\n", granthamMetric(variants[0], props, &(coeff[0]), &msa));

	closeFiles(fp, 4);
	free(msa.acids);
	//granthamFree(props);
	for(v=0; v<3; v++){
		if(nVariants[v]){
			free(variants[v]);
		}
	}
	return 0;
}

void closeFiles(FILE *fp[], int f){
	int i;
	for(i=0; i<f; i++){
		fclose(fp[i]);
	}
}
