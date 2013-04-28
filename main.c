/*
 * main.c
 *
 *  Created on: 23/03/2013
 *      Author: arran
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "pso/pso.h"
#include "pso/pso.c"
#include "grantham.h"
#include "main.h"

int main(int argc, char *argv[]){
	if(argc!=5){
		fprintf(stderr, "Usage: ./grantham <msa> <deleterious> <neutral> <classify>\n");
		return 1;
	}

	FILE *fp[4];

	granthamAAProperties = granthamInit();

	int f;
	for(f=0; f<4; f++){
		if(!(fp[f]=fopen(argv[f+1], "r"))){
			fprintf(stderr, "Could not open file: %s\n", argv[f+1]);
			closeFiles(fp, f);
			return 1;
		}
	}

	getMSA(fp[0], &granthamMSA);

	int v;
	for(v=0; v<3; v++){
		granthamNumVariants[v] = getVariants(fp[v+1], &(granthamVariants[v]), &granthamMSA, v==2, argv[v+2]);
	}

	double *c = optimiseCoefficients();
	double solution[] = GRANTHAM_SCALED;
	fprintf(stderr, "Clustering: %f\n", granthamCluster(&(solution[0])));
	fprintf(stderr, "c: %f	p: %f	v: %f	k: %f\n", solution[0], solution[1], solution[2], solution[3]);

	closeFiles(fp, 4);
	free(granthamMSA.acids);
	free(c);
	//granthamFree(granthamAAProperties);
	for(v=0; v<3; v++){
		if(granthamNumVariants[v]){
			free(granthamVariants[v]);
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

double* optimiseCoefficients(){
	pso_obj_fun_t obj_fun = granthamPSO;
	pso_settings_t settings;
	pso_set_default_settings(&settings);

	settings.size = 30;
	settings.steps = 30;

	settings.dim = 4;
	settings.nhood_strategy = PSO_NHOOD_RANDOM;
	settings.nhood_size = 10;
	settings.w_strategy = PSO_W_CONST;
	settings.x_lo = 0;
	settings.x_hi = GRANTHAM_PSO_SCALE;
	settings.goal = -99999; //deliberately unachievable
	// Seed with a constant to allow reproducibility of results
	settings.seed = 314159265;
	settings.print_every = 0;

	pso_result_t solution;
	solution.gbest = malloc(settings.dim * sizeof(double));

	pso_solve(obj_fun, NULL, &solution, &settings);

	return solution.gbest;
}
