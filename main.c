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



	matthews_t *matt = assessModel();
	fprintf(stdout, "TP: %i	TN: %i	FP: %i	FN: %i\n", matt->tp, matt->tn, matt->fp, matt->fn);
	fprintf(stdout, "MCC: %f\n", matt->coefficient);
	fprintf(stdout, "Chi-squared: %f\n", matt->chiSquare);
	fprintf(stdout, "Sens.: %f	Spec.: %f\n", (double) matt->tp / (matt->tp + matt->fn), (double) matt->tn / (matt->tn + matt->fp));
	//free(matt);



	//closeFiles(fp, 4);
	//free(granthamMSA.acids);
	//free(c);
	//granthamFree(granthamAAProperties);
	for(v=0; v<3; v++){
		if(granthamNumVariants[v]){
			//free(granthamVariants[v]);
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

	settings.size = 20;
	settings.steps = 20;

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

// Iterate through all known variants, optimise the model with the variant excluded, and classify the variant
// Return the Matthews correlation coefficient
matthews_t* assessModel(){
	int i, j, r;
	matthews_t *matt = (matthews_t*) calloc(sizeof(matthews_t), 0);
	variant_t *classify = (variant_t*) malloc(sizeof(variant_t));
	double *outcome = (double*) malloc(6*sizeof(double));

	for(i=0; i<2; i++){
		granthamNumVariants[i]--; //exclude the last variant in this set

		for(j=0; j<granthamNumVariants[i]+1; j++){
			double *c = optimiseCoefficients();
			double scaled[] = GRANTHAM_SCALED;
			memcpy(classify, granthamVariants[i] + granthamNumVariants[i], sizeof(variant_t)); //note that granthamNumVariants[i] has already been decremented so this passes the excluded variant

			fprintf(stderr, "%c%i%c	", classify->wt, classify->pos+1, classify->variant);

			if(granthamClassify(classify, &(scaled[0]), outcome)){
				if(!i){
					matt->tp++;
					fprintf(stderr, "TP");
				}
				else {
					matt->fp++;
					fprintf(stderr, "FP");
				}
			}
			else {
				if(!i){
					matt->fn++;
					fprintf(stderr, "FN");
				}
				else {
					matt->tn++;
					fprintf(stderr, "TN");
				}
			}

			fprintf(stderr, "	%f %c %f	%f	%f	%f	%f	%i	%i	%f\n", outcome[0], outcome[0]>outcome[1] ? '>' : '<', outcome[1], outcome[2], outcome[3], outcome[4], outcome[5], granthamNumVariants[0], granthamNumVariants[1], scaled[3]);

			// Rotate the variants so the next one is excluded
			for(r=granthamNumVariants[i]-1; r>=0; r--){
				memcpy(granthamVariants[i] + r + 1, granthamVariants[i] + r, sizeof(variant_t));
			}
			memcpy(granthamVariants[i], classify, sizeof(variant_t));

			free(c);
		}

		granthamNumVariants[i]++;
	}

	free(classify);
	free(outcome);

#define NO_ZERO(a,b) (a==0 && b==0 ? 1 : a+b)

	matt->coefficient = (matt->tp*matt->tn - matt->fp*matt->fn) / sqrt(NO_ZERO(matt->tp, matt->fp) * NO_ZERO(matt->tp, matt->fn) * NO_ZERO(matt->tn, matt->fp) * NO_ZERO(matt->tn, matt->fn));
	matt->chiSquare = pow(matt->coefficient, 2) * (granthamNumVariants[0] + granthamNumVariants[1]);
	return matt;
}
