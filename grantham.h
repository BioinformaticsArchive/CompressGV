/*
 * grantham.h
 *
 *  Created on: 24/03/2013
 *      Author: arran
 */

#ifndef GRANTHAM_H_
#define GRANTHAM_H_



#endif /* GRANTHAM_H_ */

#include <math.h>
#include <stddef.h>
#include "uthash.h"

#define GRANTHAM_PSO_SCALE 1e6
#define GRANTHAM_K_MAX 10
#define GRANTHAM_K_ONLY 1

typedef char aa_t;

typedef enum {
	c, p, v, k
} prop_t;

typedef struct {
	aa_t aa;
	double properties[3];
	UT_hash_handle hh;
} aaProp_t;

typedef struct {
	unsigned int pos;
	aa_t wt;
	aa_t variant;
	aa_t *msa;
} variant_t;

typedef struct {
	int no_of_species;
	int length;
	aa_t *acids;
} msa_t;

typedef struct {
	aaProp_t *properties;
	msa_t *msa;
	variant_t *variants[2];
	int nVariants[2];
	double coeff[4];
} granthamParam_t;

msa_t granthamMSA;
variant_t *granthamVariants[3];
int granthamNumVariants[3];
aaProp_t *granthamAAProperties;

aaProp_t* getAcidProperties(){
	aaProp_t *hash = NULL, *props;

	props = (aaProp_t*) malloc(20*sizeof(aaProp_t));
	memset(props, 0, 20*sizeof(aaProp_t));
	props[0] = (aaProp_t) {'A', {0, 8.1, 31}};
	props[1] = (aaProp_t) {'C', {2.75, 5.5, 55}};
	props[2] = (aaProp_t) {'D', {1.38, 13, 54}};
	props[3] = (aaProp_t) {'E', {0.92, 12.3, 83}};
	props[4] = (aaProp_t) {'F', {0, 5.2, 132}};
	props[5] = (aaProp_t) {'G', {0.74, 9, 3}};
	props[6] = (aaProp_t) {'H', {0.58, 10.4, 96}};
	props[7] = (aaProp_t) {'I', {0, 5.2, 111}};
	props[8] = (aaProp_t) {'K', {0.33, 11.3, 119}};
	props[9] = (aaProp_t) {'L', {0, 4.9, 111}};
	props[10] = (aaProp_t) {'M', {0, 5.7, 105}};
	props[11] = (aaProp_t) {'N', {1.33, 11.6, 56}};
	props[12] = (aaProp_t) {'P', {0.39, 8, 32.5}};
	props[13] = (aaProp_t) {'Q', {0.89, 10.5, 85}};
	props[14] = (aaProp_t) {'R', {0.65, 10.5, 124}};
	props[15] = (aaProp_t) {'S', {1.42, 9.2, 32}};
	props[16] = (aaProp_t) {'T', {0.71, 8.6, 61}};
	props[17] = (aaProp_t) {'V', {0, 5.9, 84}};
	props[18] = (aaProp_t) {'W', {0.13, 5.4, 170}};
	props[19] = (aaProp_t) {'Y', {0.2, 6.2, 136}};

	int i;
	for(i=0; i<20; i++){
		HASH_ADD(hh, hash, aa, sizeof(aa_t), &(props[i]));
	}

	return hash;
}

aaProp_t* granthamInit(){
	return getAcidProperties();
}

void granthamFree(aaProp_t *props){
	// See the uthash User Guide for more information - http://troydhanson.github.io/uthash/userguide.html#_delete_item
	aaProp_t *currProp, *tmp;
	HASH_ITER(hh, props, currProp, tmp){
		HASH_DEL(props, currProp);
		free(currProp);
	}
}

#ifndef fmin
	#define fmin(a,b) ((a)<(b) ? (a) : (b))
#endif

#ifndef fmax
	#define fmax(a,b) ((a)>(b) ? (a) : (b))
#endif

double gv(aa_t *acids, unsigned int n, double *coeff){
	double minP, maxP, gv = 0.;
	int i;
	aaProp_t *aaProp;

	for(i=0; i<3; i++){
		minP = 999.;
		maxP = 0.;
		int a;
		for(a=0; a<n; a++){
			HASH_FIND(hh, granthamAAProperties, &(acids[a]), sizeof(aa_t), aaProp);
			if(aaProp!=NULL){
				double prop = aaProp->properties[i];
				minP = fmin(minP, prop);
				maxP = fmax(maxP, prop);
			}
			else if(acids[a]!='-' && acids[a]!='X'){
				fprintf(stderr, "Could not find properties for amino acid: %c\n", acids[a]);
				return -1.;
			}
		}
		gv += coeff[i] * pow(maxP-minP, 2);
	}
	return sqrt(gv);
}

void granthamCoefficients(double *coeff){
	coeff[0] = 1.833;
	coeff[1] = 0.1018;
	coeff[2] = 0.000399;
	coeff[3] = 1.;
}

bool getMSA(FILE *fp, msa_t *msa){
	msa->length = 0;
	int c, i=0, j=0, mem=512;
	bool *tmpB, *is_aa_pos = (bool*) malloc(mem*sizeof(bool));

	// Determine the length of the human protein
	// Note which columns of the text file denote an amino acid in the human protein
	while((c=getc(fp))!=EOF){
		if(i+1==mem){
			mem += 256;
			tmpB = realloc(is_aa_pos, mem);
			if(tmpB!=NULL){
				is_aa_pos = tmpB;
			}
			else {
				fprintf(stderr, "Error allocating memory for first line of MSA\n");
				free(is_aa_pos);
				return false;
			}
		}

		if(c!='-' && c!='\n'){
			msa->length++;
			is_aa_pos[i] = true;
		}
		else {
			is_aa_pos[i] = false;
		}
		i++;
	}
	rewind(fp);

	// Parse the text file for the MSA
	msa->acids = (aa_t*) malloc(msa->length*sizeof(aa_t));
	msa->no_of_species = 1;
	i = j = 0;
	aa_t* tmpA;
	bool newLine = false;

	while((c=getc(fp))!=EOF){
		if(!(c>=65 && c<=89) && c!=10 && c!=45){
			fprintf(stderr, "Invalid character in MSA line %i pos %i: '%c'\n", msa->no_of_species, i+1, c);
			free(is_aa_pos);
			// Not freeing acids memory, it is expected that it will be done by the function creating the msa struct
			return false;
		}
		if(c=='\n'){
			if(newLine){
				continue;
			}
			newLine = true;
			msa->no_of_species++;
			i = j = 0;
			tmpA = realloc(msa->acids, msa->length*msa->no_of_species*sizeof(aa_t));
			if(tmpA!=NULL){
				msa->acids = tmpA;
			}
			else {
				fprintf(stderr, "Error allocating memory for MSA species number %i\n", msa->no_of_species);
				free(is_aa_pos);
				// Not freeing acids memory, it is expected that it will be done by the function creating the msa struct
				return false;
			}
			continue;
		}

		newLine = false;

		if(is_aa_pos[i]){
			msa->acids[(msa->no_of_species-1)*msa->length+j] = c;
			j++;
		}
		i++;
	}

	// One too many lines have been read; memory usage is insignificant so no point reducing the size of msa->acids
	if(newLine){
		msa->no_of_species--;
	}

	free(is_aa_pos);
	return true;
}

#define STATE_WT 0
#define STATE_POS 1
#define STATE_VAR 2

int getVariants(FILE *fp, variant_t** varsPtr, msa_t *msa, bool canBeEmpty, char* file_name){
	int n=0, state=STATE_WT, digits=0, pos=0;
	variant_t *vars, *tmp, *currVar;

	char c;
	while((c=getc(fp))!=EOF){
		pos++;
		if(!(c>=65 && c<=89) && c!=10 && !(c>=48 && c<=57)){
			fprintf(stderr, "Invalid character in variant file [%s] line %i pos %i: '%c'\n", file_name, n + (state==STATE_WT), pos, c);
			return -1;
		}

		switch(state){
		case STATE_WT:
			if(c=='\n'){
				break;
			}

			n++;

			if(!(c>=65 && c<=89)){
				fprintf(stderr, "Invalid syntax in variant file [%s] line %i: no wild-type amino acid found\n", file_name, n);
				return -1;
			}

			if(n==1){
				vars = (variant_t*) malloc(sizeof(variant_t));
			}
			else {
				tmp = realloc(vars, n*sizeof(variant_t));
				if(tmp!=NULL){
					vars = tmp;
				}
				else {
					fprintf(stderr, "Error allocating memory for variant number %i in file '%s'\n", n, file_name);
					return -1;
				}
			}

			currVar = &(vars[n-1]);
			currVar->wt = c;
			currVar->pos = 0;

			state = STATE_POS;
			digits = 0;

			break;
		case STATE_POS:
			if((c>=48 && c<=57)){
				digits++;
				currVar->pos = 10*currVar->pos + c - 48;
			}
			else if((c>=65 && c<=89)){
				currVar->pos--; //Biologists don't know about array indexing so start at 1
				if(!digits){
					fprintf(stderr, "Invalid syntax in variant file [%s] line %i: no position found\n", file_name, n);
					return -1;
				}
				else if(msa->acids[currVar->pos]!=currVar->wt){
					fprintf(stderr, "Error in variant file [%s] line %i: wild-type amino acid is '%c' but MSA at position %i is '%c'\n", file_name, n, currVar->wt, currVar->pos+1, msa->acids[currVar->pos]);
					return -1;
				}

				int s;
				currVar->msa = (aa_t*) malloc(msa->no_of_species * sizeof(aa_t));
				for(s=0; s<msa->no_of_species; s++){
					currVar->msa[s] = msa->acids[s*msa->length + currVar->pos];
				}

				currVar->variant = c;
				state = STATE_VAR;
			}
			else {
				fprintf(stderr, "Invalid syntax in variant file [%s] line %i: encountered invalid character ('%c') while parsing position\n", file_name, n, c);
				return -1;
			}
			break;
		case STATE_VAR:
			if(c!=10){
				fprintf(stderr, "Invalid syntax in variant file [%s] line %i: extra character ('%c') after variant amino acid\n", file_name, n, c);
				return -1;
			}
			state = STATE_WT;
			pos = 0;
			break;
		}
	}

	if(!n && !canBeEmpty){
		fprintf(stderr, "Variant file [%s] is empty\n", file_name);
		return -1;
	}

	*varsPtr = vars;
	return n;
}

double granthamMetric(variant_t *var, double *coeff){
	aa_t snp[2] = {var->wt, var->variant};
	return gv(&(snp[0]), 2, coeff) * pow(coeff[3], -gv(var->msa, granthamMSA.no_of_species, coeff));
}

double granthamCluster(double *coeff, int returnCutoff, double *returnAll){
	double mean[2] = {0., 0.};
	double sd2[2] = {0., 0.};
	double metric, delta;

	int i, j;
	for(i=0; i<2; i++){
		for(j=0; j<granthamNumVariants[i]; j++){
			metric = granthamMetric(granthamVariants[i] + j, coeff);

			//See http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
			delta = metric - mean[i];
			mean[i] += delta / (j+1);
			sd2[i] += delta * (metric-mean[i]) / (granthamNumVariants[i]-1);
		}
	}

	if(returnAll!=NULL){
		returnAll[0] = mean[0];
		returnAll[1] = mean[1];
		returnAll[2] = sd2[0];
		returnAll[3] = sd2[1];
	}

	if(returnCutoff){
		return (mean[0] + mean[1])/2;
	}
	else {
		double denom = sd2[0] + sd2[1];
		return denom ? (mean[0] - mean[1])/(sqrt(sd2[0]) + sqrt(sd2[1])) : -1;
	}
}

#define GRANTHAM_K_SCALE (GRANTHAM_K_MAX-1)*c[3]/GRANTHAM_PSO_SCALE+1

#if GRANTHAM_K_ONLY==1
#define GRANTHAM_SCALED {1.833, 0.1018, 0.000399, GRANTHAM_K_SCALE}
#else
#define GRANTHAM_SCALED {c[0]/GRANTHAM_PSO_SCALE, c[1]/GRANTHAM_PSO_SCALE, c[2]/GRANTHAM_PSO_SCALE, GRANTHAM_K_SCALE}
#endif

double granthamPSO(double *c, int dim, void *params) {
	double coeff[4] = GRANTHAM_SCALED;
	return -granthamCluster(&(coeff[0]), false, NULL);
}

bool granthamClassify(variant_t *var, double *coeff, double *outcome){
	double metric = granthamMetric(var, coeff);
	double cutoff = granthamCluster(coeff, true, outcome==NULL ? NULL : outcome + 2);
	if(outcome!=NULL){
		outcome[0] = metric;
		outcome[1] = cutoff;
	}
	return metric > cutoff;
}
