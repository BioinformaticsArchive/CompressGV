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

typedef char aa_t;

typedef enum {
	c, p, v
} prop_t;

typedef struct {
	aa_t aa;
	float properties[3];
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
		HASH_ADD_INT(hash, aa, &(props[i]));
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

float gv(aaProp_t *properties, aa_t *acids, unsigned int n, float *coeff){
	float minP, maxP, gv = 0.;
	int i;
	for(i=0; i<3; i++){
		minP = 999.;
		maxP = 0.;
		int a;
		for(a=0; a<n; a++){
			float prop = properties[acids[a]].properties[i];
			minP = fmin(minP, prop);
			maxP = fmax(maxP, prop);
		}
		gv += coeff[i] * pow(maxP-minP, 2);
	}
	return coeff[3] * sqrt(gv);
}

void granthamCoefficients(float *coeff){
	coeff[0] = 1.833;
	coeff[1] = 0.1018;
	coeff[2] = 0.000399;
	coeff[3] = 50.723;
}
