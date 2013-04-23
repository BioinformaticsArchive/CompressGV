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

	aaProp_t *props = granthamInit();

	int f;
	for(f=0; f<4; f++){
		if(!(fp[f]=fopen(argv[f+1], "r"))){
			fprintf(stderr, "Could not open file: %s\n", argv[f+1]);
			closeFiles(fp, f);
			return 1;
		}
	}

	getMSA(fp[0], &msa);

	int i, j;
	for(i=0; i<msa.no_of_species; i++){
		for(j=0; j<msa.length; j++){
			fprintf(stderr, "%c", msa.acids[i*msa.length+j]);
		}
		fprintf(stderr, "\n\n");
	}

	closeFiles(fp, 4);
	free(msa.acids);
	//granthamFree(props);
	return 0;
}

void closeFiles(FILE *fp[], int f){
	int i;
	for(i=0; i<f; i++){
		fclose(fp[i]);
	}
}

bool getMSA(FILE *fp, msa_t *msa){
	msa->length = 0;
	int c, i=0, j=0, mem=512;
	bool *tmpB, *is_aa_pos = (bool*) malloc(mem*sizeof(bool));

	// Determine the length of the human protein
	// Note which columns of the text file denote an amino acid in the human protein
	do {
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

		c = getc(fp);
		if(c!='-' && c!='\n'){
			msa->length++;
			is_aa_pos[i] = true;
		}
		else {
			is_aa_pos[i] = false;
		}
		i++;
	}
	while(c!='\n');
	rewind(fp);

	// Parse the text file for the MSA
	msa->acids = (aa_t*) malloc(msa->length*sizeof(aa_t));
	msa->no_of_species = 1;
	i = j = 0;
	aa_t* tmpA;
	bool newLine = false;

	do {
		c = getc(fp);
		if((c<65 || c>90) && c!=10 && c!=45){
			continue;
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
				// Not freeing memory, it is expected that it will be done by the function creating the msa struct
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
	while(c!=EOF);

	// One too many lines have been read; memory usage is insignificant so no point reducing the size of msa->acids
	if(newLine){
		msa->no_of_species--;
	}

	free(is_aa_pos);
	return true;
}
