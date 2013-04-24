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

	int v;
	for(v=0; v<3; v++){
		nVariants[v] = getVariants(fp[v+1], &(variants[v]), &msa, v==2, argv[v+2]);
	}

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

int getVariants(FILE *fp, variant_t** vars, msa_t *msa, bool canBeEmpty, char* file_name){
	int n=0, state=STATE_WT, digits=0, pos=0;
	variant_t *tmp, currVar;

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
				*vars = (variant_t*) malloc(sizeof(variant_t));
			}
			else {
				tmp = realloc(*vars, n*sizeof(variant_t));
				if(tmp!=NULL){
					*vars = tmp;
				}
				else {
					fprintf(stderr, "Error allocating memory for variant number %i in file '%s'\n", n, file_name);
					return -1;
				}
			}

			currVar = (*vars)[n-1];
			currVar.wt = c;
			currVar.pos = 0;

			state = STATE_POS;
			digits = 0;

			break;
		case STATE_POS:
			if((c>=48 && c<=57)){
				digits++;
				currVar.pos = 10*currVar.pos + c - 48;
			}
			else if((c>=65 && c<=89)){
				currVar.pos--; //Biologists don't know about array indexing so start at 1
				if(!digits){
					fprintf(stderr, "Invalid syntax in variant file [%s] line %i: no position found\n", file_name, n);
					return -1;
				}
				else if(msa->acids[currVar.pos]!=currVar.wt){
					fprintf(stderr, "Error in variant file [%s] line %i: wild-type amino acid is '%c' but MSA at position %i is '%c'\n", file_name, n, currVar.wt, currVar.pos+1, msa->acids[currVar.pos]);
					return -1;
				}

				int s;
				currVar.msa = (aa_t*) malloc(msa->no_of_species * sizeof(aa_t));
				for(s=0; s<msa->no_of_species; s++){
					currVar.msa[s] = msa->acids[s*msa->length + currVar.pos];
				}

				currVar.variant = c;
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

	return n;
}
