/*
 * main.h
 *
 *  Created on: 23/03/2013
 *      Author: arran
 
---------------------------------------------------------------------------------------

Copyright 2013 Arran Schlosberg.

This file is part of https://github.com/aschlosberg/CompressGV (CompressGV)

    CompressGV is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CompressGV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CompressGV. If not, see <http://www.gnu.org/licenses/>.

---------------------------------------------------------------------------------------
 
 */

#ifndef MAIN_H_
#define MAIN_H_



#endif /* MAIN_H_ */

typedef struct {
	int tp, tn, fp, fn; // can not be unsigned otherwise calculation of negative coefficients fail
	double coefficient, chiSquare;
} matthews_t;

void closeFiles(FILE *fp[], int f);

double* optimiseCoefficients();

matthews_t* assessModel();

int runTests(){
	double *coeff = (double*) malloc(4*sizeof(double));
	granthamCoefficients(coeff);

	//from Grantham (1974)
	int known[] = {110,145,74,58,99,124,56,142,155,144,112,89,68,46,121,65,80,135,177,102,103,71,112,96,125,97,97,77,180,29,43,86,26,96,54,91,101,98,92,96,32,138,5,22,36,198,99,113,153,107,172,138,15,61,38,27,68,42,95,114,110,169,77,76,91,103,108,93,87,147,58,69,59,89,103,92,149,47,42,65,78,85,65,81,128,64,60,94,113,112,195,86,91,111,106,126,107,84,148,109,29,50,55,192,84,96,133,97,152,121,21,88,135,153,147,159,98,87,80,127,94,98,127,184,21,33,198,94,109,149,102,168,134,10,61,22,205,100,116,158,102,177,140,28,40,194,83,99,143,85,160,122,36,37,174,154,139,202,154,170,196,215,24,68,32,81,40,87,115,46,53,61,29,101,130,94,23,42,142,174,101,56,95,110,45,160,181,126,152,67};
	aa_t cols[] = {'R','L','P','T','A','V','G','I','F','Y','C','H','Q','N','K','D','E','M','W'};
	aa_t rows[] = {'S','R','L','P','T','A','V','G','I','F','Y','C','H','Q','N','K','D','E','M'};

	int c, r, i=0, err=0;
	for(r=0; r<19; r++){
		for(c=r; c<19; c++){
			aa_t acids[] = {cols[c], rows[r]};
			double g = 50.723 * gv(&(acids[0]), 2, coeff);
			if(floor(g)!=known[i] && ceil(g)!=known[i]){
				fprintf(stderr, "r: %i (%i)	c: %i (%i)	known: %i	calc: %f\n", r, rows[r], c, cols[c], known[i], g);
				err++;
			}
			i++;
		}
	}
	fprintf(stderr, "%i errors in %i checks\n", err, i);

	free(coeff);
	return err;
}
