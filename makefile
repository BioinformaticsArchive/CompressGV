MKDIR = mkdir -p

all:
	${MKDIR} ./bin
	gcc -ggdb -o ./bin/grantham main.c -lm -lgsl -lgslcblas -lz
