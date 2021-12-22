/**************************************************************************************
* 
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2020/21
* 
* Progetto dell'algoritmo di Regressione
* in linguaggio assembly x86-32 + SSE
* 
* Fabrizio Angiulli, aprile 2019
* 
**************************************************************************************/

/*
* 
* Software necessario per l'esecuzione:
* 
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
* 
* entrambi sono disponibili come pacchetti software 
* installabili mediante il packaging tool del sistema 
* operativo; per esempio, su Ubuntu, mediante i comandi:
* 
*    sudo apt-get install nasm
*    sudo apt-get install gcc
* 
* potrebbe essere necessario installare le seguenti librerie:
* 
*    sudo apt-get install lib32gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
* 
* Per generare il file eseguibile:
* 
* nasm -f elf64 fss64.nasm && gcc -m64 -msse -O0 -no-pie sseutils64.o fss64.o fss46c.c -o fss64c -lm && ./fss64c $pars
* 
* oppure
* 
* ./runfss64
* 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

#define	ftype		float
#define	itype		int
#define	MATRIX		ftype*
#define	VECTOR		ftype*














































































































































































































































/*
* 
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate 
* 	mediante un array (double*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere 
* 	memorizzate mediante array di array (double**).
* 
* 	In entrambi i casi il candidato dovrà inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente è che le matrici siano in row-major order.
* 
*/

void* get_block(int size, int elements) { 
	return _mm_malloc(elements*size,16);
}

void free_block(void* p) { 
	_mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(ftype),rows*cols);
}

void dealloc_matrix(void* mat) {
	free_block(mat);
}

/*
* 
* 	load_data
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*8 byte: matrix data in row-major order --> numeri floating-point a precisione doppia
* 
*****************************************************************************
*	Se lo si ritiene opportuno, è possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
MATRIX load_data(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}

	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);
	
	MATRIX data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(ftype), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*k = cols;
	
	return data;
}


int exists(char* fname){
	FILE* f;
	if(f=fopen(fname, "rb")){
		fclose(f);
		return 1;
	}
	return 0;
}

int main(int argc, char** argv) {

	char fname[256];
	int silent = 0;

	char* xhfilename;
	char* solfilename;
	int np;
	int dim;
	clock_t t;
	double err;

	//
	// Visualizza la sintassi del passaggio dei parametri da riga comandi
	//

	if(argc <= 1){
		printf("%s -xh <xh> -sol <sol> -np <np> [-s]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\txh: file ds2 contenente il vettore calcolato\n");
		printf("\tsol: file ds2 contenente il vettore corretto\n");
		printf("\tnp: il numero di pesci, default 25\n");
		printf("\t-s: silenzioso\n");
		exit(0);
	}

	//
	// Legge i valori dei parametri da riga comandi
	//
	int par = 1;
	while (par < argc) {
		if (strcmp(argv[par],"-s") == 0) {
			silent = 1;
			par++;
		} else if (strcmp(argv[par],"-xh") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing computed result file name!\n");
				exit(1);
			}
			xhfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-sol") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing correct result file name!\n");
				exit(1);
			}
			solfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-np") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing np value!\n");
				exit(1);
			}
			np = atoi(argv[par]);
			par++;
		} else{
			printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
			par++;
		}
	}

	//
	// Legge i dati e verifica la correttezza dei parametri
	//
	if(solfilename == NULL || strlen(solfilename) == 0){
		printf("Missing correct result file name!\n");
		exit(1);
	}
	if(xhfilename == NULL || strlen(xhfilename) == 0){
		printf("Missing computed result file name!\n");
		exit(1);
	}
	int x;
	MATRIX sol = load_data(solfilename, &x, &dim);
	int d;
	if(!exists(xhfilename)){
		err = -1;
		if(!silent)
			printf("File %s does not exist!\n", xhfilename);
	}
	else{
		MATRIX xh = load_data(xhfilename, &x, &d);
		if(x>1 || d != dim){
			err = -1;
			if(!silent)
				printf("File %s has wrong size, found %dx%d, expected 1x%d!\n", xhfilename, x, d, dim);
		}
		else{
			err = 0;
			for(int i=0; i<dim; i++)
				err = err + (xh[i] - sol[i])*(xh[i] - sol[i]);
			err = sqrt(err);
		}
	}
exit:
	printf("%f\n", err);
	dealloc_matrix(sol);
	return 0;
}
