/**************************************************************************************
* 
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2020/21
* 
* Progetto dell'algoritmo Fish School Search 221 231 a
* in linguaggio assembly x86-64 + SSE
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
*    sudo apt-get install lib64gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
* 
* Per generare il file eseguibile:
* 
* nasm -f elf64 fss64.nasm && gcc -m64 -msse -O0 -no-pie sseutils64.o fss64.o fss64c.c -o fss64c -lm && ./fss64c $pars
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

#define	type		double
#define	MATRIX		type*
#define	VECTOR		type*
#define EPSILON 0.000000001

typedef struct {
	MATRIX x; //posizione dei pesci
	VECTOR xh; //punto associato al minimo di f, soluzione del problema
	VECTOR c; //coefficienti della funzione
	VECTOR r; //numeri casuali
	int np; //numero di pesci, quadrato del parametro np
	int d; //numero di dimensioni del data set
	int iter; //numero di iterazioni
	type stepind; //parametro stepind
	type stepvol; //parametro stepvol
	type wscale; //parametro wscale
	int display;
	int silent;
} params;

typedef struct {
	VECTOR w;
	VECTOR deltaf;
	VECTOR deltax;
	VECTOR baricentro;
	int rand;
	type wbranco;
	type stepindIni; //parametro stepind
	type stepvolIni; //parametro stepvol
}var;

int np;
int d;
int iter;
int effettuato;
int allineamento=4;
int allineamentoPerfetto=0;
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
	return _mm_malloc(elements*size,32); 
}

void free_block(void* p) { 
	_mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(type),rows*cols);
}

void dealloc_matrix(MATRIX mat) {
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
* 	primi 4 byte: numero di righe (N) --> numero intero
* 	successivi 4 byte: numero di colonne (M) --> numero intero
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri doubleing-point a precisione singola
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
	status = fread(data, sizeof(type), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*k = cols;
	
	return data;
}

/*
* 	save_data
* 	=========
* 
*	Salva su file un array lineare in row-major order
*	come matrice di N righe e M colonne
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 64 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 64 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri interi o doubleing-point a precisione singola
*/
void save_data(char* filename, void* X, int n, int k) {
	FILE* fp;
	int i;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&k, 4, 1, fp);
		fwrite(&n, 4, 1, fp);
		for (i = 0; i < n; i++) {
			fwrite(X, sizeof(type), k, fp);
			//printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
			X += sizeof(type)*k;
		}
	}
	else{
		int x = 0;
		fwrite(&x, 4, 1, fp);
		fwrite(&x, 4, 1, fp);
	}
	fclose(fp);
}

// PROCEDURE ASSEMBLY

extern void prova(params* input);
extern type distEuclidea(VECTOR v1, VECTOR v2, int dim);
extern type pesoTot(VECTOR v, int dim);
extern type prodScalare(VECTOR v1, VECTOR v2,int dim);
extern void subVettori(VECTOR v1,VECTOR v2, VECTOR ris, int dim);
extern void addVettori(VECTOR v1,VECTOR v2, VECTOR ris, int dim);
extern void prodVet_x_Scalare(VECTOR v1, type s, VECTOR ris, int dim);
extern void prodVet_x_ScalareUn(VECTOR v1, type s, VECTOR ris, int dim);
///////////

VECTOR copyAlnVector(VECTOR v, int inizio, int dim){
    if(inizio % allineamento == 0)
        return v+inizio;
	VECTOR ret=get_block(sizeof(type),dim);
    int unrolling=16;
    int i=0;
	for(i=0;i+unrolling<dim;i+=unrolling){
		ret[i]=v[i+inizio];
        ret[i+1]=v[i+inizio+1];
        ret[i+2]=v[i+inizio+2];
        ret[i+3]=v[i+inizio+3];
        ret[i+4]=v[i+inizio+4];
        ret[i+5]=v[i+inizio+5];
        ret[i+6]=v[i+inizio+6];
        ret[i+7]=v[i+inizio+7];
	    ret[i+8]=v[i+inizio+8];
        ret[i+9]=v[i+inizio+9];
        ret[i+10]=v[i+inizio+10];
        ret[i+11]=v[i+inizio+11];
        ret[i+12]=v[i+inizio+12];
        ret[i+13]=v[i+inizio+13];
        ret[i+14]=v[i+inizio+14];
        ret[i+15]=v[i+inizio+15];
	}
    for(;i<dim;i++)
        ret[i]=v[i+inizio];
	return ret;	
}

/*
void addVettori(VECTOR v1,VECTOR v2, VECTOR ris, int dim){
	for(int i =0; i<dim; i++)
		ris[i]=v1[i]+v2[i];
}


type distEuclidea(VECTOR v1, VECTOR v2, int dim){
	type v=0;
	for(int i=0;i<dim;i++){
		v+= ((v2[i]-v1[i])*(v2[i]-v1[i]));
	}
	return (type)sqrt(v);
}

/*
type pesoTot(VECTOR v, int dim){
	type tmp=0;
	for(int i=0; i<dim;i++){
		tmp+=v[i];
	}
	return tmp;
}
*/
/*
void prodVet_x_Scalare(VECTOR v1, type s, VECTOR ris, int dim){	
	for(int i =0; i<dim; i++){
		ris[i]=v1[i]*s;
	}
}

void prodVet_x_ScalareUn(VECTOR v1, type s, VECTOR ris, int dim){	
	for(int i =0; i<dim; i++){
		ris[i]=v1[i]*s;
	}
}*/
/*
type prodScalare(VECTOR v1, VECTOR v2,int dim){
	type ris=0.0;
	for(int i =0; i<dim; i++){
		ris+=v1[i]*v2[i];
	}
	return ris;
}
*/
/*
void subVettori(VECTOR v1,VECTOR v2, VECTOR ris, int dim){
	for(int i =0; i<dim; i++){
		ris[i]=v1[i]-v2[i];
	}
}*/



type funzione(VECTOR vettore,params* input,int dim){
	type x2 = prodScalare(vettore,vettore,dim);
	type ex2 =(type)exp(x2);
	type cx = prodScalare(input->c,vettore,dim);
	
	return ex2+x2-cx;
}

type funzioneMatrix(MATRIX matrice,params* input,int inizio,int dim){
    //VECTOR vettore=matrice+inizio*dim*sizeof(type);
    VECTOR vettore=copyAlnVector(matrice,inizio,dim);
	type ret=funzione(vettore,input,dim);
	return ret;

}

type getRand(params* input,var* vars){
	return input->r[vars->rand++];
}

void addMatriceVettore(MATRIX matrix,VECTOR vector,int riga,int d){
    if((riga*d%allineamento)==0){
        addVettori(matrix+riga*d,vector,matrix+riga*d,d);
    }
    else{
        int unroll=4;
        int i=0;
	    for( i=0;i+unroll<d;i+=unroll){
		    matrix[riga*d+i]+=vector[i];
            matrix[riga*d+i+1]+=vector[i+1];
            matrix[riga*d+i+2]+=vector[i+2];
            matrix[riga*d+i+3]+=vector[i+3];
	    }
        for(;i<d;i++)
            matrix[riga*d+i]+=vector[i];
    }
}

void addMatriceVettoreBroad(MATRIX matrix,VECTOR vector,int d){
    int unroll=4;
    int pesce=0;    
    for(pesce=0;pesce+unroll<np;pesce+=unroll){
		addMatriceVettore(matrix,vector,pesce,d);
        addMatriceVettore(matrix,vector,pesce+1,d);
        addMatriceVettore(matrix,vector,pesce+2,d);
        addMatriceVettore(matrix,vector,pesce+3,d);
    }
    for(;pesce<np;pesce++)
        addMatriceVettore(matrix,vector,pesce,d);
    
}

void prodTrasMatVet(MATRIX matrix,VECTOR vector, VECTOR ris, int righe,int dim){	
    VECTOR parz=get_block(sizeof(type),d);    
    int unroll=4;    
    for(int i=0;i<dim;i++){
		ris[i]=0;
	}
    int pesce=0;
    if(allineamentoPerfetto){
         for(pesce=0;pesce+unroll<righe;pesce+=4){
            prodVet_x_Scalare(matrix+pesce*dim,vector[pesce],parz,dim);
            addVettori(ris,parz,ris,dim);
            prodVet_x_Scalare(matrix+(pesce+1)*dim,vector[pesce+1],parz,dim);
            addVettori(ris,parz,ris,dim);
            prodVet_x_Scalare(matrix+(pesce+2)*dim,vector[pesce+2],parz,dim);
            addVettori(ris,parz,ris,dim);
            prodVet_x_Scalare(matrix+(pesce+3)*dim,vector[pesce+3],parz,dim);
            addVettori(ris,parz,ris,dim);
        }
    
    }
    else{	
        for(pesce=0;pesce+unroll<righe;pesce+=4){
           
            if((pesce*dim%allineamento)==0){
                prodVet_x_Scalare(matrix+pesce*dim,vector[pesce],parz,dim);
                addVettori(ris,parz,ris,dim);
                prodVet_x_ScalareUn(matrix+(pesce+1)*dim,vector[pesce+1],parz,dim);
                addVettori(ris,parz,ris,dim);
                prodVet_x_ScalareUn(matrix+(pesce+2)*dim,vector[pesce+2],parz,dim);
                addVettori(ris,parz,ris,dim);
                prodVet_x_ScalareUn(matrix+(pesce+3)*dim,vector[pesce+3],parz,dim);
                addVettori(ris,parz,ris,dim);
            }
            else{
                prodVet_x_ScalareUn(matrix+pesce*dim,vector[pesce],parz,dim);
                addVettori(ris,parz,ris,dim);
                 prodVet_x_ScalareUn(matrix+(pesce+1)*dim,vector[pesce+1],parz,dim);
                addVettori(ris,parz,ris,dim);
                 prodVet_x_ScalareUn(matrix+(pesce+2)*dim,vector[pesce+2],parz,dim);
                addVettori(ris,parz,ris,dim);
                 prodVet_x_ScalareUn(matrix+(pesce+3)*dim,vector[pesce+3],parz,dim);
                addVettori(ris,parz,ris,dim);
            }
           
        }
    }
    for(;pesce<righe;pesce++){
        if((pesce*dim%allineamento)==0){
            prodVet_x_Scalare(matrix+pesce*dim,vector[pesce],parz,dim);
        }
        else{
            prodVet_x_ScalareUn(matrix+pesce*dim,vector[pesce],parz,dim);
        }
    }
        addVettori(ris,parz,ris,dim);
    free_block(parz);
               
    
}

type minimoVettore(VECTOR vector,int dim){
	int index=0;
	type min=vector[0];
	for(int i=0;i<dim;i++){
		if(vector[i]<min){
			index=i;
			min=(vector[i]);
		}
	}
	return vector[index];
}

void zeroRowMatrix(MATRIX matrix, int riga, int dim){
	for(int i=0;i<dim; i++){
		matrix[riga*dim+i]=0;
	}//for
}//zero

void replaceMatrixRowVector(MATRIX matrix,VECTOR vector,int riga, int dim){
	for(int i=0;i<dim; i++){
		matrix[riga*dim+i]=vector[i];
	}//for
}

void generaPosizioneCasuale(VECTOR pos,MATRIX posPesci,int pesce,int dim,params* input,var* vars ){
	for(int i=0;i<dim; i++){
		pos[i]=posPesci[pesce*dim+i]+(getRand(input,vars)*2-1)*input->stepind;
	}//for
}



///////////


void movimentoIndividuale(params* input,var* vars,int pesce){
    VECTOR newPosition=get_block(sizeof(type),d);
    generaPosizioneCasuale(newPosition,input->x,pesce,d,input,vars);
    type deltaf= funzione(newPosition,input,d)-funzioneMatrix(input->x ,input,pesce*d,d);  
    if(deltaf<0){
        effettuato=1;
        vars->deltaf[pesce]=deltaf;
        VECTOR x_i=copyAlnVector(input->x,pesce*d,d);
        //VECTOR x_i=input->x+pesce*d*sizeof(type);
        VECTOR deltax_i=get_block(sizeof(type),d);
        subVettori(newPosition,x_i,deltax_i,d);
        replaceMatrixRowVector(vars->deltax,deltax_i,pesce,d);
        replaceMatrixRowVector(input->x,newPosition,pesce,d);
        //if((pesce*d%allineamento)!=0)
            //free_block(x_i);
        //free_block(deltax_i);
    }//if
    else{
        zeroRowMatrix(vars->deltax,pesce,d);
        vars->deltaf[pesce]=0;
    }//else
    free_block(newPosition);
}//

void alimentazione(params* input, var* vars){
    vars-> wbranco=pesoTot(vars->w,np);
    //printf("%f ", vars->wbranco);
    if(effettuato){
        type max=-minimoVettore(vars->deltaf,np);
        if(max>EPSILON){
        	    VECTOR ris=get_block(sizeof(type),np);
            	prodVet_x_Scalare(vars->deltaf,(type)1.0/max,ris,np);
            	subVettori(vars->w,ris,vars->w, np);
            	free_block(ris);
        }//if
    }//if
}//alimentazione

void movimentoIstintivo(params* input, var* vars){
    if(effettuato){
        VECTOR I= get_block(sizeof(type),d);
        VECTOR num= get_block(sizeof(type),d);
        prodTrasMatVet(vars->deltax,vars->deltaf,num,np,d);
        type denom=(type)1.0/pesoTot(vars->deltaf,np);
        if(denom>EPSILON || denom<EPSILON){
            prodVet_x_Scalare(num,denom,I,d);
            addMatriceVettoreBroad(input->x,I,d);
        }//if
        free_block(I);
        free_block(num);	      
    }//if
}//mov istintivo

void baricentro(params* input, var* vars){
    type denom=(type)1.0/pesoTot(vars->w,np);
    prodTrasMatVet(input->x,vars->w,vars->baricentro,np,d);
    prodVet_x_Scalare(vars->baricentro,denom,vars->baricentro,d);
}//baricentro

void movimentoVolitivo(params* input, var* vars){ 
    int segno=1;
    if(pesoTot(vars->w,np)>vars->wbranco){
        segno=-1;
        //printf("pesoAumentato.");
    }
    VECTOR ris=get_block(sizeof(type),d);
    VECTOR volVec=get_block(sizeof(type),d);
    for(int pesce =0;pesce <np ;pesce++){
        type rnd=getRand(input,vars);
        VECTOR x_i=copyAlnVector(input->x,pesce*d,d);
        type distanza=distEuclidea(x_i,vars->baricentro,d);
        type scalare=(input->stepvol*rnd*(type)segno)/distanza;
        subVettori(x_i,vars->baricentro,ris,d);
        prodVet_x_Scalare(ris,scalare,volVec,d);
        addMatriceVettore(input->x,volVec,pesce,d);
        if((pesce*d%allineamento)!=0){
            free_block(x_i);
        }  
    }//for 
    free_block(ris);    
    free_block(volVec);      
}//movimentoVolitivo

void minimo(params* input){
	type valore_minimo = funzioneMatrix(input->x, input, 0, d); 
	int index = 0;
    
	for(int i=0; i<np; i++){
		type valore_tmp = funzioneMatrix(input->x, input, i*d, d); 
		//printf("valore_tmp:%f ",valore_tmp);
        if(valore_tmp<valore_minimo){
			valore_minimo=valore_tmp;
			index=i;
		}
	}	
	for(int i=0; i<d; i++){
		input->xh[i]=input->x[index*d+i];
	}	
}

void init(params* input, var* vars){
    d=input->d;
    np=input->np;
    if(d%allineamento==0)
        allineamentoPerfetto=1;
    iter=input->iter;
    vars->w=get_block(sizeof(type),np);
    vars->deltax=alloc_matrix(np,d);
    vars->deltaf=get_block(sizeof(type),np);
    vars->stepindIni=input->stepind;
    vars->stepvolIni=input->stepvol;
    vars->baricentro=get_block(sizeof(type),input->d);
    input->xh=get_block(sizeof(type),input->d);
    vars->rand=0;
    vars->wbranco=0;
    for(int i=0;i<input->np;i++){
        vars->w[i]=input->wscale/2;
        vars->wbranco+= vars->w[i];
    }
    for(int i=0;i<input->d;i++){
    	vars->baricentro[i]=0;
    	input->xh[i]=0;
    }
  
}

void aggiornaParametri(params* input, var* vars){
	input->stepind=input->stepind-(vars->stepindIni/iter);
	input->stepvol=input->stepvol-(vars->stepvolIni/iter);
}


void fss(params* input){
    int it =0;   
    var* vars=get_block(sizeof(var),1);
    init(input,vars);
    while (it<iter){
     	effettuato=0;
        for(int pesce=0;pesce<np;pesce++){
            movimentoIndividuale(input,vars,pesce);
        }
        
        alimentazione(input,vars);
        movimentoIstintivo(input,vars);
        baricentro(input,vars);
        movimentoVolitivo(input,vars);
        aggiornaParametri(input,vars);
    	it+=1;
    }
    minimo(input);	
}

int main(int argc, char** argv) {

	char fname[256];
	char* coefffilename = NULL;
	char* randfilename = NULL;
	char* xfilename = NULL;
	int i, j, k;
	clock_t t;
	double time;
	
	//
	// Imposta i valori di default dei parametri
	//

	params* input = malloc(sizeof(params));

	input->x = NULL;
	input->xh = NULL;
	input->c = NULL;
	input->r = NULL;
	input->np = 25;
	input->d = 2;
	input->iter = 350;
	input->stepind = 1;
	input->stepvol = 0.1;
	input->wscale = 10;
	
	input->silent = 0;
	input->display = 0;

	//
	// Visualizza la sintassi del passaggio dei parametri da riga comandi
	//

	if(argc <= 1){
		printf("%s -c <c> -r <r> -x <x> -np <np> -si <stepind> -sv <stepvol> -w <wscale> -it <itmax> [-s] [-d]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\tc: il nome del file ds2 contenente i coefficienti\n");
		printf("\tr: il nome del file ds2 contenente i numeri casuali\n");
		printf("\tx: il nome del file ds2 contenente le posizioni iniziali dei pesci\n");
		printf("\tnp: il numero di pesci, default 25\n");
		printf("\tstepind: valore iniziale del parametro per il movimento individuale, default 1\n");
		printf("\tstepvol: valore iniziale del parametro per il movimento volitivo, default 0.1\n");
		printf("\twscale: valore iniziale del peso, default 10\n");
		printf("\titmax: numero di iterazioni, default 350\n");
		printf("\nOptions:\n");
		printf("\t-s: modo silenzioso, nessuna stampa, default 0 - false\n");
		printf("\t-d: stampa a video i risultati, default 0 - false\n");
		exit(0);
	}

	//
	// Legge i valori dei parametri da riga comandi
	//

	int par = 1;
	while (par < argc) {
		if (strcmp(argv[par],"-s") == 0) {
			input->silent = 1;
			par++;
		} else if (strcmp(argv[par],"-d") == 0) {
			input->display = 1;
			par++;
		} else if (strcmp(argv[par],"-c") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing coefficient file name!\n");
				exit(1);
			}
			coefffilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-r") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing random numbers file name!\n");
				exit(1);
			}
			randfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-x") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing initial fish position file name!\n");
				exit(1);
			}
			xfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-np") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing np value!\n");
				exit(1);
			}
			input->np = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-si") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing stepind value!\n");
				exit(1);
			}
			input->stepind = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-sv") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing stepvol value!\n");
				exit(1);
			}
			input->stepvol = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-w") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing wscale value!\n");
				exit(1);
			}
			input->wscale = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-it") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing iter value!\n");
				exit(1);
			}
			input->iter = atoi(argv[par]);
			par++;
		} else{
			printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
			par++;
		}
	}

	//
	// Legge i dati e verifica la correttezza dei parametri
	//

	if(coefffilename == NULL || strlen(coefffilename) == 0){
		printf("Missing coefficient file name!\n");
		exit(1);
	}

	if(randfilename == NULL || strlen(randfilename) == 0){
		printf("Missing random numbers file name!\n");
		exit(1);
	}

	if(xfilename == NULL || strlen(xfilename) == 0){
		printf("Missing initial fish position file name!\n");
		exit(1);
	}

	int x,y;
	input->c = load_data(coefffilename, &input->d, &y);
	input->r = load_data(randfilename, &x, &y);
	input->x = load_data(xfilename, &x, &y);

	if(input->np < 0){
		printf("Invalid value of np parameter!\n");
		exit(1);
	}

	if(input->stepind < 0){
		printf("Invalid value of si parameter!\n");
		exit(1);
	}

	if(input->stepvol < 0){
		printf("Invalid value of sv parameter!\n");
		exit(1);
	}

	if(input->wscale < 0){
		printf("Invalid value of w parameter!\n");
		exit(1);
	}

	if(input->iter < 0){
		printf("Invalid value of it parameter!\n");
		exit(1);
	}

	//
	// Visualizza il valore dei parametri
	//

	if(!input->silent){
		printf("Coefficient file name: '%s'\n", coefffilename);
		printf("Random numbers file name: '%s'\n", randfilename);
		printf("Initial fish position file name: '%s'\n", xfilename);
		printf("Dimensions: %d\n", input->d);
		printf("Number of fishes [np]: %d\n", input->np);
		printf("Individual step [si]: %f\n", input->stepind);
		printf("Volitive step [sv]: %f\n", input->stepvol);
		printf("Weight scale [w]: %f\n", input->wscale);
		printf("Number of iterations [it]: %d\n", input->iter);
	}

	// COMMENTARE QUESTA RIGA!
	//prova(input);
	//

	//
	// Fish School Search
	//

	t = clock();
	fss(input);
	t = clock() - t;
	time = ((double)t)/CLOCKS_PER_SEC;

	if(!input->silent)
		printf("FSS time = %.3f secs\n", time);
	else
		printf("%.3f\n", time);

	//
	// Salva il risultato di xh
	//
	sprintf(fname, "xh64_%d_%d_%d.ds2", input->d, input->np, input->iter);
	save_data(fname, input->xh, 1, input->d);
	if(input->display){
		if(input->xh == NULL)
			printf("xh: NULL\n");
		else{
			printf("xh: [");
			for(i=0; i<input->d-1; i++)
				printf("%f,", input->xh[i]);
			printf("%f]\n", input->xh[i]);
		}
	}

	if(!input->silent)
		printf("\nDone.\n");

	return 0;
}
