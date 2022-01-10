/**********
* 
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2020/21
* 
* Progetto dell'algoritmo Fish School Search 221 231 a
* in linguaggio assembly x86-32 + SSE
* 
* Fabrizio Angiulli, aprile 2019
* 
**********/

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
* nasm -f elf32 fss32.nasm && gcc -m32 -msse -O0 -no-pie sseutils32.o fss32.o fss32c.c -o fss32c -lm && ./fss32c $pars
* 
* oppure
* 
* ./runfss32
* 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

#define	type		float
#define	MATRIX		type*
#define	VECTOR		type*
#define RANDMAX 85050
#define EPSILON 0.0000001

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
} var;

int np;
int d;
int iter;

/*
* 
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate 
* 	mediante un array (float*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere 
* 	memorizzate mediante array di array (float**).
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
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri floating-point a precisione singola
* 
*********
*	Se lo si ritiene opportuno, è possibile cambiare la codifica in memoria
* 	della matrice. 
*********
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
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri interi o floating-point a precisione singola
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
extern void addVettori(VECTOR v1,VECTOR v2, VECTOR ris, int dim);
extern void subVettori(VECTOR v1,VECTOR v2, VECTOR ris, int dim);
extern void prodVet_x_Scalare(VECTOR v1, type s, VECTOR ris,int dim);
extern type prodScalare(VECTOR v1, VECTOR v2,int dim);
extern type pesoTot(VECTOR v, int dim);
extern VECTOR copyAlnVector(VECTOR v, int inizio, int dim);
extern type distEuclidea(VECTOR v1, VECTOR v2, int dim);
/*
type prodScalare(VECTOR v1, VECTOR v2,int inizio1,int inizio2,int dim){
	type ris=0.0;
	for(int i =0; i<dim; i++){
		ris+=v1[inizio1+i]*v2[inizio2+i];
	}
	return ris;
}
*/
/*
void subVettori(VECTOR v1,VECTOR v2, VECTOR ris,int inizio1, int inizio2, int dim){
	for(int i =0; i<dim; i++){
		ris[i]=v1[i+inizio1]-v2[i+inizio2];
	}
	//return ris;
}*/

/*void prodVet_x_Scalare(VECTOR v1, type s, VECTOR ris, int inizio,int dim){	
	for(int i =0; i<dim; i++){
		ris[i]=v1[i+inizio]*s;
		//printf(" v[i+inizio] %f, x %f = %f \n \n \n",v1[i],s,ris[i]);
	}
}*/
/*
VECTOR copyAlnVector(VECTOR v, int inizio, int dim){
	VECTOR ret=get_block(sizeof(type),dim);
	for(int i=0;i<dim;i++){
		ret[i]=v[i+inizio];
	}
	return ret;	
}
*/


type funzione(VECTOR vettore,params* input,int inizio,int dim){
	VECTOR v=copyAlnVector(vettore,inizio,dim);
	type x2 = prodScalare(v,v,dim);
	type ex2 =(type)exp(x2);
	type cx = prodScalare(input->c,v,dim);
	//free_block(v);
	//printf("%f,%f, %f \n", ex2, x2, cx);
	return ex2+x2-cx;
}


/*
type distEuclidea(VECTOR v1, VECTOR v2, int inizio1, int inizio2, int dim){
	type v=0;
	for(int i=0;i<dim;i++){
		v+= ((v2[i+inizio2]-v1[i+inizio1])*(v2[i+inizio2]-v1[i+inizio1]));
	}
	return (type)sqrtf(v);
}
*/
/*
type pesoTot(VECTOR v, int dim){
	type tmp=0;
	for(int i=0; i<dim;i++){
		tmp+=v[i];
	}
	return tmp;
}
*/

void minimo(params* input){	
	type valore_minimo = funzione(input->x, input, 0, input->d); 
	int index = 0;
	for(int i=0; i<np; i++){
		type valore_tmp = funzione(input->x, input, i*d, d);
		//printf("%f \n",valore_tmp); 
		if(valore_tmp<valore_minimo){
			valore_minimo=valore_tmp;
			index=i;
		}
	}	
	for(int i=0; i<d; i++){
		input->xh[i]=input->x[index*d+i];
	}	
}


void movimentoIndividuale(params* input,var* vars,int pesce){
    VECTOR newPosition= get_block(sizeof(type),d);
    for(int i=0;i<d;i++){
    	newPosition[i]=0;
    }
    for(int i=0;i<d;i++){
        newPosition[i]=input->x[pesce*d+i] + ((2* (input->r[vars->rand]))-1)*input->stepind;
        vars->rand=(vars->rand+1)%RANDMAX;
       // printf("pesce %d, vecchio valore di x%d, %f, nuovo %f \n",pesce,i,input->x[pesce*input->d+i],newPosition[i]);     
    }
    type deltaf=(funzione(newPosition,input,(int)0,d))-funzione((VECTOR)input->x,input,pesce*d,d);

    if(deltaf<0){
    	//printf("effettuato movimento Individuale \n");
        VECTOR ris=get_block(sizeof(type),d);
	VECTOR v2 =copyAlnVector((VECTOR)input->x,pesce*d,d);
        subVettori(newPosition,v2,ris,d);

        for(int i=0;i<d;i++){
    		input->x[pesce*d+i]=newPosition[i]; //ciclo di copia  
    	}  
        vars->deltaf[pesce]=deltaf; //assegnamento
       	//printf("deltaf %f \n",deltaf);
        for(int i=0;i<d;i++){        
        	vars->deltax[pesce*d+i]= ris[i];
        	
        }
        free_block(ris);
        //free_block(v2);
    }else{
    	vars->deltaf[pesce]=0;
        for(int i=0; i< d;i++){
    	    vars->deltax[pesce*d+i]=0;
    	}
    }
    free_block(newPosition);
}	

void alimentazione(params* input, var* vars){
    vars->wbranco=0;
    type max=fabs((double)vars->deltaf[0]);
    for(int pesce=0;pesce<np;pesce+=1){
    	vars->wbranco+=vars->w[pesce];
        //printf("deltaf x max %f \n",vars->deltaf[pesce]);
        if(fabs((double)vars->deltaf[pesce])>max)
            max=fabs((double)vars->deltaf[pesce]);
    }
    if(!(max<EPSILON&&max>-EPSILON)){
   	//printf("max deltaf %f \n",max);
   	for(int pesce=0;pesce<np;pesce+=1){
   		vars->w[pesce]=vars->w[pesce]- vars->deltaf[pesce]/max; //- al posto di +
   	}//for
    }//if
}//alimentazione


void movimentoIstintivo(params* input, var* vars){
//CALCOLO LA SOMMATORIA da 1 a n di deltax(vettore) per deltaf(scalare) 
		VECTOR num = get_block(sizeof(type),input->d);
		VECTOR I= get_block(sizeof(type),input->d);
		//ciclo di inizializzazione di I
		for(int i=0; i<d; i++){
			I[i]=0;
			num[i]=0;
		}//ciclo azzeramento
		for(int i=0; i<input->np; i++){
			VECTOR v=get_block(sizeof(type),d);
			//copyAlnVector(vars->deltax,v,i*d,d);
			//prodVet_x_Scalare(v,vars->deltaf[i],num, d);//il risultato va in contatore
			for(int k=0; k<d; k++){
				//printf("%f ",vars->deltax[k]);
				}
			//addVettori(VECTOR v1,VECTOR v2, VECTOR ris,int inizio1, int inizio2, int dim)  <----
			addVettori(I,num,I,d);
			free_block(v);
			//printf("I: %f,%f,%f,%f,%f,%f,%f \n", I[0],I[1],I[2],I[3],I[4],I[5],I[6]);
		}
	
		type sommaDeltaf=0.0;
		for(int i=0; i<np; i++){
			if(vars->deltaf[i]!=0)
			//printf("deltaf considerati :%f \n", vars->deltaf[i]);
			sommaDeltaf+=vars->deltaf[i];
		}
		//printf("somma deltaf %f \n", sommaDeltaf);
		if(sommaDeltaf>EPSILON || sommaDeltaf<-EPSILON){
			//printf("movimento istintivo \n");
			prodVet_x_Scalare(I, (type)1/sommaDeltaf, I, d);
			//printf("I: %f,%f,%f,%f,%f,%f,%f \n", I[0],I[1],I[2],I[3],I[4],I[5],I[6]);
//per ogni pesce
			for(int pesce=0; pesce<input->np; pesce++){
				VECTOR ris=get_block(sizeof(type),d);
				VECTOR v1=copyAlnVector(input->x,pesce*d,d);
				addVettori(v1,I,ris,d);
				//printf("Nuove coordinate pesce %d dopo mov Instintivo: ",pesce);
				for(int i=0;i < d; i++){
					input->x[pesce*d+i]=ris[i];
					//printf("%f " ,ris[i] );
				}
				//free_block(v1);
			}
		}//if
	free_block(num);
	free_block(I); 
	
}//mov istintivo

void baricentro(params* input, var* vars){
	//int dim = input->np*d;
	
	VECTOR prod_tmp= get_block(sizeof(type),input->d);
	VECTOR num= get_block(sizeof(type),input->d);
	for(int i=0; i<input->d;i++){
		num[i]=0;
		prod_tmp[i]=0;
	}
	
	type denominatore=0;
	//mi muovo a blocchi di d su x 
	for(int i=0; i<np; i++){
		//prodVet_x_Scalare(VECTOR v1, int s, VECTOR ris, int inizio,int fine)
		VECTOR v=copyAlnVector(input->x,i*d,d);
		prodVet_x_Scalare(v,vars->w[i],prod_tmp,d);//<-- Problema su x
		for(int i=0; i<input->d;i++){
		
		//printf("%f ",prod_tmp[i]);
		}
		//printf("\n");
		addVettori(num, prod_tmp, num, d);//num(vettore)=valore corrente+ prod_tmp
	}

	for(int i =0; i<input->np; i++){
		denominatore += vars->w[i];
	}
	prodVet_x_Scalare(num, (type)1/denominatore, vars->baricentro, d);
	//printf("Baricentro: %f,%f,%f,%f,%f,%f,%f,%f \n",vars->baricentro[0],vars->baricentro[1],vars->baricentro[2],vars->baricentro[3],vars->baricentro[4],vars->baricentro[5],vars->baricentro[6],vars->baricentro[7]);
	free_block(prod_tmp);
	free_block(num);
}

void movimentoVolitivo(params* input, var* vars){
	VECTOR diff = get_block(sizeof(type),input->d);
	type pesoTotAtt = pesoTot(vars->w,input->np);
	type pesoTotPre = vars->wbranco;
	//printf(" wAtt;%f,  wPre%f \n",pesoTotAtt,pesoTotPre );
	if(pesoTotAtt>pesoTotPre){ //segno "-" nell'equazione (il banco si avvicina al baricentro)
		for(int pesce=0;pesce<np; pesce++){
			VECTOR v1=copyAlnVector(input->x,pesce*d,d);
			subVettori(v1,vars->baricentro,diff,d);
			//type distEuclidea(VECTOR v1, VECTOR v2, int inizio1, int inizio2, int dim){
			type dist = distEuclidea(v1, vars->baricentro,d);
			for(int i=0;i<d;i++){
				if(dist>EPSILON || dist<-EPSILON){
				input->x[pesce*d+i]=input->x[pesce*d+i]-input->stepvol* input->r[vars->rand] *(diff[i]/dist);
				//void prodVet_x_Scalare(VECTOR v1, float s, VECTOR ris, int inizio,int dim){	
				vars->rand++;}
			}
	        }
	}
	else{ //segno "+" nell'equazione (il banco si allontana dal baricentro)
		for(int pesce=0;pesce<np; pesce++){
			VECTOR v1=copyAlnVector(input->x,pesce*d,d );
			subVettori(v1,vars->baricentro,diff, d);
			type dist = distEuclidea(v1, vars->baricentro,d);
			for(int i=0;i<d;i++){
				if(dist>EPSILON || dist<-EPSILON){
				input->x[pesce*d+i]=input->x[pesce*d+i]+input->stepvol* input->r[vars->rand] *(diff[i]/dist);
				vars->rand++;
				}
			}
	       		//free_block(v1);
	        }
	}
	free_block(diff);
	
}

void aggiornaParametri(params* input, var* vars){
	input->stepind=input->stepind-(vars->stepindIni/iter);
	input->stepvol=input->stepvol-(vars->stepvolIni/iter);
}

void init(params* input, var* vars){
    vars->w=get_block(sizeof(type),input->np);
    vars->deltax=alloc_matrix(input->np,input->d);
    vars->deltaf=get_block(sizeof(type),input->np);
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
    d=input->d;
    np=input->np;
    iter=input->iter;
}


void fss(params* input){
    int it =0;   
    var* vars=get_block(sizeof(var),1);
    init(input,vars);
    while (it<iter){
        
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
	float time;
	
	//
	// Imposta i valori di default dei parametri
	//

	params* input = get_block(sizeof(params),1);

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
		printf("\t-d: stgot ampa a video i risultati, default 0 - false\n");
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
	/*VECTOR v1=get_block(sizeof(type),8);
	VECTOR v2=get_block(sizeof(type),8);
	for(int i=0;i<8;i++){
		v1[i]=2.0;
		v2[i]=1.0*i ;
		//v2[7]=0;
		printf("%f \n",v1[i]*v2[i]);
		}
		
	printf("%f \n",prodScalare(v1,v2,0,0,7));
  */

	t = clock();
	fss(input);
	t = clock() - t;
	time = ((float)t)/CLOCKS_PER_SEC;

	if(!input->silent)
		printf("FSS time = %.3f secs\n", time);
	else
		printf("%.3f\n", time);

	//
	// Salva il risultato di xh
	//
	sprintf(fname, "xh32_%d_%d_%d.ds2", input->d, input->np, input->iter);
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
