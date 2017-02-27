#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <new>
#include <vector>
#include <string.h>
#include <stddef.h>
#include <iostream>
#include <omp.h>
#include "node_struct.h"
using namespace std;

//estas variables hace falta porque se meten como externas en la hoja de funciones
int *L,*Lp,*Ln;

//DEFINICIONES DE FUNCIONES
nodo * leer_red_new(int *N,char *nombre); //lee la red
nodo * simetrizar(nodo *M,int N); //simetriza la red pa porsi
void loops(int NODO,nodo *M,int cont,int signo,int final,int * rastro,int * estado,int contmax); //saca los loops

//---------------------------------------------------------------------------//


/////////////////////////////////////////////////////////////////////////////
//                  PROGRAMA PRINSIPAL  ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]){

//propiedades de las redes--------------
int N;//numero de nodos de la red
char nombre[30];//almacena el nombre del archivo de red
nodo *M; //la red
nodo *Msim; //la red simetrizada
//----------------------------------
//Calculo paralelo -----------------	
int iCPU = omp_get_num_procs();
//printf("#el n procesadores es %d\n",iCPU);
iCPU=1;
int num_proc = iCPU;
//----------------------------------
//cosas de loops--------------------
int contmax,cont_ciclos; //longitud maxima de los loops
contmax=4; //la longitud maxima de los looops q vas buscando, si no hay pon N.
L=new int [N+1]; Lp=new int [N+1]; Ln=new int [N+1];
//---------------------------------

//------------ ENTRADA DE DATOS --------------------------------------------//
if(argc < 2 ){
	        printf("ERROR: debes dar el nombre y la longitud de los loops\n");
		        return 1;
}
	
	strcpy(nombre,argv[1]);//copiamos el nombre de la red, necesario apra leer del archivo. Tiene q tener extension .red 
	contmax=atoi(argv[2]);
	//Esto lee el archivo y almacena la red en M
	M=leer_red_new(&N,nombre);
	//Esto simetriza la red por si la que leemos originalmente no viene simetrizada
	Msim=simetrizar(M,N);
	
//--------------Buscamos los loops-----------------------------------------//
	for (int i=0;i<N+1;i++) L[i]=0;//Ponemos los contadores de loops a cero
				#pragma omp parallel for schedule(dynamic,1) num_threads(num_proc) 
     				for (int i=0;i<N;i++){ 
       				//printf("Desde %d\n",i);
				//int this_thread = omp_get_thread_num(), num_threads = omp_get_num_threads();
				////printf("Hola, soy la hebra %d de un total de %d\n",this_thread,num_threads);
				int cont=0;
				int signo=1;//el signo de un loop particular
				int *rastro;
				int *estado;
				estado=new int[N];
				for(int j=0;j<N;j++) estado[j]=1;
	
				////printf("encontramos loops desde %d \n",i);
				rastro=(int *)malloc(sizeof(int)*1000);	
	  			loops(i,Msim,cont,signo,i,rastro,estado,contmax);
	
				delete []estado;
				free(rastro);

    				 }//End quitar nodos --------------------------------------fin cilco principal---- a partir de aqui no hay mas M!! ahora si,pqno la borro
				//return 0;
				//return 0;
				L[2]=L[2]*2; //estos solo aparecen en una direccion por culpa del algoritmo de busqueda, asi q multiplicamos por dos
				for (int i=1;i<contmax+1;i++) printf("%d %d \n",i,L[i]);	
	//--------------------------------------------------------------------------------------------------------
	
	return 0;
}

//funciones
//--------------------LEER RED-----------------------------------------------//
nodo * leer_red_new(int *N,char *nombre){
int nodoi,nodoj,indice,cosa,cont,ord_c;
float qq;
double xvar1;
char signo[5];
int signoint;
char nombre_arch[20];
char basura[4];
vector <vector <int> > aux_in; //auxiliar de entrada para cada nodo
vector <vector <int> > sig_in;
vector <vector <int> > aux_out; //auxiliar de salida para cada nodo
vector <vector <int> > sig_out;
int n;
FILE *fichero;

*N=0;// aun no sabemos cuantos nodos tiene
sprintf(nombre_arch,"%s.red",nombre); //generamos el nombre del archivo a usar
fichero=fopen(nombre_arch,"r");
	while(!feof(fichero)) { //sacamos el numero de nodos
	fscanf(fichero,"%s\n",basura);
	if( (strcmp(basura,"+")) or (strcmp(basura,"-"))) {
	   fscanf(fichero,"%lf\n", &xvar1);
	   nodoi=(int)xvar1;
	   //printf("leido %lf %d\n",xvar1,nodoi);
	   if (nodoi>*N){*N=nodoi;}
	   }
	}
	aux_in.resize(*N+1);//creo espacio para guardar las listas
	sig_in.resize(*N+1);
	aux_out.resize(*N+1);//creo espacio para guardar las listas
	sig_out.resize(*N+1);
	rewind(fichero);
	//printf("#ya tenemos N=%d\n",*N);
	while(!feof(fichero)) { // creamos la red
		fscanf(fichero,"%lf\n", &xvar1);
        	nodoi=(int)xvar1-1;
		//printf("nodoi=%i   ",nodoi);
       		fscanf(fichero,"%lf\n", &xvar1);
      		nodoj=(int)xvar1-1;
      		//printf("nodoj=%i   ",nodoj);
      		fscanf(fichero,"%s\n",signo);
      		//printf("signo %s\n",signo);
      		if (strcmp(signo,"-")== 0) signoint=-1;
      		else if (strcmp(signo,"+")== 0) signoint=1;
      		else if (strcmp(signo,"0")== 0) signoint=-1;//esto tengo q mirarlo!!!
      		else (signoint=8);
      		//printf("signo=%d\n",signoint);
      		//printf("itnento meter en auxin[%d] auxout[%d] de %d:\n",nodoi,nodoj,*N+1);     		
      		aux_in[nodoi].push_back(nodoj);
      		sig_in[nodoi].push_back(signoint);
      		aux_out[nodoj].push_back(nodoi);
      		sig_out[nodoj].push_back(signoint);
      		
      		signoint=0;
      		//printf("el nodo %d se come a ",nodoi);
      		//for (int i=0;i<aux_in[nodoi].size();i++) printf("%d ",aux_in[nodoi][i]);
      		//printf("\n");
      		
   	}
//printf("va bien la asignacion? \n");
//return;
/*for(int i=0;i<*N;i++){*/
/*printf("M[%d]in: ",i);*/
/*	for (int j=0;j<aux_in[i].size();j++){*/
/*		//printf("%d*",sig_in[i][j]);*/
/*		printf("%d ",aux_in[i][j]);*/
/*	}*/
/*	printf("\n");*/
/*printf("M[%d]out: ",i);*/
/*	for (int j=0;j<aux_out[i].size();j++){*/
/*		//printf("%d*",sig_out[i][j]);*/
/*		printf("%d ",aux_out[i][j]);*/
/*	}*/
/*	printf("\n");*/
/*}*/
n=*N;
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//printf("#N=%d n=%d\n",*N,n);
	nodo *my_node;
	my_node = new nodo [*N];
	for (int i=0;i<*N;i++) { my_node[i].kout=0; my_node[i].label=1; my_node[i].myself=i;}//pongo todas las conectividades out a cero
//printf("lo metemos dentro de la final\n");
for (int i=0;i<*N;i++){
	my_node[i].kin=aux_in[i].size();
	my_node[i].in_nodos=new int [my_node[i].kin];
	my_node[i].in_sig=new double [my_node[i].kin];
//printf("Kin[%d]=%d %d \n ",i,aux_in[i].size(),my_node[i].kin);
//printf("M[%d]in=",i);
	for (int j=0;j<my_node[i].kin;j++){
		my_node[i].in_nodos[j]=aux_in[i][j];
		my_node[i].in_sig[j]=sig_in[i][j];
//		printf("%d*(%.0f) ",my_node[i].in_nodos[j],my_node[i].in_sig[j]);
	}
//	printf("\n");
	
	my_node[i].kout=aux_out[i].size();
//printf("kout=%d %d \n",aux_out[i].size(),my_node[i].kout);
//printf("M[%d]out=",i);
	my_node[i].out_nodos=new int [my_node[i].kout];
	my_node[i].out_sig=new double [my_node[i].kout];
	for (int j=0;j<my_node[i].kout;j++){
		my_node[i].out_nodos[j]=aux_out[i][j];
		my_node[i].out_sig[j]=sig_out[i][j];
//		printf("%d*(%.0f) ",my_node[i].out_nodos[j],my_node[i].out_sig[j]);
	}
//	printf("\n");
}
//printf("hemos copiado --------------------------------------------  \n");
/*for (int i=0;i<*N;i++){*/
/*printf("Nodo %d , Kin=%d se come a :",i,my_node[i].kin);*/
/*	for (int j=0;j<my_node[i].kin;j++){*/
/*		printf("(%.0f)*%d ",my_node[i].in_sig[j],my_node[i].in_nodos[j]+1);*/
/*		*/
/*	}*/
/*	printf("\n");*/
/*printf("es comido por %d bichos:",my_node[i].kout);*/
/*for (int j=0;j<my_node[i].kout;j++){*/
/*		printf("(%.0f)*%d ",my_node[i].out_sig[j],my_node[i].out_nodos[j]+1);*/
/*		*/
/*	}*/
/*	printf("\n");*/
/*}*/
/*printf("terminamos?\n");*/

	aux_in.clear();
	aux_out.clear();
	
	//for (int i=0;i<N;i++) my_node[i].label=i;
	
	return my_node;
}//------------------------------------------------------------------------//
//----------------------------------------------------------------------------
nodo * simetrizar(nodo *M,int N)
{	//concretamente pone en una direccion todas las conexiones q habia antes
	vector < vector <int> > aux_in;
	vector < vector <int> > aux_out;
	aux_in.resize(N);
	aux_out.resize(N);
	
	int repe;
	printf("entramos a simemtrizar N=%d 1\n",N);
	

printf("simetrizamos\n");

	//printf("copiamos los vectores a los auxiliares\n");
	for (int i=0;i<N;i++){
		for (int j=0;j<M[i].kin;j++){
			aux_in[i].push_back(M[i].in_nodos[j]);	
			//printf("M[%d].in_modos[%d]=%d aux_in[%d][%d]=%d\n",i,j,M[i].in_nodos[j],i,j,aux_in[i][j]);	
		}
		for (int j=0;j<M[i].kout;j++){
			aux_out[i].push_back(M[i].out_nodos[j]);
			//printf("M[%d].out_modos[%d]=%d aux_out[%d][%d]=%d\n",i,j,M[i].out_nodos[j],i,j,aux_out[i][j]);
		}
	}

	printf("ponemos los q faltan\n");
	int nodo2;
	for (int i=0;i<N;i++){
		for (int j=0;j<M[i].kin;j++){
			//tenemos que ver si j se come a i
			repe=0;
			nodo2=M[i].in_nodos[j];
			for(int k=0;k<aux_in[nodo2].size();k++){
				if(i==aux_in[nodo2][k]) repe=1;
			}
			if (repe==0) aux_in[nodo2].push_back(i); //El q tiene toda la informacion es aux_in
		}
		for (int j=0;j<M[i].kout;j++){
			//tenemos que ver si por el q es comido es comido por i ya
			repe=0;
			nodo2=M[i].out_nodos[j];
			for (int k=0;k<aux_out[nodo2].size();k++){
				if (i==aux_out[nodo2][k]) repe=1;
			}
			if (repe==0) aux_out[nodo2].push_back(i);
		}
	}
	printf("voy a borrar\n");

	nodo *my_node;
	my_node = new nodo [N];
	//return my_node;
	
	
	for (int i=0;i<N;i++){
		my_node[i].kin=aux_in[i].size();
		my_node[i].in_nodos=new int [my_node[i].kin];
		my_node[i].in_sig=new double [my_node[i].kin];
		for (int j=0;j<my_node[i].kin;j++){
			my_node[i].in_nodos[j]=aux_in[i][j];
		}
		my_node[i].kout=aux_out[i].size();
		my_node[i].out_nodos=new int [my_node[i].kout];
		my_node[i].out_sig=new double [my_node[i].kout];
		for (int j=0;j<my_node[i].kout;j++){
			my_node[i].out_nodos[j]=aux_out[i][j];
		}
	}
//	printf("matriz simetrica\n");
/*	for (int i=0;i<N;i++){*/
/*	printf("Nodo %d , se come a %d bischos :",i,my_node[i].kin);*/
/*		for (int j=0;j<my_node[i].kin;j++){*/
/*			printf("%d ",my_node[i].in_nodos[j]);*/
/*		*/
/*		}*/
/*		printf("\n");*/
/*	printf("%d es comido por %d bichos:",i,my_node[i].kout);*/
/*		for (int j=0;j<my_node[i].kout;j++){*/
/*			printf("%d ",my_node[i].out_nodos[j]);*/
/*		*/
/*		}*/
/*		printf("\n");*/
/*	}*/
/*	printf("liberamos memoria\n");*/
// liberamos memoria ---------------------------------------------
	for (int i=0;i<N;i++) {
		if (aux_in[i].size() > 0) aux_in[i].clear();
		if (aux_out[i].size() > 0) aux_out[i].clear();
	}
	aux_in.clear();
	aux_out.clear();
	return my_node;
}
//--------------------------------------------------------------------------//
void loops(int NODO,nodo *M,int cont,int signo,int final,int * rastro,int * estado,int contmax){//saca los loops
int nodo2;
//printf("nivel %d : %d   ",cont,NODO);
//M[NODO].label=0;
estado[NODO]=0;
rastro[cont]=NODO;
cont=cont+1;

for (int i=0;i<M[NODO].kout;i++){
	nodo2=M[NODO].out_nodos[i];
	//printf("%d ",nodo2);
	if (nodo2 == final) {
		//printf("hemos encontrado uno de tamaño %d: \n",cont);
		#pragma omp atomic
		L[cont]++;
		//ciclos[cont_ciclos]=new int [cont];
		//DSECOMENTAR SI QUIERES QUE SAQUE LOS LOOPS 1 a 1
		for(int o=0; o<cont; o++){
			//ciclos[cont_ciclos][o]=rastro[o];
			printf("%d ",rastro[o]+1);
		}
		printf("%d",final+1);
		printf("\n");
		//sleep(1);
		//cont_ciclos ++;
	}
	else {
		if ((estado[nodo2] == 1) and (nodo2 > final) and (cont < contmax)) {//lo visitamos ;(nodo2 > final) ->eswto es como si eliminasemos los nodos de los q sabemos todos los loops
			//M[nodo2].label=0;
			//printf("level %d :",cont);
			signo=signo*M[NODO].out_sig[i];
			loops(nodo2,M,cont,signo,final,rastro,estado,contmax);
			//printf("%d*%d    ",M[NODO].out_sig[i],nodo2);
		}
		else {	

		}
	}

}
//M[NODO].label=1;
estado[NODO]=1;

return;
}

