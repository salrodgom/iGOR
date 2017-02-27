#ifndef _NODE_STRUCT_H
#define _NODE_STRUCT_H

 struct node { //Esto genera dos vectores de punteros
  int kin;
  int kout;
  int label;
  node ** in_nodes;
  node ** out_nodes;
} ;

 struct nodo { //Esto genera dos vectores de enteros
  int kin;
  int kout;
  int label;
  int myself;
  double s;  
  int * in_nodos;
  double * in_sig;
  int * out_nodos;
  double * out_sig;
} ;

 struct edge {

   int in;
   int out;
   int importance;
   int *loop;
   int label;
   int label2;

};


#endif
