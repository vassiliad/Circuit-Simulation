#include <stdio.h>
#include <stdlib.h>

int *nodes=NULL;
FILE **files = NULL;
static int num_nodes = 0;
extern int unique_hash;

void plot_node(int id, char *str) 
{
  char temp[1230];

  nodes = ( int* ) realloc(nodes, sizeof(int)*(num_nodes+1));
  files = ( FILE**) realloc(files, sizeof(FILE*) * (num_nodes+1));

  sprintf(temp, "plot_v_%s", str);
  files[num_nodes] = fopen(temp, "w");
  nodes[num_nodes] = id;
  num_nodes++;
}

void print_plots(double x, double *sol, int *P)
{
  int i =0;
	if ( P )
		for (i=0; i<num_nodes; i++ )
			fprintf(files[i], "%10g %10g\n", x ,sol[P[i]]);
	else
		for (i=0; i<num_nodes; i++ )
			fprintf(files[i], "%10g %10g\n", x ,sol[i]);

}

void plot_finalize() {
  int i;
  for (i=0;i<num_nodes; i++)
    fclose(files[i]);
	
	free(files);
	free(nodes);
}

