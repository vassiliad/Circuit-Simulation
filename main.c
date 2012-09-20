#include <stdio.h>
#include "options.h"
#include "components.h"
#include "hash_table.h"
#include "mna.h"
#include "dc_instruction.h"
#include "transient.h"

extern FILE* yyin;
int yyparse();
int yylex_destroy();

enum IterativeMethods method_iter    = BiCG;
enum SolutionMethods method_choice  = NonIterative;
enum NonIterativeMethods method_noniter = LUDecomp;
enum TransientMethods method_tran = Tr;
int sparse_use = 0;

double itol = 1e-6;
extern int mna_size;

int main(int argc, char* argv[])
{
	int debug;
  yyin = fopen(argv[1], "r");
	
  if ( argc !=4 ) {
    printf("usage:\n%s <file.spice> <output_file> DEBUG\n", argv[0]);
    return 0;
  }
  
  f = fopen(argv[2],"w");
	name_of_file = argv[2];
	
	debug = atoi(argv[3]);
	
		
  
  if ( yyin == NULL ) {
    printf("[-] Could not open %s\n", argv[1]);
    return 0;
  }

  hash_initialize();

  if ( yyparse() == 0 ) {
    printf("[+] Parsed %s succesfully.\n", argv[1]);
    		if (debug){
			method_choice = NonIterative;
			method_noniter = LUDecomp;
			sparse_use=0;
		}
    mna_analysis();
    solve_dc();
    if ( do_transient ) {
			printf("[+] Performing Transient analysis\n");
      transient_analysis();
		}
		
		if ( do_dc_instruction) {
			printf("[+] Performing DC instruction\n");
			dc_instruction();
		}
		
		mna_free();
  }

  components_cleanup();
  fclose(yyin);
  yylex_destroy();
  hash_cleanup();




  return mna_size; 
}
