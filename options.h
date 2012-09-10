#ifndef OPTIONS_H
#define OPTIONS_H

#define ITERATIVE_UNIMPLENTED assert(0 && "Iterative Methods are not implemented");
#define SPARSE_UNIMPLEMENTED assert(0 && "Sparse matrices are not implemented");

enum IterativeMethods   {BiCG, CG};
enum SolutionMethods    {Iterative, NonIterative};
enum NonIterativeMethods{LUDecomp, CholDecomp};
enum TransientMethods {Tr, Be};


extern int  sparse_use;
extern enum TransientMethods method_tran;
extern enum IterativeMethods method_iter;
extern enum SolutionMethods method_choice;
extern enum NonIterativeMethods method_noniter;
extern double itol;
#endif
