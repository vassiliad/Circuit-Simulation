#ifndef MNA_H
#define MNA_H
double g_read(int row, int col);
double c_read(int row, int col);
void g_write(int row, int col, double val);
void c_write(int row, int col, double val);
void c_add(int row, int col, double val);
void g_add(int row, int col, double val);

void mna_free();
void mna_analysis();
void solve_dc();


FILE *f;
char *name_of_file;
#endif
