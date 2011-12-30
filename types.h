#ifndef __TYPES_H
#define __TYPES_H
#include <stdio.h>
#include <stdlib.h>

typedef struct NUMBER_T number_t;
typedef struct TRANSIENT_SPEC_T transient_spec_t;

/*gia V,I,R,C,L*/
struct T1_t {
		int id, original_id;
		int plus, minus;
		double val;
		int is_ground;
		transient_spec_t *transient;
};

/*gia tis diodous*/
struct T2_t {
		int id;
		int plus, minus;
		double area;
		int area_used;
		char *model_name;
};

/*gia ta transistor CMOS*/
struct T3_t {
	int d,g,s,b, id;
	double l,w ;
	char *model_name;
};

/* gia ta BJT */
struct T4_t {
	int id;
	int c,b,e, area_used;
	double area;
	char *model_name;
};

struct description_t {
		int type;
		struct T1_t t1;	// V C I R L
		struct T2_t t2; // D
		struct T3_t t3; // M
		struct T4_t t4; // Q
};

struct components_t {	// diplh sundemenh lista pou krataei posa stoixeia uparxoun ta kuklwma
		struct description_t data;
		struct components_t *next, *prev;
};

enum INSTRUCTION_TYPE { Dc, Plot };
enum SourceType { Voltage, Current };
struct instruction_t {
	struct instruction_t *next, *prev;
	int type;
	union {
		struct _DC_t {
			int sourceType;
			int source;
			double begin;
			double end;
			double inc;
		} dc;

		struct _PLOT_t {
			int *list;
			int num;
			FILE **output;
		} plot;
	};
};

struct ENTRIES_T {
	struct instruction_t *instructions;
	struct components_t  *components;
};

typedef struct ENTRIES_T entries_t;

enum NumberType { Integer, Double };
struct NUMBER_T {
	int type;
	union {
		int integer;
		double dbl;
	};
};

enum Options { SPD, ITER, ITOL };
enum IterType{ NoIter=-1, CG, BiCG };
struct option_t {
	struct option_t *next, *prev;
  int type;
  union {
    double itol;
    int iter_type;
  };
};


typedef struct EXP_T exp_t;
struct EXP_T
{
	double i1, i2;
	double td1, td2;
	double tc1, tc2;
};

typedef struct SIN_T sin_t;
struct SIN_T
{
	double i1;
	double ia, fr, td, df, ph;
};

typedef struct PULSE_T pulse_t;
struct PULSE_T
{
	double i1, i2;
	double td, tr, tf, pw, per;
};

typedef struct PAIR_T pair_t;
struct PAIR_T
{
	double i,t;
};

typedef struct PWL_T pwl_t;
struct PWL_T
{
	pair_t *pairs;
	int size;
};

enum TransientSpecType { Exp, Sin, Pulse, Pwl };
struct TRANSIENT_SPEC_T
{
	enum TransientSpecType type;
	union 
	{
		exp_t exp;
		sin_t _sin;
		pwl_t pwl;
		pulse_t pulse;
	};

};
#endif
