#ifndef COMPONENTS_H
#define COMPONENTS_H

enum TransientType { Sin, Pwl, Pulse, Exp };

typedef struct V_T i_t;
typedef struct R_T l_t;
typedef struct R_T c_t;

typedef struct PAIR_T
{
  double t, i;
} pair_t;

typedef struct SIN_T
{
  double i1, ia, fr, td, df, ph;
} sin_t;

typedef struct EXP_T
{
  double i1, i2, td1, tc1, td2, tc2;
} exp_t;

typedef struct PULSE_T
{
  double i1, i2, td, tr, tf, pw, per;
} pulse_t;

typedef struct PWL_T
{
  pair_t *pairs;
  int size;
} pwl_t;

typedef struct TRANSIENT_T
{
  enum TransientType type;
  union {
    sin_t tsin;
    exp_t texp;
    pulse_t tpulse;
    pwl_t tpwl;
  };
} transient_t;

typedef struct V_T
{
  int plus, minus;
  double val;
  int id;
  transient_t *transient;

  struct V_T* next;
} v_t;

typedef struct R_T
{
  int plus, minus;
  double val;
  int id;

  struct R_T *next;
} r_t;

extern int voltages;
extern int currents;
extern int resistors;
extern int capacitors;
extern int inductors;


extern i_t *p_i;
extern v_t *p_v;
extern r_t *p_r;
extern l_t *p_l;
extern c_t *p_c;


void new_v(int plus, int minus, double value, transient_t *transient);
void new_i(int plus, int minus, double value, transient_t *transient);
void new_r(int plus, int minus, double value);
void new_c(int plus, int minus, double value);
void new_l(int plus, int minus, double value);

void components_cleanup();



#endif
