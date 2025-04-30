/*
 * simulink_experiment_debug_type1_types.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "simulink_experiment_debug_type1".
 *
 * Model version              : 13.1
 * Simulink Coder version : 9.8 (R2022b) 13-May-2022
 * C source code generated on : Wed Apr 30 14:02:04 2025
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_simulink_experiment_debug_type1_types_h_
#define RTW_HEADER_simulink_experiment_debug_type1_types_h_
#include "rtwtypes.h"
#ifndef struct_tag_tumO1kQxtCUPQMu2jqDpTD
#define struct_tag_tumO1kQxtCUPQMu2jqDpTD

struct tag_tumO1kQxtCUPQMu2jqDpTD
{
  int32_T isInitialized;
  real_T rg_val;
  real_T L_val;
  real_T g_val;
  real_T K_val;
  real_T tau_val;
  real_T dt;
  real_T max_V;
  real_T K[4];
  real_T xh[4];
  real_T t_prev;
  real_T u;
  real_T P[16];
  real_T xhat_prev[4];
  real_T Q[16];
  real_T R[4];
};

#endif                                 /* struct_tag_tumO1kQxtCUPQMu2jqDpTD */

#ifndef typedef_studentControllerInterface_si_T
#define typedef_studentControllerInterface_si_T

typedef struct tag_tumO1kQxtCUPQMu2jqDpTD studentControllerInterface_si_T;

#endif                             /* typedef_studentControllerInterface_si_T */

/* Parameters (default storage) */
typedef struct P_simulink_experiment_debug_t_T_ P_simulink_experiment_debug_t_T;

/* Forward declaration for rtModel */
typedef struct tag_RTM_simulink_experiment_d_T RT_MODEL_simulink_experiment__T;

#endif                 /* RTW_HEADER_simulink_experiment_debug_type1_types_h_ */
