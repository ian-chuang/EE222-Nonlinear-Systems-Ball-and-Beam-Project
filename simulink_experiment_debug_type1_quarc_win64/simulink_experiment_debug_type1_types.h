/*
 * simulink_experiment_debug_type1_types.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "simulink_experiment_debug_type1".
 *
 * Model version              : 13.2
 * Simulink Coder version : 9.8 (R2022b) 13-May-2022
 * C source code generated on : Wed Apr 30 13:53:53 2025
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
#ifndef struct_tag_sR4cbkrteQY0cAdcJwPknRE
#define struct_tag_sR4cbkrteQY0cAdcJwPknRE

struct tag_sR4cbkrteQY0cAdcJwPknRE
{
  real_T p_ball_ref;
  real_T v_ball_ref;
  real_T a_ball_ref;
};

#endif                                 /* struct_tag_sR4cbkrteQY0cAdcJwPknRE */

#ifndef typedef_sR4cbkrteQY0cAdcJwPknRE_simul_T
#define typedef_sR4cbkrteQY0cAdcJwPknRE_simul_T

typedef struct tag_sR4cbkrteQY0cAdcJwPknRE sR4cbkrteQY0cAdcJwPknRE_simul_T;

#endif                             /* typedef_sR4cbkrteQY0cAdcJwPknRE_simul_T */

#ifndef struct_tag_AtXt6DhukLcigBVlCE7LlH
#define struct_tag_AtXt6DhukLcigBVlCE7LlH

struct tag_AtXt6DhukLcigBVlCE7LlH
{
  sR4cbkrteQY0cAdcJwPknRE_simul_T workspace;
};

#endif                                 /* struct_tag_AtXt6DhukLcigBVlCE7LlH */

#ifndef typedef_anonymous_function_1_simulink_T
#define typedef_anonymous_function_1_simulink_T

typedef struct tag_AtXt6DhukLcigBVlCE7LlH anonymous_function_1_simulink_T;

#endif                             /* typedef_anonymous_function_1_simulink_T */

#ifndef struct_tag_MmHWQKjO5LbhrLk0ODk58G
#define struct_tag_MmHWQKjO5LbhrLk0ODk58G

struct tag_MmHWQKjO5LbhrLk0ODk58G
{
  anonymous_function_1_simulink_T fun;
};

#endif                                 /* struct_tag_MmHWQKjO5LbhrLk0ODk58G */

#ifndef typedef_s_MmHWQKjO5LbhrLk0ODk58G_simu_T
#define typedef_s_MmHWQKjO5LbhrLk0ODk58G_simu_T

typedef struct tag_MmHWQKjO5LbhrLk0ODk58G s_MmHWQKjO5LbhrLk0ODk58G_simu_T;

#endif                             /* typedef_s_MmHWQKjO5LbhrLk0ODk58G_simu_T */

#ifndef struct_tag_Nld0ZqjQXuSxtFp4wj4VmD
#define struct_tag_Nld0ZqjQXuSxtFp4wj4VmD

struct tag_Nld0ZqjQXuSxtFp4wj4VmD
{
  s_MmHWQKjO5LbhrLk0ODk58G_simu_T workspace;
};

#endif                                 /* struct_tag_Nld0ZqjQXuSxtFp4wj4VmD */

#ifndef typedef_anonymous_function_3_simulink_T
#define typedef_anonymous_function_3_simulink_T

typedef struct tag_Nld0ZqjQXuSxtFp4wj4VmD anonymous_function_3_simulink_T;

#endif                             /* typedef_anonymous_function_3_simulink_T */

#ifndef struct_tag_JDojttTZGKUkoWcqH66oEB
#define struct_tag_JDojttTZGKUkoWcqH66oEB

struct tag_JDojttTZGKUkoWcqH66oEB
{
  int32_T isInitialized;
  real_T t_prev;
  real_T dt;
  real_T rg_val;
  real_T L_val;
  real_T g_val;
  real_T K_val;
  real_T tau_val;
  real_T x_hat[4];
  real_T u;
  char_T controller[6];
  char_T observer[3];
  real_T x_eq[5];
  real_T xh[4];
  real_T P[16];
  real_T xhat_prev[4];
  real_T Q[16];
  real_T R[4];
  real_T Q_tvlqr[16];
  real_T R_tvlqr;
  real_T X_prior[4];
  real_T V_servo;
  real_T initialState[4];
};

#endif                                 /* struct_tag_JDojttTZGKUkoWcqH66oEB */

#ifndef typedef_studentControllerInterface_si_T
#define typedef_studentControllerInterface_si_T

typedef struct tag_JDojttTZGKUkoWcqH66oEB studentControllerInterface_si_T;

#endif                             /* typedef_studentControllerInterface_si_T */

#ifndef struct_tag_xNPby8EIILqaDjffu3lQ3E
#define struct_tag_xNPby8EIILqaDjffu3lQ3E

struct tag_xNPby8EIILqaDjffu3lQ3E
{
  anonymous_function_3_simulink_T nonlin;
  real_T f_1;
  real_T cEq_1;
  real_T f_2;
  real_T cEq_2;
  int32_T nVar;
  int32_T mIneq;
  int32_T mEq;
  int32_T numEvals;
  boolean_T SpecifyObjectiveGradient;
  boolean_T SpecifyConstraintGradient;
  boolean_T isEmptyNonlcon;
  boolean_T hasLB[3];
  boolean_T hasUB[3];
  boolean_T hasBounds;
  int32_T FiniteDifferenceType;
};

#endif                                 /* struct_tag_xNPby8EIILqaDjffu3lQ3E */

#ifndef typedef_s_xNPby8EIILqaDjffu3lQ3E_simu_T
#define typedef_s_xNPby8EIILqaDjffu3lQ3E_simu_T

typedef struct tag_xNPby8EIILqaDjffu3lQ3E s_xNPby8EIILqaDjffu3lQ3E_simu_T;

#endif                             /* typedef_s_xNPby8EIILqaDjffu3lQ3E_simu_T */

/* Parameters (default storage) */
typedef struct P_simulink_experiment_debug_t_T_ P_simulink_experiment_debug_t_T;

/* Forward declaration for rtModel */
typedef struct tag_RTM_simulink_experiment_d_T RT_MODEL_simulink_experiment__T;

#endif                 /* RTW_HEADER_simulink_experiment_debug_type1_types_h_ */
