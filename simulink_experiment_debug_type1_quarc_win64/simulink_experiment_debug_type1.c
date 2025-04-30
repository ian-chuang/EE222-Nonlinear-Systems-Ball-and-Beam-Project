/*
 * simulink_experiment_debug_type1.c
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

#include "simulink_experiment_debug_type1.h"
#include "simulink_experiment_debug_type1_types.h"
#include "rtwtypes.h"
#include <string.h>
#include <math.h>
#include "rt_nonfinite.h"
#include "simulink_experiment_debug_type1_private.h"
#include <emmintrin.h>
#include "simulink_experiment_debug_type1_dt.h"

/* Block signals (default storage) */
B_simulink_experiment_debug_t_T simulink_experiment_debug_typ_B;

/* Block states (default storage) */
DW_simulink_experiment_debug__T simulink_experiment_debug_ty_DW;

/* Real-time model */
static RT_MODEL_simulink_experiment__T simulink_experiment_debug_ty_M_;
RT_MODEL_simulink_experiment__T *const simulink_experiment_debug_ty_M =
  &simulink_experiment_debug_ty_M_;

/* Forward declaration for local functions */
static studentControllerInterface_si_T *studentControllerInterface_stud
  (studentControllerInterface_si_T *obj);
static void simulink_e_finDiffEvalAndChkErr(real_T
  obj_nonlin_workspace_fun_worksp, real_T obj_nonlin_workspace_fun_work_0,
  real_T obj_nonlin_workspace_fun_work_1, int32_T dim, real_T delta, const
  real_T xk[3], boolean_T *evalOK, real_T *b_fplus, real_T *b_cEqPlus, real_T
  b_xk[3]);
static void simu_computeFiniteDifferences_j(real_T
  obj_nonlin_workspace_fun_worksp, real_T obj_nonlin_workspace_fun_work_0,
  real_T obj_nonlin_workspace_fun_work_1, real_T obj_f_1, real_T obj_cEq_1,
  real_T cEqCurrent, const real_T xk[3], real_T JacCeqTrans[3], boolean_T
  *evalOK, real_T b_xk[3], s_xNPby8EIILqaDjffu3lQ3E_simu_T *b_obj);
static real_T simulink_experiment_debug__norm(const real_T x[3]);
static real_T simulink_experiment_debug_xnrm2(int32_T n, const real_T x[12],
  int32_T ix0);
static void simulink_experiment_deb_xzlarfg(int32_T n, real_T alpha1, const
  real_T x[12], int32_T ix0, real_T *b_alpha1, real_T b_x[12], real_T *tau);
static void simulink_experiment_de_xzlarf_j(int32_T m, int32_T n, int32_T iv0,
  real_T tau, const real_T C[12], int32_T ic0, real_T work[3], real_T b_C[12]);
static void simulink_e_linearLeastSquares_j(const real_T lhs[12], real_T rhs[4],
  real_T b_lhs[12], real_T b_dx[3]);
static void simulink_experiment_debu_fsolve(real_T fun_workspace_p_ball_ref,
  real_T fun_workspace_v_ball_ref, real_T fun_workspace_a_ball_ref, const real_T
  x[3], real_T b_x[3]);
static void simulink_experiment_debug_t_inv(const real_T x[64], real_T y[64]);
static real_T simulink_experiment_deb_xnrm2_j(int32_T n, const real_T x[32],
  int32_T ix0);
static void simulink_experiment_de_mldivide(const real_T A[32], const real_T B
  [32], real_T Y[16]);
static void studentControllerInterface_step(studentControllerInterface_si_T *obj,
  real_T t, real_T p_ball, real_T theta, real_T *V_servo, real_T x_hat[16]);
static void rate_monotonic_scheduler(void);
time_T rt_SimUpdateDiscreteEvents(
  int_T rtmNumSampTimes, void *rtmTimingData, int_T *rtmSampleHitPtr, int_T
  *rtmPerTaskSampleHits )
{
  rtmSampleHitPtr[1] = rtmStepTask(simulink_experiment_debug_ty_M, 1);
  rtmSampleHitPtr[2] = rtmStepTask(simulink_experiment_debug_ty_M, 2);
  UNUSED_PARAMETER(rtmNumSampTimes);
  UNUSED_PARAMETER(rtmTimingData);
  UNUSED_PARAMETER(rtmPerTaskSampleHits);
  return(-1);
}

/*
 *         This function updates active task flag for each subrate
 *         and rate transition flags for tasks that exchange data.
 *         The function assumes rate-monotonic multitasking scheduler.
 *         The function must be called at model base rate so that
 *         the generated code self-manages all its subrates and rate
 *         transition flags.
 */
static void rate_monotonic_scheduler(void)
{
  /* To ensure a deterministic data transfer between two rates,
   * data is transferred at the priority of a fast task and the frequency
   * of the slow task.  The following flags indicate when the data transfer
   * happens.  That is, a rate interaction flag is set true when both rates
   * will run, and false otherwise.
   */

  /* tid 1 shares data with slower tid rate: 2 */
  if (simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[1] == 0) {
    simulink_experiment_debug_ty_M->Timing.RateInteraction.TID1_2 =
      (simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[2] == 0);

    /* update PerTaskSampleHits matrix for non-inline sfcn */
    simulink_experiment_debug_ty_M->Timing.perTaskSampleHits[5] =
      simulink_experiment_debug_ty_M->Timing.RateInteraction.TID1_2;
  }

  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[2])++;
  if ((simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[2]) > 4) {/* Sample time: [0.01s, 0.0s] */
    simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[2] = 0;
  }
}

static studentControllerInterface_si_T *studentControllerInterface_stud
  (studentControllerInterface_si_T *obj)
{
  studentControllerInterface_si_T *b_obj;
  int32_T i;
  static const char_T tmp[6] = { 'T', 'V', '-', 'L', 'Q', 'R' };

  static const real_T tmp_0[16] = { 0.5, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0,
    0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.01 };

  static const real_T tmp_1[16] = { 0.01, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0,
    0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.01 };

  b_obj = obj;
  b_obj->t_prev = -1.0;
  b_obj->rg_val = 0.0254;
  b_obj->L_val = 0.4255;
  b_obj->g_val = 9.81;
  b_obj->K_val = 1.5;
  b_obj->tau_val = 0.025;
  b_obj->x_hat[0] = 0.0;
  b_obj->x_hat[1] = 0.0;
  b_obj->x_hat[2] = 0.0;
  b_obj->x_hat[3] = 0.0;
  b_obj->u = 0.0;
  for (i = 0; i < 6; i++) {
    b_obj->controller[i] = tmp[i];
  }

  b_obj->observer[0] = 'E';
  b_obj->observer[1] = 'K';
  b_obj->observer[2] = 'F';
  for (i = 0; i < 5; i++) {
    b_obj->x_eq[i] = 0.0;
  }

  for (i = 0; i < 16; i++) {
    b_obj->P[i] = tmp_0[i];
  }

  b_obj->xhat_prev[0] = 0.0;
  b_obj->xhat_prev[1] = 0.0;
  b_obj->xhat_prev[2] = -0.99483767363676778;
  b_obj->xhat_prev[3] = 0.0;
  for (i = 0; i < 16; i++) {
    b_obj->Q[i] = tmp_1[i];
  }

  b_obj->R[0] = 0.003;
  b_obj->R[1] = 0.0;
  b_obj->R[2] = 0.0;
  b_obj->R[3] = 0.003;
  b_obj->isInitialized = 0;

  /*  obj.controller = 'TV-LQR'; */
  /*  obj.observer = 'ELO'; */
  /* %%%%% EKF %%%%%%%% */
  /*  discrete_f = @(x, u) x + obj.dt * [ x(2);... */
  /*  ((5 * obj.g_val * obj.rg_val)/(7 * obj.L_val)) * sin(x(3)) - (5/7) * ((obj.L_val/2) - x(1)) * ((obj.rg_val/obj.L_val) * x(4))^2 * cos(x(3))^2;... */
  /*  x(4); -1 * x(4)/obj.tau_val + (obj.K_val/obj.tau_val) * u ]; */
  /*   */
  /*  syms x1 x2 x3 x4 u_sym dt_sym real */
  /*   */
  /*  x_sym = [x1; x2; x3; x4]; */
  /*   */
  /*  f1_sym = x1 + dt_sym * x2; */
  /*  f2_sym = x2 + dt_sym * (((5 * obj.g_val * obj.rg_val)/(7 * obj.L_val)) * sin(x3) - (5/7) * ((obj.L_val/2) - x1) * ((obj.rg_val/obj.L_val) * x4)^2 * cos(x3)^2); */
  /*  f3_sym = x3 + dt_sym * x4; */
  /*  f4_sym = x4 + dt_sym * (-1 * x4/obj.tau_val + (obj.K_val/obj.tau_val) * u_sym); */
  /*  f_sym = [f1_sym; f2_sym; f3_sym; f4_sym]; */
  /*   */
  /*  J_sym = jacobian(f_sym, x_sym); */
  /*   */
  /*  J_handle = matlabFunction(J_sym, 'Vars', {[x1; x2; x3; x4], u_sym, dt_sym}); */
  /*   */
  /*  f_Jacobian = @(x, u) J_handle(x, u, obj.dt); */
  /*   */
  /*  measurement_z = @(x) [x(1); x(3)]; */
  /*   */
  /*  z_Jacobian = @(x) [1 0 0 0; 0 0 1 0]; */
  /*   */
  /*  obj.ekf = extendedKalmanFilter(discrete_f, measurement_z, obj.initialState, ... */
  /*      'StateTransitionJacobianFcn', f_Jacobian, ... */
  /*      'MeasurementJacobianFcn', z_Jacobian); */
  /*   */
  /*  obj.ekf.ProcessNoise = diag([1e-3, 1e-2, 1e-3, 1e-3]); */
  /*  obj.ekf.MeasurementNoise = diag([1e-6, 1e-6]); */
  /* %%%%% MHE %%%%%%%% */
  /*  obj.setupDynamics(); */
  /*  obj.setupMHE(); */
  return b_obj;
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T tmp;
  real_T tmp_0;
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    tmp = fabs(u0);
    tmp_0 = fabs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = 1.0;
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = (rtNaN);
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

static void simulink_e_finDiffEvalAndChkErr(real_T
  obj_nonlin_workspace_fun_worksp, real_T obj_nonlin_workspace_fun_work_0,
  real_T obj_nonlin_workspace_fun_work_1, int32_T dim, real_T delta, const
  real_T xk[3], boolean_T *evalOK, real_T *b_fplus, real_T *b_cEqPlus, real_T
  b_xk[3])
{
  real_T a;
  real_T f_sym_continuous_idx_0;
  real_T temp;
  real_T x1;
  real_T x2;
  real_T x3;
  real_T x4;
  real_T x_idx_2;
  real_T y;
  int32_T idx;
  boolean_T b;
  boolean_T c;
  b_xk[0] = xk[0];
  b_xk[1] = xk[1];
  b_xk[2] = xk[2];
  *b_fplus = 0.0;
  *evalOK = true;
  temp = b_xk[dim - 1];
  b_xk[dim - 1] += delta;
  x3 = b_xk[0];
  x4 = b_xk[1];
  x_idx_2 = b_xk[2];
  f_sym_continuous_idx_0 = obj_nonlin_workspace_fun_worksp;
  x2 = obj_nonlin_workspace_fun_work_0;

  /* F_FN */
  /*     F_SYM_CONTINUOUS = F_FN(IN1,U_SYM) */
  /*     This function was generated by the Symbolic Math Toolbox version 9.2. */
  /*     30-Apr-2025 09:26:50 */
  x1 = f_sym_continuous_idx_0;
  y = x4 * x4;
  a = x3;
  a = cos(a);
  a *= a;
  x3 = sin(x3);
  f_sym_continuous_idx_0 = x2;
  x2 = (x1 * 0.7142857142857143 - 0.15196428571428569) * (y * a) *
    0.0035634305945448849 + x3 * 0.41828772872251141;
  x3 = x4;
  x4 = x_idx_2 * 60.0 - x4 * 40.0;
  *b_cEqPlus = f_sym_continuous_idx_0 - obj_nonlin_workspace_fun_work_0;
  x2 -= obj_nonlin_workspace_fun_work_1;
  *b_cEqPlus += x2;
  *b_cEqPlus += 0.0 * x3;
  *b_cEqPlus += 0.0 * x4;
  idx = 1;
  while ((*evalOK) && (idx <= 1)) {
    b = rtIsInf(*b_cEqPlus);
    c = !b;
    b = rtIsNaN(*b_cEqPlus);
    b = !b;
    *evalOK = (c && b);
    idx = 2;
  }

  b_xk[dim - 1] = temp;
}

static void simu_computeFiniteDifferences_j(real_T
  obj_nonlin_workspace_fun_worksp, real_T obj_nonlin_workspace_fun_work_0,
  real_T obj_nonlin_workspace_fun_work_1, real_T obj_f_1, real_T obj_cEq_1,
  real_T cEqCurrent, const real_T xk[3], real_T JacCeqTrans[3], boolean_T
  *evalOK, real_T b_xk[3], s_xNPby8EIILqaDjffu3lQ3E_simu_T *b_obj)
{
  real_T xk_0[3];
  real_T obj_f_1_0;
  real_T obj_nonlin_workspace_fun_work_2;
  real_T obj_nonlin_workspace_fun_work_3;
  real_T obj_nonlin_workspace_fun_work_4;
  real_T signVec;
  real_T varargin_2;
  int32_T b_idx;
  int32_T c_obj_numEvals;
  int32_T i;
  int32_T idx;
  boolean_T exitg1;
  boolean_T guard1;
  b_obj->nonlin.workspace.fun.workspace.p_ball_ref =
    obj_nonlin_workspace_fun_worksp;
  b_obj->nonlin.workspace.fun.workspace.v_ball_ref =
    obj_nonlin_workspace_fun_work_0;
  b_obj->nonlin.workspace.fun.workspace.a_ball_ref =
    obj_nonlin_workspace_fun_work_1;
  b_obj->f_1 = obj_f_1;
  b_obj->cEq_1 = obj_cEq_1;
  b_xk[0] = xk[0];
  b_xk[1] = xk[1];
  b_xk[2] = xk[2];
  obj_nonlin_workspace_fun_work_2 =
    b_obj->nonlin.workspace.fun.workspace.p_ball_ref;
  obj_nonlin_workspace_fun_work_3 =
    b_obj->nonlin.workspace.fun.workspace.v_ball_ref;
  obj_nonlin_workspace_fun_work_4 =
    b_obj->nonlin.workspace.fun.workspace.a_ball_ref;
  obj_f_1_0 = b_obj->f_1;
  varargin_2 = b_obj->cEq_1;
  *evalOK = true;
  c_obj_numEvals = 0;
  b_idx = 1;
  exitg1 = false;
  while ((!exitg1) && (b_idx - 1 < 3)) {
    idx = b_idx - 1;
    xk_0[0] = b_xk[0];
    xk_0[1] = b_xk[1];
    xk_0[2] = b_xk[2];
    signVec = 1.0 - (real_T)(xk_0[idx] < 0.0) * 2.0;
    obj_f_1_0 = xk_0[idx];
    obj_f_1_0 = fabs(obj_f_1_0);
    if (!(obj_f_1_0 >= 1.0)) {
      obj_f_1_0 = 1.0;
    }

    signVec = 1.4901161193847656E-8 * signVec * obj_f_1_0;
    for (i = 0; i < 3; i++) {
      xk_0[i] = b_xk[i];
    }

    simulink_e_finDiffEvalAndChkErr(obj_nonlin_workspace_fun_work_2,
      obj_nonlin_workspace_fun_work_3, obj_nonlin_workspace_fun_work_4, idx + 1,
      signVec, xk_0, evalOK, &obj_f_1_0, &varargin_2, b_xk);
    obj_f_1_0 = 0.0;
    c_obj_numEvals++;
    guard1 = false;
    if (!*evalOK) {
      signVec = -signVec;
      for (i = 0; i < 3; i++) {
        xk_0[i] = b_xk[i];
      }

      simulink_e_finDiffEvalAndChkErr(obj_nonlin_workspace_fun_work_2,
        obj_nonlin_workspace_fun_work_3, obj_nonlin_workspace_fun_work_4, idx +
        1, signVec, xk_0, evalOK, &obj_f_1_0, &varargin_2, b_xk);
      obj_f_1_0 = 0.0;
      c_obj_numEvals++;
      if (!*evalOK) {
        exitg1 = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      JacCeqTrans[idx] = (varargin_2 - cEqCurrent) / signVec;
      b_idx++;
    }
  }

  b_obj->nonlin.workspace.fun.workspace.p_ball_ref =
    obj_nonlin_workspace_fun_work_2;
  b_obj->nonlin.workspace.fun.workspace.v_ball_ref =
    obj_nonlin_workspace_fun_work_3;
  b_obj->nonlin.workspace.fun.workspace.a_ball_ref =
    obj_nonlin_workspace_fun_work_4;
  b_obj->f_1 = obj_f_1_0;
  b_obj->cEq_1 = varargin_2;
  b_obj->numEvals = c_obj_numEvals;
}

static real_T simulink_experiment_debug__norm(const real_T x[3])
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  scale = 3.3121686421112381E-170;
  absxk = x[0];
  absxk = fabs(absxk);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }

  absxk = x[1];
  absxk = fabs(absxk);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  absxk = x[2];
  absxk = fabs(absxk);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  return scale * sqrt(y);
}

static real_T simulink_experiment_debug_xnrm2(int32_T n, const real_T x[12],
  int32_T ix0)
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  int32_T b;
  int32_T kend;
  y = 0.0;
  if (n == 1) {
    absxk = x[ix0 - 1];
    y = fabs(absxk);
  } else {
    scale = 3.3121686421112381E-170;
    kend = n;
    kend--;
    kend += ix0;
    for (b = ix0; b <= kend; b++) {
      absxk = x[b - 1];
      absxk = fabs(absxk);
      if (absxk > scale) {
        t = scale / absxk;
        y = y * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
    }

    y = scale * sqrt(y);
  }

  return y;
}

real_T rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T a;
  real_T b;
  real_T y;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = sqrt(a * a + 1.0) * b;
  } else if (a > b) {
    b /= a;
    y = sqrt(b * b + 1.0) * a;
  } else if (rtIsNaN(b)) {
    y = (rtNaN);
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

static void simulink_experiment_deb_xzlarfg(int32_T n, real_T alpha1, const
  real_T x[12], int32_T ix0, real_T *b_alpha1, real_T b_x[12], real_T *tau)
{
  __m128d tmp;
  real_T xnorm;
  real_T y;
  int32_T a;
  int32_T c;
  int32_T knt;
  int32_T nm1;
  int32_T scalarLB;
  int32_T vectorUB;
  memcpy(&b_x[0], &x[0], 12U * sizeof(real_T));
  *b_alpha1 = alpha1;
  *tau = 0.0;
  nm1 = n - 2;
  xnorm = simulink_experiment_debug_xnrm2(nm1 + 1, b_x, ix0);
  if (xnorm != 0.0) {
    xnorm = rt_hypotd_snf(*b_alpha1, xnorm);
    if (*b_alpha1 >= 0.0) {
      xnorm = -xnorm;
    }

    y = fabs(xnorm);
    if (y < 1.0020841800044864E-292) {
      knt = -1;
      do {
        knt++;
        c = nm1;
        c += ix0;
        scalarLB = ((((c - ix0) + 1) / 2) << 1) + ix0;
        vectorUB = scalarLB - 2;
        for (a = ix0; a <= vectorUB; a += 2) {
          tmp = _mm_loadu_pd(&b_x[a - 1]);
          tmp = _mm_mul_pd(tmp, _mm_set1_pd(9.9792015476736E+291));
          _mm_storeu_pd(&b_x[a - 1], tmp);
        }

        for (a = scalarLB; a <= c; a++) {
          b_x[a - 1] *= 9.9792015476736E+291;
        }

        xnorm *= 9.9792015476736E+291;
        *b_alpha1 *= 9.9792015476736E+291;
        y = fabs(xnorm);
      } while ((y < 1.0020841800044864E-292) && (knt + 1 < 20));

      xnorm = simulink_experiment_debug_xnrm2(nm1 + 1, b_x, ix0);
      xnorm = rt_hypotd_snf(*b_alpha1, xnorm);
      if (*b_alpha1 >= 0.0) {
        xnorm = -xnorm;
      }

      y = xnorm - *b_alpha1;
      *tau = y / xnorm;
      y = *b_alpha1 - xnorm;
      *b_alpha1 = 1.0 / y;
      c = nm1;
      c += ix0;
      scalarLB = ((((c - ix0) + 1) / 2) << 1) + ix0;
      vectorUB = scalarLB - 2;
      for (a = ix0; a <= vectorUB; a += 2) {
        tmp = _mm_loadu_pd(&b_x[a - 1]);
        tmp = _mm_mul_pd(tmp, _mm_set1_pd(*b_alpha1));
        _mm_storeu_pd(&b_x[a - 1], tmp);
      }

      for (a = scalarLB; a <= c; a++) {
        b_x[a - 1] *= *b_alpha1;
      }

      for (a = 0; a <= knt; a++) {
        xnorm *= 1.0020841800044864E-292;
      }

      *b_alpha1 = xnorm;
    } else {
      y = xnorm - *b_alpha1;
      *tau = y / xnorm;
      y = *b_alpha1 - xnorm;
      *b_alpha1 = 1.0 / y;
      c = nm1;
      c += ix0;
      scalarLB = ((((c - ix0) + 1) / 2) << 1) + ix0;
      vectorUB = scalarLB - 2;
      for (a = ix0; a <= vectorUB; a += 2) {
        tmp = _mm_loadu_pd(&b_x[a - 1]);
        tmp = _mm_mul_pd(tmp, _mm_set1_pd(*b_alpha1));
        _mm_storeu_pd(&b_x[a - 1], tmp);
      }

      for (a = scalarLB; a <= c; a++) {
        b_x[a - 1] *= *b_alpha1;
      }

      *b_alpha1 = xnorm;
    }
  }
}

static void simulink_experiment_de_xzlarf_j(int32_T m, int32_T n, int32_T iv0,
  real_T tau, const real_T C[12], int32_T ic0, real_T work[3], real_T b_C[12])
{
  real_T A[12];
  real_T x[12];
  real_T b;
  real_T c;
  real_T yjy;
  int32_T colbottom;
  int32_T coltop;
  int32_T e;
  int32_T exitg1;
  int32_T ix;
  int32_T jy;
  int32_T lastc;
  int32_T lastv;
  int32_T mm1;
  int32_T nm1;
  boolean_T exitg2;
  memcpy(&b_C[0], &C[0], 12U * sizeof(real_T));
  if (tau != 0.0) {
    lastv = m;
    colbottom = m;
    lastc = (iv0 + colbottom) - 2;
    while ((lastv > 0) && (b_C[lastc] == 0.0)) {
      lastv--;
      lastc--;
    }

    lastc = n;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      colbottom = lastc;
      colbottom = (colbottom - 1) << 2;
      coltop = ic0 + colbottom;
      colbottom = lastv;
      colbottom = (coltop + colbottom) - 1;
      do {
        exitg1 = 0;
        if (coltop <= colbottom) {
          if (b_C[coltop - 1] != 0.0) {
            exitg1 = 1;
          } else {
            coltop++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = 0;
  }

  if (lastv > 0) {
    for (colbottom = 0; colbottom < 12; colbottom++) {
      c = b_C[colbottom];
      x[colbottom] = c;
      A[colbottom] = c;
    }

    if (lastc != 0) {
      mm1 = lastv;
      nm1 = lastc - 1;
      colbottom = nm1;
      coltop = colbottom + 1;
      coltop = (uint8_T)coltop - 1;
      for (colbottom = 0; colbottom <= coltop; colbottom++) {
        jy = colbottom;
        work[jy] = 0.0;
      }

      jy = 0;
      colbottom = nm1 << 2;
      nm1 = ic0 + colbottom;
      for (colbottom = ic0; colbottom <= nm1; colbottom += 4) {
        ix = iv0 - 1;
        c = 0.0;
        e = (colbottom + mm1) - 1;
        for (coltop = colbottom; coltop <= e; coltop++) {
          yjy = A[coltop - 1];
          b = x[ix];
          yjy *= b;
          c += yjy;
          ix++;
        }

        work[jy] += c;
        jy++;
      }
    }

    c = -tau;
    if (!(c == 0.0)) {
      mm1 = ic0;
      jy = 0;
      coltop = (uint8_T)lastc - 1;
      for (lastc = 0; lastc <= coltop; lastc++) {
        yjy = work[jy];
        if (yjy != 0.0) {
          yjy *= c;
          ix = iv0 - 1;
          colbottom = mm1;
          nm1 = (lastv + mm1) - 1;
          for (e = colbottom; e <= nm1; e++) {
            b_C[e - 1] += b_C[ix] * yjy;
            ix++;
          }
        }

        jy++;
        mm1 += 4;
      }
    }
  }
}

static void simulink_e_linearLeastSquares_j(const real_T lhs[12], real_T rhs[4],
  real_T b_lhs[12], real_T b_dx[3])
{
  __m128d tmp;
  __m128d tmp_0;
  real_T b_lhs_0[12];
  real_T b_work[3];
  real_T c_jpvt[3];
  real_T tau[3];
  real_T vn1[3];
  real_T vn2[3];
  real_T absxk;
  real_T b_atmp;
  real_T t;
  real_T temp;
  int32_T fcol;
  int32_T ii;
  int32_T im1;
  int32_T ip1;
  int32_T ix;
  int32_T jcol;
  int32_T k;
  int32_T mcol;
  int32_T mmip1;
  int32_T nfxd;
  int32_T nmi;
  int32_T nmip1;
  memcpy(&b_lhs[0], &lhs[0], 12U * sizeof(real_T));
  c_jpvt[0] = 0.0;
  c_jpvt[1] = 0.0;
  c_jpvt[2] = 0.0;
  nfxd = 0;
  for (ii = 0; ii < 3; ii++) {
    im1 = ii;
    if (c_jpvt[im1] != 0.0) {
      nfxd++;
      if (im1 + 1 != nfxd) {
        jcol = im1 << 2;
        fcol = (nfxd - 1) << 2;
        temp = b_lhs[jcol];
        b_lhs[jcol] = b_lhs[fcol];
        b_lhs[fcol] = temp;
        jcol++;
        fcol++;
        temp = b_lhs[jcol];
        b_lhs[jcol] = b_lhs[fcol];
        b_lhs[fcol] = temp;
        jcol++;
        fcol++;
        temp = b_lhs[jcol];
        b_lhs[jcol] = b_lhs[fcol];
        b_lhs[fcol] = temp;
        jcol++;
        fcol++;
        temp = b_lhs[jcol];
        b_lhs[jcol] = b_lhs[fcol];
        b_lhs[fcol] = temp;
        c_jpvt[im1] = c_jpvt[nfxd - 1];
        c_jpvt[nfxd - 1] = (real_T)im1 + 1.0;
      } else {
        c_jpvt[im1] = (real_T)im1 + 1.0;
      }
    } else {
      c_jpvt[im1] = (real_T)im1 + 1.0;
    }

    tau[ii] = 0.0;
    b_work[ii] = 0.0;
  }

  k = nfxd - 1;
  for (ip1 = 0; ip1 <= k; ip1++) {
    jcol = ip1;
    im1 = jcol;
    ii = (im1 << 2) + im1;
    nmi = 2 - jcol;
    fcol = 4 - jcol;
    mmip1 = fcol;
    b_atmp = b_lhs[ii];
    memcpy(&b_lhs_0[0], &b_lhs[0], 12U * sizeof(real_T));
    simulink_experiment_deb_xzlarfg(mmip1, b_atmp, b_lhs_0, ii + 2, &absxk,
      b_lhs, &temp);
    tau[jcol] = temp;
    b_lhs[ii] = absxk;
    if (jcol + 1 < 3) {
      b_atmp = b_lhs[ii];
      b_lhs[ii] = 1.0;
      temp = tau[jcol];
      memcpy(&b_lhs_0[0], &b_lhs[0], 12U * sizeof(real_T));
      simulink_experiment_de_xzlarf_j(mmip1, nmi, ii + 1, temp, b_lhs_0, ii + 5,
        b_work, b_lhs);
      b_lhs[ii] = b_atmp;
    }
  }

  if (nfxd < 3) {
    b_work[0] = 0.0;
    vn1[0] = 0.0;
    vn2[0] = 0.0;
    b_work[1] = 0.0;
    vn1[1] = 0.0;
    vn2[1] = 0.0;
    b_work[2] = 0.0;
    vn1[2] = 0.0;
    vn2[2] = 0.0;
    mmip1 = nfxd;
    k = nfxd;
    for (im1 = k + 1; im1 < 4; im1++) {
      vn1[im1 - 1] = simulink_experiment_debug_xnrm2(4 - nfxd, b_lhs, (((im1 - 1)
        << 2) + mmip1) + 1);
      vn2[im1 - 1] = vn1[im1 - 1];
    }

    for (jcol = nfxd + 1; jcol < 4; jcol++) {
      im1 = jcol;
      ip1 = jcol;
      ii = (((im1 - 1) << 2) + im1) - 1;
      nmi = 3 - jcol;
      fcol = 4 - jcol;
      mmip1 = fcol + 1;
      nmip1 = nmi + 1;
      mcol = 1;
      if (nmip1 > 1) {
        ix = jcol - 1;
        absxk = vn1[jcol - 1];
        temp = fabs(absxk);
        b_atmp = temp;
        for (k = 2; k <= nmip1; k++) {
          ix++;
          absxk = vn1[ix];
          temp = fabs(absxk);
          if (temp > b_atmp) {
            mcol = k;
            b_atmp = temp;
          }
        }
      }

      ix = (im1 + mcol) - 2;
      if (ix + 1 != jcol) {
        nmip1 = ix << 2;
        mcol = (im1 - 1) << 2;
        temp = b_lhs[nmip1];
        b_lhs[nmip1] = b_lhs[mcol];
        b_lhs[mcol] = temp;
        nmip1++;
        mcol++;
        temp = b_lhs[nmip1];
        b_lhs[nmip1] = b_lhs[mcol];
        b_lhs[mcol] = temp;
        nmip1++;
        mcol++;
        temp = b_lhs[nmip1];
        b_lhs[nmip1] = b_lhs[mcol];
        b_lhs[mcol] = temp;
        nmip1++;
        mcol++;
        temp = b_lhs[nmip1];
        b_lhs[nmip1] = b_lhs[mcol];
        b_lhs[mcol] = temp;
        temp = c_jpvt[ix];
        c_jpvt[ix] = c_jpvt[jcol - 1];
        c_jpvt[jcol - 1] = temp;
        vn1[ix] = vn1[jcol - 1];
        vn2[ix] = vn2[jcol - 1];
      }

      b_atmp = b_lhs[ii];
      memcpy(&b_lhs_0[0], &b_lhs[0], 12U * sizeof(real_T));
      simulink_experiment_deb_xzlarfg(mmip1, b_atmp, b_lhs_0, ii + 2, &absxk,
        b_lhs, &temp);
      tau[jcol - 1] = temp;
      b_lhs[ii] = absxk;
      if (jcol < 3) {
        b_atmp = b_lhs[ii];
        b_lhs[ii] = 1.0;
        temp = tau[jcol - 1];
        memcpy(&b_lhs_0[0], &b_lhs[0], 12U * sizeof(real_T));
        simulink_experiment_de_xzlarf_j(mmip1, nmi, ii + 1, temp, b_lhs_0, ii +
          5, b_work, b_lhs);
        b_lhs[ii] = b_atmp;
      }

      for (ii = ip1 + 1; ii < 4; ii++) {
        mmip1 = (((ii - 1) << 2) + im1) - 1;
        if (vn1[ii - 1] != 0.0) {
          absxk = b_lhs[mmip1];
          temp = fabs(absxk);
          temp /= vn1[ii - 1];
          temp = 1.0 - temp * temp;
          if (temp < 0.0) {
            temp = 0.0;
          }

          b_atmp = vn1[ii - 1] / vn2[ii - 1];
          b_atmp = b_atmp * b_atmp * temp;
          if (b_atmp <= 1.4901161193847656E-8) {
            nmi = mmip1 + 2;
            temp = 0.0;
            if (fcol == 1) {
              absxk = b_lhs[nmi - 1];
              temp = fabs(absxk);
            } else {
              b_atmp = 3.3121686421112381E-170;
              k = fcol;
              k--;
              mmip1 = nmi + k;
              for (k = nmi; k <= mmip1; k++) {
                absxk = b_lhs[k - 1];
                absxk = fabs(absxk);
                if (absxk > b_atmp) {
                  t = b_atmp / absxk;
                  temp = temp * t * t + 1.0;
                  b_atmp = absxk;
                } else {
                  t = absxk / b_atmp;
                  temp += t * t;
                }
              }

              temp = b_atmp * sqrt(temp);
            }

            vn1[ii - 1] = temp;
            vn2[ii - 1] = vn1[ii - 1];
          } else {
            temp = sqrt(temp);
            vn1[ii - 1] *= temp;
          }
        }
      }
    }
  }

  for (ii = 0; ii < 3; ii++) {
    im1 = ii;
    b_atmp = tau[im1];
    if (b_atmp != 0.0) {
      absxk = rhs[im1];
      mmip1 = im1;
      for (jcol = mmip1 + 2; jcol < 5; jcol++) {
        t = b_lhs[((im1 << 2) + jcol) - 1];
        temp = rhs[jcol - 1];
        temp *= t;
        absxk += temp;
      }

      absxk *= b_atmp;
      if (absxk != 0.0) {
        rhs[im1] -= absxk;
        k = im1;
        jcol = ((((3 - k) / 2) << 1) + k) + 2;
        nfxd = jcol - 2;
        for (ip1 = k + 2; ip1 <= nfxd; ip1 += 2) {
          tmp = _mm_loadu_pd(&b_lhs[((im1 << 2) + ip1) - 1]);
          tmp = _mm_mul_pd(tmp, _mm_set1_pd(absxk));
          tmp_0 = _mm_loadu_pd(&rhs[ip1 - 1]);
          tmp = _mm_sub_pd(tmp_0, tmp);
          _mm_storeu_pd(&rhs[ip1 - 1], tmp);
        }

        for (ip1 = jcol; ip1 < 5; ip1++) {
          rhs[ip1 - 1] -= b_lhs[((im1 << 2) + ip1) - 1] * absxk;
        }
      }
    }
  }

  b_dx[0] = rhs[0];
  b_dx[1] = rhs[1];
  b_dx[2] = rhs[2];
  for (im1 = 2; im1 >= 0; im1--) {
    k = im1 + 1;
    fcol = im1;
    fcol <<= 2;
    nfxd = k + fcol;
    k = im1;
    k++;
    ii = k - 1;
    absxk = b_dx[ii];
    temp = b_lhs[nfxd - 1];
    temp = absxk / temp;
    b_dx[ii] = temp;
    mmip1 = im1;
    k = mmip1 - 1;
    for (ip1 = 0; ip1 <= k; ip1++) {
      jcol = ip1 + 1;
      fcol = jcol;
      ix = ii - fcol;
      fcol = (nfxd - jcol) - 1;
      b_dx[ix] -= b_dx[ii] * b_lhs[fcol];
    }
  }

  temp = b_dx[0];
  b_atmp = b_dx[1];
  absxk = b_dx[2];
  b_dx[(int32_T)c_jpvt[0] - 1] = temp;
  b_dx[(int32_T)c_jpvt[1] - 1] = b_atmp;
  b_dx[(int32_T)c_jpvt[2] - 1] = absxk;
}

static void simulink_experiment_debu_fsolve(real_T fun_workspace_p_ball_ref,
  real_T fun_workspace_v_ball_ref, real_T fun_workspace_a_ball_ref, const real_T
  x[3], real_T b_x[3])
{
  s_xNPby8EIILqaDjffu3lQ3E_simu_T c_FiniteDifferences;
  real_T y[12];
  real_T y_0[12];
  real_T b_rhs[4];
  real_T B[3];
  real_T JacCeqTrans[3];
  real_T a__3[3];
  real_T gradf[3];
  real_T x_0[3];
  real_T xp[3];
  real_T a;
  real_T b_gamma;
  real_T b_x_0;
  real_T c;
  real_T c_FiniteDifferences_cEq_1;
  real_T c_FiniteDifferences_f_1;
  real_T c_FiniteDifferences_nonlin_wo_0;
  real_T c_FiniteDifferences_nonlin_wo_1;
  real_T c_FiniteDifferences_nonlin_work;
  real_T f_sym_continuous_idx_0;
  real_T f_temp;
  real_T funDiff;
  real_T relFactor;
  real_T resnorm;
  real_T x1;
  real_T x2;
  real_T x3;
  real_T x4;
  int32_T aIdx;
  int32_T b_i;
  int32_T c_exitflag;
  int32_T exitg1;
  int32_T funcCount;
  int32_T iter;
  int32_T j;
  boolean_T b;
  boolean_T b_evalOK;
  boolean_T guard1;
  boolean_T stepSuccessful;
  funDiff = (rtInf);
  iter = 0;
  b_x_0 = x[0];
  x_0[0] = b_x_0;
  b_x[0] = b_x_0;
  b_x_0 = x[1];
  x_0[1] = b_x_0;
  b_x[1] = b_x_0;
  b_x_0 = x[2];
  x_0[2] = b_x_0;
  b_x[2] = b_x_0;
  f_sym_continuous_idx_0 = fun_workspace_p_ball_ref;
  x2 = fun_workspace_v_ball_ref;
  x3 = x_0[0];
  x4 = x_0[1];
  f_temp = x_0[2];

  /* F_FN */
  /*     F_SYM_CONTINUOUS = F_FN(IN1,U_SYM) */
  /*     This function was generated by the Symbolic Math Toolbox version 9.2. */
  /*     30-Apr-2025 09:26:50 */
  x1 = f_sym_continuous_idx_0;
  c = x4 * x4;
  a = x3;
  a = cos(a);
  a *= a;
  x3 = sin(x3);
  f_sym_continuous_idx_0 = x2;
  x2 = (x1 * 0.7142857142857143 - 0.15196428571428569) * (c * a) *
    0.0035634305945448849 + x3 * 0.41828772872251141;
  x3 = x4;
  x4 = f_temp * 60.0 - x4 * 40.0;
  f_temp = f_sym_continuous_idx_0 - fun_workspace_v_ball_ref;
  f_sym_continuous_idx_0 = x2 - fun_workspace_a_ball_ref;
  f_temp += f_sym_continuous_idx_0;
  f_temp += 0.0 * x3;
  f_temp += 0.0 * x4;
  resnorm = f_temp * f_temp;
  c_FiniteDifferences_nonlin_work = fun_workspace_p_ball_ref;
  c_FiniteDifferences_nonlin_wo_0 = fun_workspace_v_ball_ref;
  c_FiniteDifferences_nonlin_wo_1 = fun_workspace_a_ball_ref;
  x_0[0] = b_x[0];
  gradf[0] = (rtInf);
  x_0[1] = b_x[1];
  gradf[1] = (rtInf);
  x_0[2] = b_x[2];
  gradf[2] = (rtInf);
  simu_computeFiniteDifferences_j(c_FiniteDifferences_nonlin_work,
    c_FiniteDifferences_nonlin_wo_0, c_FiniteDifferences_nonlin_wo_1, 0.0,
    c_FiniteDifferences_cEq_1, f_temp, x_0, gradf, &b_evalOK, a__3,
    &c_FiniteDifferences);
  c_FiniteDifferences_nonlin_work =
    c_FiniteDifferences.nonlin.workspace.fun.workspace.p_ball_ref;
  c_FiniteDifferences_nonlin_wo_0 =
    c_FiniteDifferences.nonlin.workspace.fun.workspace.v_ball_ref;
  c_FiniteDifferences_nonlin_wo_1 =
    c_FiniteDifferences.nonlin.workspace.fun.workspace.a_ball_ref;
  c_FiniteDifferences_f_1 = c_FiniteDifferences.f_1;
  c_FiniteDifferences_cEq_1 = c_FiniteDifferences.cEq_1;
  funcCount = c_FiniteDifferences.numEvals;
  funcCount++;
  y[0] = gradf[0];
  y[4] = gradf[1];
  y[8] = gradf[2];
  b_x_0 = f_temp;
  b_gamma = 0.01;
  for (b_i = 0; b_i < 3; b_i++) {
    c_exitflag = b_i;
    j = ((c_exitflag + 1) << 2) - 4;
    y[j + 1] = 0.0;
    y[j + 2] = 0.0;
    y[j + 3] = 0.0;
    y[(c_exitflag + (c_exitflag << 2)) + 1] = 0.1;
  }

  B[0] = y[0];
  x_0[0] = B[0];
  gradf[0] = 0.0;
  B[1] = y[4];
  x_0[1] = B[1];
  gradf[1] = 0.0;
  B[2] = y[8];
  x_0[2] = B[2];
  gradf[2] = 0.0;
  c_exitflag = 0;
  for (j = 0; j < 3; j++) {
    aIdx = j + 1;
    c = 0.0;
    for (b_i = aIdx; b_i <= aIdx; b_i++) {
      c += x_0[b_i - 1] * f_temp;
    }

    gradf[c_exitflag] += c;
    c_exitflag++;
  }

  relFactor = 0.0;
  c = gradf[0];
  c = fabs(c);
  b = rtIsNaN(c);
  if (b || (c > 0.0)) {
    relFactor = c;
  }

  c = gradf[1];
  c = fabs(c);
  b = rtIsNaN(c);
  if (b || (c > relFactor)) {
    relFactor = c;
  }

  c = gradf[2];
  c = fabs(c);
  b = rtIsNaN(c);
  if (b || (c > relFactor)) {
    relFactor = c;
  }

  if (!(relFactor >= 1.0)) {
    relFactor = 1.0;
  }

  stepSuccessful = true;
  f_temp = 0.0;
  x_0[0] = b_x[0];
  c = gradf[0];
  c = fabs(c);
  b = rtIsNaN(c);
  if (b || (c > 0.0)) {
    f_temp = c;
  }

  x_0[1] = b_x[1];
  c = gradf[1];
  c = fabs(c);
  b = rtIsNaN(c);
  if (b || (c > f_temp)) {
    f_temp = c;
  }

  x_0[2] = b_x[2];
  c = gradf[2];
  c = fabs(c);
  b = rtIsNaN(c);
  if (b || (c > f_temp)) {
    f_temp = c;
  }

  if ((!(f_temp <= 1.0E-10 * relFactor)) && (funcCount < 600)) {
    a__3[0] = (rtInf);
    a__3[1] = (rtInf);
    a__3[2] = (rtInf);
    if (!(simulink_experiment_debug__norm(a__3) <
          (simulink_experiment_debug__norm(x_0) + 1.4901161193847656E-8) *
          1.0E-6)) {
      do {
        exitg1 = 0;
        c = -b_x_0;
        b_rhs[0] = c;
        b_rhs[1] = 0.0;
        b_rhs[2] = 0.0;
        b_rhs[3] = 0.0;
        memcpy(&y_0[0], &y[0], 12U * sizeof(real_T));
        simulink_e_linearLeastSquares_j(y_0, b_rhs, y, JacCeqTrans);
        f_temp = b_x[0] + JacCeqTrans[0];
        x_0[0] = f_temp;
        xp[0] = f_temp;
        f_temp = b_x[1] + JacCeqTrans[1];
        x_0[1] = f_temp;
        xp[1] = f_temp;
        f_temp = b_x[2] + JacCeqTrans[2];
        x_0[2] = f_temp;
        xp[2] = f_temp;
        f_sym_continuous_idx_0 = fun_workspace_p_ball_ref;
        x2 = fun_workspace_v_ball_ref;
        x3 = x_0[0];
        x4 = x_0[1];
        f_temp = x_0[2];

        /* F_FN */
        /*     F_SYM_CONTINUOUS = F_FN(IN1,U_SYM) */
        /*     This function was generated by the Symbolic Math Toolbox version 9.2. */
        /*     30-Apr-2025 09:26:50 */
        x1 = f_sym_continuous_idx_0;
        c = x4 * x4;
        a = x3;
        a = cos(a);
        a *= a;
        x3 = sin(x3);
        f_sym_continuous_idx_0 = x2;
        x2 = (x1 * 0.7142857142857143 - 0.15196428571428569) * (c * a) *
          0.0035634305945448849 + x3 * 0.41828772872251141;
        x3 = x4;
        x4 = f_temp * 60.0 - x4 * 40.0;
        f_temp = f_sym_continuous_idx_0 - fun_workspace_v_ball_ref;
        f_sym_continuous_idx_0 = x2 - fun_workspace_a_ball_ref;
        f_temp += f_sym_continuous_idx_0;
        f_temp += 0.0 * x3;
        f_temp += 0.0 * x4;
        x1 = f_temp * f_temp;
        b = rtIsInf(f_temp);
        b_evalOK = !b;
        b = rtIsNaN(f_temp);
        b = !b;
        b = (b_evalOK && b);
        funcCount++;
        guard1 = false;
        if ((x1 < resnorm) && b) {
          iter++;
          c = x1 - resnorm;
          c = fabs(c);
          funDiff = c / resnorm;
          resnorm = x1;
          gradf[0] = JacCeqTrans[0];
          gradf[1] = JacCeqTrans[1];
          gradf[2] = JacCeqTrans[2];
          simu_computeFiniteDifferences_j(c_FiniteDifferences_nonlin_work,
            c_FiniteDifferences_nonlin_wo_0, c_FiniteDifferences_nonlin_wo_1,
            c_FiniteDifferences_f_1, c_FiniteDifferences_cEq_1, f_temp, xp,
            gradf, &b_evalOK, a__3, &c_FiniteDifferences);
          c_exitflag = c_FiniteDifferences.numEvals;
          funcCount += c_exitflag;
          y[0] = gradf[0];
          y[4] = gradf[1];
          y[8] = gradf[2];
          b_x_0 = f_temp;
          B[0] = y[0];
          B[1] = y[4];
          B[2] = y[8];
          if (b_evalOK) {
            b_x[0] = xp[0];
            b_x[1] = xp[1];
            b_x[2] = xp[2];
            if (stepSuccessful) {
              b_gamma *= 0.1;
            }

            stepSuccessful = true;
            guard1 = true;
          } else {
            exitg1 = 1;
          }
        } else {
          b_gamma *= 10.0;
          stepSuccessful = false;
          x_0[0] = B[0];
          y[0] = x_0[0];
          x_0[1] = B[1];
          y[4] = x_0[1];
          x_0[2] = B[2];
          y[8] = x_0[2];
          guard1 = true;
        }

        if (guard1) {
          c = b_gamma;
          c = sqrt(c);
          for (b_i = 0; b_i < 3; b_i++) {
            c_exitflag = b_i;
            j = ((c_exitflag + 1) << 2) - 4;
            y[j + 1] = 0.0;
            y[j + 2] = 0.0;
            y[j + 3] = 0.0;
            y[(c_exitflag + (c_exitflag << 2)) + 1] = c;
            x_0[b_i] = B[b_i];
            gradf[c_exitflag] = 0.0;
          }

          c_exitflag = 0;
          for (j = 0; j < 3; j++) {
            aIdx = j + 1;
            c = 0.0;
            for (b_i = aIdx; b_i <= aIdx; b_i++) {
              c += x_0[b_i - 1] * b_x_0;
            }

            gradf[c_exitflag] += c;
            c_exitflag++;
          }

          f_temp = 0.0;
          x_0[0] = JacCeqTrans[0];
          c = gradf[0];
          c = fabs(c);
          b = rtIsNaN(c);
          if (b || (c > 0.0)) {
            f_temp = c;
          }

          x_0[1] = JacCeqTrans[1];
          c = gradf[1];
          c = fabs(c);
          b = rtIsNaN(c);
          if (b || (c > f_temp)) {
            f_temp = c;
          }

          x_0[2] = JacCeqTrans[2];
          c = gradf[2];
          c = fabs(c);
          b = rtIsNaN(c);
          if (b || (c > f_temp)) {
            f_temp = c;
          }

          if ((f_temp <= 1.0E-10 * relFactor) || (funcCount >= 600) || (iter >=
               400) || (simulink_experiment_debug__norm(x_0) <
                        (simulink_experiment_debug__norm(b_x) +
                         1.4901161193847656E-8) * 1.0E-6) || (funDiff <= 1.0E-6))
          {
            exitg1 = 1;
          }
        }
      } while (exitg1 == 0);
    }
  }
}

static void simulink_experiment_debug_t_inv(const real_T x[64], real_T y[64])
{
  real_T A[64];
  real_T s;
  real_T smax;
  real_T x_0;
  int32_T a;
  int32_T b_j;
  int32_T ix;
  int32_T jj;
  int32_T jm1;
  int32_T jp1j;
  int32_T jpiv;
  int32_T jpiv_offset;
  int32_T jy;
  int32_T mmj;
  int32_T pipk;
  int8_T ipiv[8];
  int8_T p[8];
  for (ix = 0; ix < 64; ix++) {
    y[ix] = 0.0;
    A[ix] = x[ix];
  }

  for (ix = 0; ix < 8; ix++) {
    ipiv[ix] = (int8_T)(ix + 1);
  }

  for (b_j = 0; b_j < 7; b_j++) {
    pipk = b_j;
    jm1 = pipk;
    mmj = 7 - pipk;
    jpiv = jm1 * 9;
    a = 1;
    jj = jpiv + 1;
    jp1j = jj + 1;
    jpiv = mmj;
    ix = jj - 1;
    x_0 = A[jj - 1];
    s = fabs(x_0);
    smax = s;
    for (jy = 2; jy <= jpiv + 1; jy++) {
      ix++;
      x_0 = A[ix];
      s = fabs(x_0);
      if (s > smax) {
        a = jy;
        smax = s;
      }
    }

    jpiv_offset = a - 1;
    jpiv = (jj + jpiv_offset) - 1;
    if (A[jpiv] != 0.0) {
      if (jpiv_offset != 0) {
        jpiv = (pipk + jpiv_offset) + 1;
        ipiv[pipk] = (int8_T)jpiv;
        jpiv = jm1;
        ix = jpiv + jpiv_offset;
        for (jy = 0; jy < 8; jy++) {
          smax = A[jpiv];
          A[jpiv] = A[ix];
          A[ix] = smax;
          jpiv += 8;
          ix += 8;
        }
      }

      jpiv = mmj;
      jpiv_offset = (jp1j + jpiv) - 1;
      for (ix = jp1j; ix <= jpiv_offset; ix++) {
        x_0 = A[ix - 1];
        s = A[jj - 1];
        smax = x_0 / s;
        A[ix - 1] = smax;
      }
    }

    jpiv = 7 - pipk;
    jy = jj + 7;
    pipk = jj + 8;
    jj = pipk;
    jpiv_offset = jpiv - 1;
    for (pipk = 0; pipk <= jpiv_offset; pipk++) {
      smax = A[jy];
      if (smax != 0.0) {
        smax = -smax;
        ix = jp1j - 1;
        jpiv = jj;
        jm1 = mmj + jj;
        for (a = jpiv + 1; a <= jm1; a++) {
          A[a - 1] += A[ix] * smax;
          ix++;
        }
      }

      jy += 8;
      jj += 8;
    }
  }

  for (ix = 0; ix < 8; ix++) {
    p[ix] = (int8_T)(ix + 1);
  }

  for (b_j = 0; b_j < 7; b_j++) {
    smax = (real_T)b_j + 1.0;
    jpiv = ipiv[(int32_T)smax - 1] - 1;
    if (jpiv + 1 > (int32_T)smax) {
      pipk = p[jpiv];
      p[jpiv] = p[(int32_T)smax - 1];
      p[(int32_T)smax - 1] = (int8_T)pipk;
    }
  }

  for (b_j = 0; b_j < 8; b_j++) {
    jy = b_j;
    jpiv = p[jy] - 1;
    y[jy + (jpiv << 3)] = 1.0;
    for (pipk = jy + 1; pipk < 9; pipk++) {
      if (y[((jpiv << 3) + pipk) - 1] != 0.0) {
        jpiv_offset = pipk;
        for (ix = jpiv_offset + 1; ix < 9; ix++) {
          y[(ix + (jpiv << 3)) - 1] -= A[(((pipk - 1) << 3) + ix) - 1] * y
            [((jpiv << 3) + pipk) - 1];
        }
      }
    }
  }

  for (b_j = 0; b_j < 8; b_j++) {
    pipk = b_j;
    pipk <<= 3;
    for (jy = 7; jy >= 0; jy--) {
      mmj = jy << 3;
      if (y[jy + pipk] != 0.0) {
        y[jy + pipk] /= A[jy + mmj];
        jpiv_offset = jy;
        jpiv = jpiv_offset - 1;
        for (jpiv_offset = 0; jpiv_offset <= jpiv; jpiv_offset++) {
          ix = jpiv_offset;
          y[ix + pipk] -= y[jy + pipk] * A[ix + mmj];
        }
      }
    }
  }
}

static real_T simulink_experiment_deb_xnrm2_j(int32_T n, const real_T x[32],
  int32_T ix0)
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  int32_T b;
  int32_T kend;
  y = 0.0;
  scale = 3.3121686421112381E-170;
  kend = n;
  kend--;
  kend += ix0;
  for (b = ix0; b <= kend; b++) {
    absxk = x[b - 1];
    absxk = fabs(absxk);
    if (absxk > scale) {
      t = scale / absxk;
      y = y * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      y += t * t;
    }
  }

  return scale * sqrt(y);
}

static void simulink_experiment_de_mldivide(const real_T A[32], const real_T B
  [32], real_T Y[16])
{
  __m128d tmp;
  __m128d tmp_0;
  real_T A_0[32];
  real_T B_0[32];
  real_T x[32];
  real_T x_0[32];
  real_T tau[4];
  real_T vn1[4];
  real_T vn2[4];
  real_T work[4];
  real_T absxk;
  real_T c;
  real_T s;
  real_T scale;
  real_T t;
  int32_T exitg1;
  int32_T ii;
  int32_T im1;
  int32_T ip1;
  int32_T itemp;
  int32_T ix;
  int32_T jy;
  int32_T kend;
  int32_T knt;
  int32_T mm1;
  int32_T mmi;
  int32_T nmi;
  int32_T nmip1;
  int32_T pj;
  int32_T pvtcol;
  int32_T scalarLB;
  int32_T vectorUB;
  int8_T jpvt[4];
  boolean_T exitg2;
  jpvt[0] = 1;
  work[0] = 0.0;
  jpvt[1] = 2;
  work[1] = 0.0;
  jpvt[2] = 3;
  work[2] = 0.0;
  jpvt[3] = 4;
  work[3] = 0.0;
  for (mm1 = 0; mm1 < 32; mm1++) {
    x_0[mm1] = A[mm1];
    B_0[mm1] = B[mm1];
    x[mm1] = x_0[mm1];
  }

  for (ip1 = 0; ip1 < 4; ip1++) {
    nmi = ip1;
    pvtcol = (nmi << 3) + 1;
    s = 0.0;
    scale = 3.3121686421112381E-170;
    kend = pvtcol + 7;
    for (jy = pvtcol; jy <= kend; jy++) {
      absxk = x[jy - 1];
      absxk = fabs(absxk);
      if (absxk > scale) {
        t = scale / absxk;
        s = s * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        s += t * t;
      }
    }

    s = scale * sqrt(s);
    vn1[nmi] = s;
    vn2[nmi] = vn1[nmi];
  }

  for (pj = 0; pj < 4; pj++) {
    mm1 = pj;
    im1 = mm1;
    ip1 = mm1;
    ii = (im1 << 3) + im1;
    nmi = 3 - mm1;
    mmi = 6 - mm1;
    kend = mmi + 1;
    nmip1 = nmi + 1;
    pvtcol = 1;
    if (nmip1 > 1) {
      ix = mm1;
      absxk = vn1[mm1];
      s = absxk;
      scale = s;
      for (jy = 2; jy <= nmip1; jy++) {
        ix++;
        absxk = vn1[ix];
        s = absxk;
        if (s > scale) {
          pvtcol = jy;
          scale = s;
        }
      }
    }

    nmip1 = (im1 + pvtcol) - 1;
    if (nmip1 != mm1) {
      pvtcol = nmip1 << 3;
      ix = im1 << 3;
      for (jy = 0; jy < 8; jy++) {
        absxk = x_0[pvtcol];
        x_0[pvtcol] = x_0[ix];
        x_0[ix] = absxk;
        pvtcol++;
        ix++;
      }

      itemp = jpvt[nmip1];
      jpvt[nmip1] = jpvt[mm1];
      jpvt[mm1] = (int8_T)itemp;
      vn1[nmip1] = vn1[mm1];
      vn2[nmip1] = vn2[mm1];
    }

    scale = x_0[ii];
    pvtcol = ii + 2;
    absxk = 0.0;
    ix = kend - 1;
    s = simulink_experiment_deb_xnrm2_j(ix + 1, x_0, pvtcol);
    if (s != 0.0) {
      t = rt_hypotd_snf(scale, s);
      if (scale >= 0.0) {
        t = -t;
      }

      s = fabs(t);
      if (s < 1.0020841800044864E-292) {
        knt = -1;
        do {
          knt++;
          nmip1 = ix;
          itemp = pvtcol + nmip1;
          scalarLB = ((((itemp - pvtcol) + 1) / 2) << 1) + pvtcol;
          vectorUB = scalarLB - 2;
          for (jy = pvtcol; jy <= vectorUB; jy += 2) {
            tmp_0 = _mm_loadu_pd(&x_0[jy - 1]);
            tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(9.9792015476736E+291));
            _mm_storeu_pd(&x_0[jy - 1], tmp_0);
          }

          for (jy = scalarLB; jy <= itemp; jy++) {
            x_0[jy - 1] *= 9.9792015476736E+291;
          }

          t *= 9.9792015476736E+291;
          scale *= 9.9792015476736E+291;
          s = fabs(t);
        } while ((s < 1.0020841800044864E-292) && (knt + 1 < 20));

        s = simulink_experiment_deb_xnrm2_j(ix + 1, x_0, pvtcol);
        t = rt_hypotd_snf(scale, s);
        if (scale >= 0.0) {
          t = -t;
        }

        absxk = t - scale;
        absxk /= t;
        s = scale - t;
        scale = 1.0 / s;
        itemp = pvtcol + nmip1;
        scalarLB = ((((itemp - pvtcol) + 1) / 2) << 1) + pvtcol;
        vectorUB = scalarLB - 2;
        for (jy = pvtcol; jy <= vectorUB; jy += 2) {
          tmp_0 = _mm_loadu_pd(&x_0[jy - 1]);
          tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(scale));
          _mm_storeu_pd(&x_0[jy - 1], tmp_0);
        }

        for (jy = scalarLB; jy <= itemp; jy++) {
          x_0[jy - 1] *= scale;
        }

        for (jy = 0; jy <= knt; jy++) {
          t *= 1.0020841800044864E-292;
        }

        scale = t;
      } else {
        absxk = t - scale;
        absxk /= t;
        s = scale - t;
        scale = 1.0 / s;
        nmip1 = ix;
        itemp = pvtcol + nmip1;
        scalarLB = ((((itemp - pvtcol) + 1) / 2) << 1) + pvtcol;
        vectorUB = scalarLB - 2;
        for (jy = pvtcol; jy <= vectorUB; jy += 2) {
          tmp_0 = _mm_loadu_pd(&x_0[jy - 1]);
          tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(scale));
          _mm_storeu_pd(&x_0[jy - 1], tmp_0);
        }

        for (jy = scalarLB; jy <= itemp; jy++) {
          x_0[jy - 1] *= scale;
        }

        scale = t;
      }
    }

    tau[mm1] = absxk;
    x_0[ii] = scale;
    if (mm1 + 1 < 4) {
      scale = x_0[ii];
      x_0[ii] = 1.0;
      t = tau[mm1];
      jy = ii + 9;
      if (t != 0.0) {
        nmip1 = kend;
        mm1 = ii + nmip1;
        while ((kend + 1 > 0) && (x_0[mm1] == 0.0)) {
          kend--;
          mm1--;
        }

        exitg2 = false;
        while ((!exitg2) && (nmi > 0)) {
          nmip1 = nmi;
          nmip1 = (nmip1 - 1) << 3;
          itemp = jy + nmip1;
          nmip1 = kend;
          nmip1 += itemp;
          do {
            exitg1 = 0;
            if (itemp <= nmip1) {
              if (x_0[itemp - 1] != 0.0) {
                exitg1 = 1;
              } else {
                itemp++;
              }
            } else {
              nmi--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        kend = -1;
        nmi = 0;
      }

      if (kend + 1 > 0) {
        for (mm1 = 0; mm1 < 32; mm1++) {
          s = x_0[mm1];
          x[mm1] = s;
          A_0[mm1] = s;
        }

        if (nmi != 0) {
          mm1 = kend;
          ix = nmi - 1;
          nmip1 = ix;
          itemp = nmip1 + 1;
          itemp--;
          for (nmip1 = 0; nmip1 <= itemp; nmip1++) {
            knt = nmip1;
            work[knt] = 0.0;
          }

          knt = 0;
          nmip1 = ix << 3;
          pvtcol = jy + nmip1;
          for (nmip1 = jy; nmip1 <= pvtcol; nmip1 += 8) {
            ix = ii;
            c = 0.0;
            itemp = nmip1 + mm1;
            for (scalarLB = nmip1; scalarLB <= itemp; scalarLB++) {
              s = A_0[scalarLB - 1];
              absxk = x[ix];
              s *= absxk;
              c += s;
              ix++;
            }

            work[knt] += c;
            knt++;
          }
        }

        s = -t;
        if (!(s == 0.0)) {
          mm1 = jy;
          jy = 0;
          itemp = nmi - 1;
          for (nmi = 0; nmi <= itemp; nmi++) {
            absxk = work[jy];
            if (absxk != 0.0) {
              absxk *= s;
              ix = ii;
              nmip1 = mm1;
              pvtcol = kend + mm1;
              for (knt = nmip1; knt <= pvtcol; knt++) {
                x_0[knt - 1] += x_0[ix] * absxk;
                ix++;
              }
            }

            jy++;
            mm1 += 8;
          }
        }
      }

      x_0[ii] = scale;
    }

    for (ii = ip1 + 2; ii < 5; ii++) {
      kend = ((ii - 1) << 3) + im1;
      if (vn1[ii - 1] != 0.0) {
        absxk = x_0[kend];
        s = fabs(absxk);
        scale = s / vn1[ii - 1];
        scale = 1.0 - scale * scale;
        if (scale < 0.0) {
          scale = 0.0;
        }

        s = vn1[ii - 1] / vn2[ii - 1];
        s = s * s * scale;
        if (s <= 1.4901161193847656E-8) {
          pvtcol = kend + 2;
          memcpy(&x[0], &x_0[0], sizeof(real_T) << 5U);
          s = 0.0;
          scale = 3.3121686421112381E-170;
          nmip1 = mmi;
          kend = pvtcol + nmip1;
          for (jy = pvtcol; jy <= kend; jy++) {
            absxk = x[jy - 1];
            absxk = fabs(absxk);
            if (absxk > scale) {
              t = scale / absxk;
              s = s * t * t + 1.0;
              scale = absxk;
            } else {
              t = absxk / scale;
              s += t * t;
            }
          }

          s = scale * sqrt(s);
          vn1[ii - 1] = s;
          vn2[ii - 1] = vn1[ii - 1];
        } else {
          scale = sqrt(scale);
          vn1[ii - 1] *= scale;
        }
      }
    }
  }

  mmi = 0;
  absxk = x_0[0];
  s = fabs(absxk);
  t = 1.7763568394002505E-14 * s;
  exitg2 = false;
  while ((!exitg2) && (mmi < 4)) {
    absxk = x_0[(mmi << 3) + mmi];
    s = fabs(absxk);
    if (!(s <= t)) {
      mmi++;
    } else {
      exitg2 = true;
    }
  }

  for (mm1 = 0; mm1 < 16; mm1++) {
    Y[mm1] = 0.0;
  }

  memcpy(&A_0[0], &x_0[0], sizeof(real_T) << 5U);
  for (ip1 = 0; ip1 < 4; ip1++) {
    nmi = ip1;
    t = tau[nmi];
    if (t != 0.0) {
      itemp = nmi;
      for (im1 = 0; im1 < 4; im1++) {
        jy = im1;
        scale = B_0[(jy << 3) + nmi];
        for (mm1 = itemp + 2; mm1 < 9; mm1++) {
          s = A_0[((nmi << 3) + mm1) - 1];
          absxk = B_0[((jy << 3) + mm1) - 1];
          s *= absxk;
          scale += s;
        }

        scale *= t;
        if (scale != 0.0) {
          B_0[nmi + (jy << 3)] -= scale;
          nmip1 = nmi;
          scalarLB = ((((7 - nmip1) / 2) << 1) + nmip1) + 2;
          vectorUB = scalarLB - 2;
          for (pj = nmip1 + 2; pj <= vectorUB; pj += 2) {
            tmp_0 = _mm_loadu_pd(&A_0[((nmi << 3) + pj) - 1]);
            tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(scale));
            tmp = _mm_loadu_pd(&B_0[((jy << 3) + pj) - 1]);
            tmp_0 = _mm_sub_pd(tmp, tmp_0);
            _mm_storeu_pd(&B_0[(pj + (jy << 3)) - 1], tmp_0);
          }

          for (pj = scalarLB; pj < 9; pj++) {
            B_0[(pj + (jy << 3)) - 1] -= A_0[((nmi << 3) + pj) - 1] * scale;
          }
        }
      }
    }
  }

  itemp = mmi - 1;
  for (im1 = 0; im1 < 4; im1++) {
    jy = im1;
    for (pj = 0; pj <= itemp; pj++) {
      mm1 = pj;
      Y[(jpvt[mm1] + (jy << 2)) - 1] = B_0[(jy << 3) + mm1];
    }

    for (nmi = mmi; nmi >= 1; nmi--) {
      pj = jpvt[nmi - 1] - 1;
      absxk = Y[(jy << 2) + pj];
      s = x_0[(((nmi - 1) << 3) + nmi) - 1];
      s = absxk / s;
      Y[pj + (jy << 2)] = s;
      nmip1 = nmi;
      pvtcol = nmip1 - 2;
      for (ip1 = 0; ip1 <= pvtcol; ip1++) {
        mm1 = ip1;
        Y[(jpvt[mm1] + (jy << 2)) - 1] -= x_0[((nmi - 1) << 3) + mm1] * Y[(jy <<
          2) + pj];
      }
    }
  }
}

static void studentControllerInterface_step(studentControllerInterface_si_T *obj,
  real_T t, real_T p_ball, real_T theta, real_T *V_servo, real_T x_hat[16])
{
  __m128d tmp_0;
  __m128d tmp_1;
  studentControllerInterface_si_T *obj_0;
  studentControllerInterface_si_T *varargin_3;
  real_T Z[64];
  real_T tmp[64];
  real_T M[32];
  real_T N[32];
  real_T N_0[32];
  real_T A[16];
  real_T G[16];
  real_T W11[16];
  real_T W21[16];
  real_T W22[16];
  real_T K_gain[8];
  real_T a[8];
  real_T b_y[8];
  real_T mt2[5];
  real_T x_pred[4];
  real_T x_prev[4];
  real_T mt1[3];
  real_T obj_1[3];
  real_T a22;
  real_T p_ball_0;
  real_T t10;
  real_T t11;
  real_T t12;
  real_T t13;
  real_T t14;
  real_T t15;
  real_T t16;
  real_T t17;
  real_T t18;
  real_T t19;
  real_T t2;
  real_T t21;
  real_T t22;
  real_T t23;
  real_T t24;
  real_T t25;
  real_T t27;
  real_T t28;
  real_T t30;
  real_T t32;
  real_T t33;
  real_T t6;
  real_T t7;
  real_T t8;
  real_T t9;
  real_T t_prev;
  real_T theta_0;
  real_T u_prev;
  real_T x4;
  real_T y_idx_0;
  real_T y_idx_1;
  int32_T r1;
  int32_T ret;
  char_T a_1[6];
  char_T b_0[6];
  char_T a_0[3];
  char_T b[3];
  boolean_T equal;
  static const char_T tmp_2[6] = { 'T', 'V', '-', 'L', 'Q', 'R' };

  static const int8_T tmp_3[8] = { 1, 0, 0, 0, 0, 1, 0, 0 };

  /*  This is the main function called every iteration. You have to implement */
  /*  the controller in this function, bu you are not allowed to */
  /*  change the signature of this function.  */
  /*  Input arguments: */
  /*    t: current time */
  /*    p_ball: position of the ball provided by the ball position sensor (m) */
  /*  */
  /*    theta: servo motor angle provided by the encoder of the motor (rad) */
  /*  Output: */
  /*    V_servo: voltage to the servo input.    */
  t_prev = obj->t_prev;
  obj->t_prev = t;
  obj->dt = t - t_prev;

  /*             %% Sample Controller: Simple Proportional Controller */
  /*  Extract reference trajectory at the current timestep. */
  *V_servo = 0.0;
  for (r1 = 0; r1 < 16; r1++) {
    x_hat[r1] = 0.0;
  }

  /*  Compute state estimate */
  a_0[0] = obj->observer[0];
  a_0[1] = obj->observer[1];
  a_0[2] = obj->observer[2];
  b[0] = 'E';
  b[1] = 'L';
  b[2] = 'O';
  ret = memcmp(&a_0[0], &b[0], 3);
  equal = (ret == 0);
  if (equal) {
    obj_0 = obj;
    x_pred[0] = obj->x_hat[0];
    x_pred[1] = obj->x_hat[1];
    x_pred[2] = obj->x_hat[2];
    x_pred[3] = obj->x_hat[3];
    u_prev = obj->u;
    y_idx_0 = p_ball;
    y_idx_1 = theta;

    /*  Co = ctrb(A_eval', C_eval'); */
    /*  rank_Co = rank(Co); */
    /*  disp(['Rank of Controllability Matrix: ', num2str(rank_Co)]); */
    /*  if rank_Co < 4 */
    /*      poles = [-4, -4.5, -5, -5.5]; */
    /*  else */
    /*      poles = [-1, -1.5, -2, -2.5]; */
    /*  end */
    /*              L = place_fn(A_eval', C_eval', poles)'; */
    if (x_pred[3] == 0.0) {
      x_pred[3] = 2.2204460492503131E-16;
    }

    /* LUENBERGER_GAINS_FUNC */
    /*     OUT1 = LUENBERGER_GAINS_FUNC(IN1,U_SYM) */
    /*     This function was generated by the Symbolic Math Toolbox version 9.2. */
    /*     30-Apr-2025 09:26:50 */
    t_prev = x_pred[0];
    a22 = x_pred[2];
    x4 = x_pred[3];
    t2 = a22;
    t2 = cos(t2);
    a22 = sin(a22);
    theta_0 = t_prev * t_prev;
    p_ball_0 = x4 * x4;
    t6 = rt_powd_snf(x4, 3.0);
    t8 = rt_powd_snf(x4, 5.0);
    t7 = p_ball_0 * p_ball_0;
    t9 = rt_powd_snf(p_ball_0, 3.0);
    t10 = t2 * t2;
    t11 = rt_powd_snf(t2, 3.0);
    t13 = rt_powd_snf(t2, 5.0);
    t15 = a22 * a22;
    t12 = t10 * t10;
    t14 = rt_powd_snf(t10, 3.0);
    t17 = t11 * x4 * 5.9032181842099944E+49;
    t18 = t10 * 5.6998592900689022E+50;
    t22 = a22 * t6 * t11 * 1.528458123128785E+47;
    t23 = t11 * t_prev * x4 * 2.7747206506274942E+50;
    t24 = t7 * t10 * t15 * 3.8211453078219642E+45;
    t25 = a22 * p_ball_0 * t10 * 2.9516090921049968E+48;
    t27 = t7 * t10 * t15 * t_prev * 3.592146000302668E+46;
    t28 = a22 * t6 * t11 * t_prev * 1.4368584001210671E+48;
    t30 = a22 * p_ball_0 * t10 * t_prev * 1.3873603253137469E+49;
    t32 = a22 * theta_0 * t6 * t11 * 3.3768705055724263E+48;
    t33 = theta_0 * t7 * t10 * t15 * 8.4421762639310657E+46;
    t16 = p_ball_0 * t12 * 3.8211453078219642E+45;
    t19 = p_ball_0 * t12 * t_prev * 3.592146000302668E+46;
    t21 = theta_0 * p_ball_0 * t12 * 8.4421762639310657E+46;
    t23 = -t23;
    t27 = -t27;
    t28 = -t28;
    t30 = -t30;
    t19 = -t19;
    t17 = ((((((((((((t16 + t17) + t18) + t19) + t21) + t22) + t24) + t25) + t23)
              + t27) + t28) + t32) + t33) + t30;
    t17 = 1.0 / t17;
    t18 = (p_ball_0 * t12 * 2.1634396204131021E+56 + t10 * 1.032992910381229E+61)
      + t6 * t13 * 1.5234080974899679E+55;
    t22 = (t7 * t14 * 9.8610004271776722E+50 + t11 * x4 * 1.069847907183821E+60)
      + a22 * p_ball_0 * t10 * 5.3492395359191073E+58;
    t24 = (a22 * t6 * t11 * 2.7700445303230439E+57 + theta_0 * p_ball_0 * t12 *
           1.529986581251914E+57) + a22 * t7 * t12 * 7.61704048744984E+53;
    t25 = (a22 * t8 * t13 * 3.9444001708710687E+52 + theta_0 * t7 * t14 *
           2.178621775375107E+52) + t7 * t10 * t15 * 6.92511132580761E+55;
    t32 = (t9 * t12 * t15 * 9.8610004271776722E+50 - p_ball_0 * t12 * t_prev *
           6.5100929032268949E+56) - t6 * t13 * t_prev * 7.1605550998353374E+55;
    t14 = (t7 * t14 * t_prev * -9.2700356542210782E+51 - t11 * t_prev * x4 *
           5.0286623134374681E+60) + a22 * theta_0 * t6 * t11 *
      6.1199463250076569E+58;
    t33 = (a22 * theta_0 * t8 * t13 * 8.7144871015004267E+53 + theta_0 * t7 *
           t10 * t15 * 1.529986581251914E+57) + theta_0 * t9 * t12 * t15 *
      2.178621775375107E+52;
    t16 = (a22 * p_ball_0 * t10 * t_prev * -2.5143311567187341E+59 - a22 * t6 *
           t11 * t_prev * 2.6040371612907581E+58) - a22 * t7 * t12 * t_prev *
      3.5802775499176687E+54;
    t8 = (a22 * t8 * t13 * t_prev * -3.7080142616884311E+53 - t7 * t10 * t15 *
          t_prev * 6.5100929032268949E+56) - t9 * t12 * t15 * t_prev *
      9.2700356542210782E+51;
    t9 = ((p_ball_0 * t12 * 1.023672660835295E+34 - t11 * x4 *
           6.4506649997227849E+36) - a22 * t6 * t11 * 1.6702027624154809E+34) +
      theta_0 * p_ball_0 * t12 * 2.261632139884469E+35;
    t13 = ((p_ball_0 * t12 * t_prev * -9.6232447552084133E+34 + t11 * t_prev *
            x4 * 3.0320399528661739E+37) - a22 * theta_0 * t6 * t11 *
           3.6900313861272912E+35) + a22 * t6 * t11 * t_prev *
      1.570108354797162E+35;
    t21 = ((p_ball_0 * t12 * 2.9323761254803291E-8 + t10 * 0.004374115599995931)
           + t11 * x4 * 0.000453017477023024) + a22 * p_ball_0 * t10 *
      2.26508738511512E-5;
    t12 = (((a22 * t6 * t11 * 1.1729504501921319E-6 + theta_0 * p_ball_0 * t12 *
             6.4785906133359738E-7) + t7 * t10 * t15 * 2.9323761254803291E-8) -
           p_ball_0 * t12 * t_prev * 2.7566403059744572E-7) - t11 * t_prev * x4 *
      0.0021293418426464122;
    theta_0 = (((a22 * theta_0 * t6 * t11 * 2.59143624533439E-5 + theta_0 * t7 *
                 t10 * t15 * 6.4785906133359738E-7) - a22 * p_ball_0 * t10 *
                t_prev * 0.00010646709213232059) - a22 * t6 * t11 * t_prev *
               1.102656122389783E-5) - t7 * t10 * t15 * t_prev *
      2.7566403059744572E-7;
    mt1[1] = ((((((((t18 + t22) + t24) + t25) + t32) + t14) + t33) + t16) + t8) *
      t17 / 1.0138814E+8;
    mt1[2] = (t9 + t13) * t17 * 1.125899906842624E+15 / (t10 * x4 *
      0.0010830283699848921 - t10 * t_prev * x4 * 0.0050906151350641211);
    mt2[0] = ((((t10 * x4 * 0.13324633414470369 + t2 * 2090.0399940196412) + t2 *
                a22 * p_ball_0 * 5.4115204738123222) - t10 * t_prev * x4 *
               0.62630474333585762) - t2 * a22 * p_ball_0 * t_prev *
              25.436053930962728) / ((t21 + t12) + theta_0);
    mt2[1] = 0.0;
    mt2[2] = 0.0;
    mt2[3] = 0.0;
    mt2[4] = 1.0;
    a[0] = 9.0;
    a[1] = mt1[1];
    a[2] = mt1[2];
    for (r1 = 0; r1 < 5; r1++) {
      a[r1 + 3] = mt2[r1];
    }

    /*  disp(L) */
    /* H_FN */
    /*     H_SYM_CONTINUOUS = H_FN(IN1,U_SYM) */
    /*     This function was generated by the Symbolic Math Toolbox version 9.2. */
    /*     30-Apr-2025 09:26:50 */
    a22 = x_pred[2];
    theta_0 = t_prev;
    t2 = a22;

    /*  disp(hat_y) */
    p_ball_0 = y_idx_0;
    p_ball_0 -= theta_0;
    y_idx_0 = p_ball_0;
    p_ball_0 = y_idx_1;
    p_ball_0 -= t2;
    y_idx_1 = p_ball_0;

    /* F_FN */
    /*     F_SYM_CONTINUOUS = F_FN(IN1,U_SYM) */
    /*     This function was generated by the Symbolic Math Toolbox version 9.2. */
    /*     30-Apr-2025 09:26:50 */
    theta_0 = x_pred[1];
    p_ball_0 = x4 * x4;
    t2 = cos(t2);
    t2 *= t2;
    a22 = sin(a22);
    x_prev[0] = theta_0;
    x_prev[1] = (t_prev * 0.7142857142857143 - 0.15196428571428569) * (p_ball_0 *
      t2) * 0.0035634305945448849 + a22 * 0.41828772872251141;
    x_prev[2] = x4;
    x_prev[3] = u_prev * 60.0 - x4 * 40.0;
    for (r1 = 0; r1 <= 2; r1 += 2) {
      tmp_0 = _mm_loadu_pd(&a[r1]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(y_idx_0));
      tmp_0 = _mm_add_pd(tmp_0, _mm_set1_pd(0.0));
      tmp_1 = _mm_loadu_pd(&a[r1 + 4]);
      tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(y_idx_1));
      tmp_0 = _mm_add_pd(tmp_1, tmp_0);
      tmp_1 = _mm_loadu_pd(&x_prev[r1]);
      tmp_0 = _mm_add_pd(tmp_1, tmp_0);
      _mm_storeu_pd(&x_prev[r1], tmp_0);
    }

    t2 = obj_0->dt;
    y_idx_0 = x_pred[0];
    t_prev = x_prev[0];
    t_prev *= t2;
    y_idx_0 += t_prev;
    x_pred[0] = y_idx_0;
    y_idx_0 = x_pred[1];
    t_prev = x_prev[1];
    t_prev *= t2;
    y_idx_0 += t_prev;
    x_pred[1] = y_idx_0;
    y_idx_0 = x_pred[2];
    t_prev = x_prev[2];
    t_prev *= t2;
    y_idx_0 += t_prev;
    x_pred[2] = y_idx_0;
    y_idx_0 = x_pred[3];
    t_prev = x_prev[3];
    t_prev *= t2;
    y_idx_0 += t_prev;
    x_pred[3] = y_idx_0;
    obj->x_hat[0] = x_pred[0];
    obj->x_hat[1] = x_pred[1];
    obj->x_hat[2] = x_pred[2];
    obj->x_hat[3] = x_pred[3];
  } else {
    a_0[0] = obj->observer[0];
    a_0[1] = obj->observer[1];
    a_0[2] = obj->observer[2];
    b[0] = 'E';
    b[1] = 'K';
    b[2] = 'F';
    ret = memcmp(&a_0[0], &b[0], 3);
    equal = (ret == 0);
    if (equal) {
      obj_0 = obj;
      y_idx_0 = p_ball;
      y_idx_1 = theta;

      /*  function setupMHE(obj) */
      /*      import casadi.* */
      /*   */
      /*      g = 9.81; */
      /*      r_arm = 0.0254; */
      /*      L = 0.4255; */
      /*   */
      /*      a = 5 * g * r_arm / (7 * L); */
      /*      b = (5 * L / 14) * (r_arm / L)^2; */
      /*      c = (5 / 7) * (r_arm / L)^2; */
      /*   */
      /*      K = 1.5; */
      /*      tau = 0.025; */
      /*   */
      /*      vars = MX.sym('x', 6); */
      /*      x = MX.sym('x', 4); */
      /*      u = MX.sym('u'); */
      /*      dt = MX.sym('dt'); */
      /*      f = Function('f', {x, u}, {[ */
      /*          x(2) */
      /*          u * sin(x(3)) - b * x(4)^2 * cos(x(3))^2 + c * x(1) * x(4)^2 * cos(x(3))^2 */
      /*          x(4) */
      /*          (- x(4) + K * u) / tau */
      /*      ]}, {'x', 'u'}, {'xdot'}); */
      /*   */
      /*      % rk4 integration for discretization */
      /*      k1 = dt*f(x, u); */
      /*      k2 = dt*f(x+dt*k1/2, u); */
      /*      k3 = dt*f(x+dt*k2/2, u); */
      /*      k4 = dt*f(x+dt*k3, u); */
      /*      xf = x + (k1 + 2*k2 + 2*k3 + k4)/6; */
      /*   */
      /*      xf = x + dt*f(x+(dt/2)*f(x, u), u); % also works if we need it */
      /*      %a bit faster. doesn't really hurt performance tbh */
      /*      F = Function('F', {x, u, dt}, {xf}).map(obj.N_MHE-1); */
      /*   */
      /*      opti = Opti(); */
      /*      % X0 = opti.variable(4); */
      /*      % U0 = opti.variable(1); */
      /*      vars = opti.variable(7, obj.N_MHE); */
      /*      X = vars(1:4, :); */
      /*      W = vars(5:7, 1:(end-1)); */
      /*      %X = opti.variable(4, obj.N_MHE); */
      /*      %W = opti.variable(obj.N_MHE-1); */
      /*      U = opti.parameter(1, obj.N_MHE-1); */
      /*      Y = opti.parameter(2, obj.N_MHE); */
      /*      DT = opti.parameter(obj.N_MHE-1); */
      /*      X_prior = opti.parameter(4); */
      /*      P_est = opti.parameter(4, 4); */
      /*      disp(size(W(2, :))); */
      /*      disp(obj.N_MHE); */
      /*      perturbation = [zeros(1, obj.N_MHE-1); W(2, :); zeros(1, obj.N_MHE-1); W(3, :)]; */
      /*      disp(perturbation); */
      /*      dynamics_gap = X(:, 2:end) - F(X(:,1:end-1), U + W(1, :), DT) + perturbation; */
      /*      observation_gap = X([1,3], :) - Y; */
      /*      %dummy_param = opti.parameter(1); */
      /*      %opti.set_value(dummy_param, 0); */
      /*   */
      /*      % disp(size(dynamics_gap)); */
      /*      % disp(size(observation_gap)); */
      /*      theta_max = deg2rad(60); % physical limit */
      /*      x_max = 0.20; % physical limit */
      /*   */
      /*      cost = 0.0; */
      /*      %opti.subject_to( X0==[0;0;0;0]) */
      /*      %opti.subject_to(X(:, 1) - X0*dummy_param + U0*dummy_param == 0) */
      /*      for i=1:(obj.N_MHE) */
      /*          if(i<obj.N_MHE) */
      /*              cost = cost + bilin(obj.Q_est, W(:, i))*DT(i) + bilin(obj.R_est, observation_gap(:, i))*DT(i); */
      /*              opti.subject_to(dynamics_gap(:, i) == [0; 0; 0; 0]); */
      /*          end */
      /*          opti.subject_to(-x_max <= X(1, i)); */
      /*          opti.subject_to(          X(1, i) <= x_max); */
      /*          opti.subject_to(-theta_max <= X(3, i)); */
      /*          opti.subject_to(              X(3, i) <= theta_max); */
      /*      end */
      /*      cost = cost + bilin(obj.R_est, observation_gap(:, obj.N_MHE))*0.01; */
      /*      cost = cost + bilin(P_est, X(:, 1) - X_prior)*DT(1); */
      /*      cost = cost + sum(vars(5:7, end) .^ 2); */
      /*      opti.minimize(cost); */
      /*      opts = struct; */
      /*      opts.expand = true; */
      /*      opts.ipopt.linear_solver = 'mumps'; % default; comes preinstalled. Small problem so mumps is good */
      /*      opts.ipopt.print_level = 0; */
      /*      opts.print_time = 0; */
      /*      opts.ipopt.max_wall_time = 0.010; % 10ms is super safe - typ. is 2-3ms */
      /*      %opts = struct; */
      /*      %opts.structure_detection = 'auto'; */
      /*      %opts.debug = true; */
      /*      opti.solver('ipopt', opts); */
      /*   */
      /*      obj.opti = opti; */
      /*      obj.X_opt = X; */
      /*      obj.U_opt = U; */
      /*      obj.Y_opt = Y; */
      /*      obj.DT_opt = DT; */
      /*      obj.W_opt = W; */
      /*      obj.X_prior = X_prior; */
      /*      obj.P_est = P_est; */
      /*      obj.X_prior_num = [0,0,0,0]; */
      /*      obj.P_est_num = zeros(4, 4); */
      /*   */
      /*  end */
      /*  function xhat = MovingWindowEstimator(obj, dt, y) */
      /*      obj.history = [obj.history(:, 2:end), [0; 0; y(2); y(3)]]; % zero is unused. */
      /*      obj.history(1, end-1) = dt; % technically dt is dt_prev */
      /*      obj.history(2, end-1) = y(1); % y(1) is u_prev */
      /*      %disp(obj.history); */
      /*      obj.opti.set_value(obj.X_prior, obj.X_prior_num); */
      /*      obj.opti.set_value(obj.P_est, obj.P_est_num); */
      /*   */
      /*      obj.opti.set_value(obj.U_opt, clip(obj.history(2, 1:end-1), -10, 10)); */
      /*      obj.opti.set_value(obj.Y_opt, obj.history(3:4, :)); */
      /*      obj.opti.set_value(obj.DT_opt, obj.history(1, 1:end-1)); */
      /*      sol = obj.opti.solve_limited(); % solve_limited makes it not error if it hits time or iter limits */
      /*      Xhat = sol.value(obj.X_opt); */
      /*      What = sol.value(obj.W_opt); */
      /*      obj.opti.set_initial(obj.X_opt, Xhat); % set up warmstarting */
      /*      obj.opti.set_initial(obj.W_opt, What); */
      /*      obj.X_prior_num = Xhat(:, 2); */
      /*      xprior_cell = num2cell(obj.X_prior_num); */
      /*      Uhat = What(1, 2) + clip(obj.history(2, 2), -10, 10); */
      /*      G = [obj.B_fn(xprior_cell{:}, Uhat), [0;1;0;0], [0;0;0;1]]*obj.history(1, 1); */
      /*      A = obj.A_fn(xprior_cell{:}, Uhat)*obj.history(1, 1) + eye(4); */
      /*      P = obj.P_est_num; */
      /*      C = [1 0 0 0; 0 0 1 0]; */
      /*      obj.P_est_num = G*obj.Q_est*G' + A*P*A' - A*P*C'*inv(obj.R_est + C*P*C')*C*P*A'; */
      /*      % disp(obj.P_est_num) */
      /*      xhat = Xhat(:, end); */
      /*      % disp(Xhat); */
      /*  end */
      p_ball_0 = y_idx_0;
      theta_0 = y_idx_1;
      u_prev = obj_0->u;
      x_prev[0] = obj_0->xhat_prev[0];
      x_prev[1] = obj_0->xhat_prev[1];
      x_prev[2] = obj_0->xhat_prev[2];
      x_prev[3] = obj_0->xhat_prev[3];
      varargin_3 = obj_0;
      t2 = varargin_3->rg_val / varargin_3->L_val * x_prev[3];
      a22 = t2 * t2;
      t_prev = x_prev[2];
      t_prev = cos(t_prev);
      x4 = t_prev * t_prev;
      t2 = varargin_3->dt;
      t_prev = x_prev[2];
      t_prev = sin(t_prev);
      x_pred[0] = x_prev[1];
      x_pred[1] = 5.0 * varargin_3->g_val * varargin_3->rg_val / (7.0 *
        varargin_3->L_val) * t_prev - (varargin_3->L_val / 2.0 - x_prev[0]) *
        0.7142857142857143 * a22 * x4;
      x_pred[2] = x_prev[3];
      x_pred[3] = varargin_3->K_val / varargin_3->tau_val * u_prev + -x_prev[3] /
        varargin_3->tau_val;
      y_idx_0 = x_pred[0];
      y_idx_0 *= t2;
      y_idx_0 += x_prev[0];
      x_pred[0] = y_idx_0;
      y_idx_0 = x_pred[1];
      y_idx_0 *= t2;
      y_idx_0 += x_prev[1];
      x_pred[1] = y_idx_0;
      y_idx_0 = x_pred[2];
      y_idx_0 *= t2;
      y_idx_0 += x_prev[2];
      x_pred[2] = y_idx_0;
      y_idx_0 = x_pred[3];
      y_idx_0 *= t2;
      y_idx_0 += x_prev[3];
      x_pred[3] = y_idx_0;
      t6 = obj_0->dt;
      t10 = obj_0->g_val;
      t7 = obj_0->L_val;
      t11 = obj_0->rg_val;
      t15 = obj_0->tau_val;
      for (r1 = 0; r1 < 16; r1++) {
        A[r1] = 0.0;
      }

      A[0] = 1.0;
      A[5] = 1.0;
      A[10] = 1.0;
      A[4] = t6;
      t2 = t11 / t7;
      a22 = t2 * t2;
      t2 = x_prev[3];
      x4 = t2 * t2;
      t_prev = x_prev[2];
      t_prev = cos(t_prev);
      t_prev *= t_prev;
      A[1] = a22 * x4 * t_prev * (0.7142857142857143 * t6);
      t2 = t11 / t7;
      a22 = t2 * t2;
      t2 = x_prev[3];
      x4 = t2 * t2;
      t_prev = x_prev[2];
      t_prev = cos(t_prev);
      u_prev = x_prev[2];
      u_prev = cos(u_prev);
      t2 = x_prev[2];
      t2 = sin(t2);
      A[9] = (t7 / 2.0 - x_prev[0]) * a22 * x4 * u_prev * t2 *
        (1.4285714285714286 * t6) + 5.0 * t10 * t11 / (7.0 * t7) * t6 * t_prev;
      t2 = t11 / t7;
      a22 = t2 * t2;
      t_prev = x_prev[2];
      t_prev = cos(t_prev);
      x4 = t_prev * t_prev;
      A[13] = (t7 / 2.0 - x_prev[0]) * a22 * x_prev[3] * x4 *
        (-1.4285714285714286 * t6);
      A[14] = t6;
      A[15] = 1.0 - t6 / t15;
      for (r1 = 0; r1 < 16; r1++) {
        W11[r1] = obj_0->P[r1];
      }

      for (r1 = 0; r1 < 4; r1++) {
        for (ret = 0; ret < 4; ret++) {
          W21[r1 + (ret << 2)] = 0.0;
          t_prev = W21[(ret << 2) + r1];
          t_prev += W11[ret << 2] * A[r1];
          W21[r1 + (ret << 2)] = t_prev;
          t_prev = W21[(ret << 2) + r1];
          t_prev += W11[(ret << 2) + 1] * A[r1 + 4];
          W21[r1 + (ret << 2)] = t_prev;
          t_prev = W21[(ret << 2) + r1];
          t_prev += W11[(ret << 2) + 2] * A[r1 + 8];
          W21[r1 + (ret << 2)] = t_prev;
          t_prev = W21[(ret << 2) + r1];
          t_prev += W11[(ret << 2) + 3] * A[r1 + 12];
          W21[r1 + (ret << 2)] = t_prev;
        }

        for (ret = 0; ret < 4; ret++) {
          G[r1 + (ret << 2)] = 0.0;
          y_idx_1 = G[(ret << 2) + r1];
          y_idx_1 += W21[r1] * A[ret];
          G[r1 + (ret << 2)] = y_idx_1;
          y_idx_1 = G[(ret << 2) + r1];
          y_idx_1 += W21[r1 + 4] * A[ret + 4];
          G[r1 + (ret << 2)] = y_idx_1;
          y_idx_1 = G[(ret << 2) + r1];
          y_idx_1 += W21[r1 + 8] * A[ret + 8];
          G[r1 + (ret << 2)] = y_idx_1;
          y_idx_1 = G[(ret << 2) + r1];
          y_idx_1 += W21[r1 + 12] * A[ret + 12];
          G[r1 + (ret << 2)] = y_idx_1;
        }
      }

      for (r1 = 0; r1 < 16; r1++) {
        obj_0->P[r1] = G[r1] + obj_0->Q[r1];
      }

      y_idx_0 = p_ball_0;
      y_idx_1 = theta_0;
      theta_0 = x_pred[0];
      t2 = x_pred[2];
      for (r1 = 0; r1 < 8; r1++) {
        a[r1] = tmp_3[r1];
      }

      for (r1 = 0; r1 < 16; r1++) {
        W11[r1] = obj_0->P[r1];
      }

      for (r1 = 0; r1 < 2; r1++) {
        for (ret = 0; ret < 4; ret++) {
          b_y[r1 + (ret << 1)] = 0.0;
          t_prev = b_y[(ret << 1) + r1];
          t_prev += W11[ret << 2] * a[r1];
          b_y[r1 + (ret << 1)] = t_prev;
          t_prev = b_y[(ret << 1) + r1];
          t_prev += W11[(ret << 2) + 1] * a[r1 + 2];
          b_y[r1 + (ret << 1)] = t_prev;
          t_prev = b_y[(ret << 1) + r1];
          t_prev += W11[(ret << 2) + 2] * a[r1 + 4];
          b_y[r1 + (ret << 1)] = t_prev;
          t_prev = b_y[(ret << 1) + r1];
          t_prev += W11[(ret << 2) + 3] * a[r1 + 6];
          b_y[r1 + (ret << 1)] = t_prev;
        }

        for (ret = 0; ret < 2; ret++) {
          x_prev[r1 + (ret << 1)] = 0.0;
          x_prev[r1 + (ret << 1)] += b_y[r1] * a[ret];
          x_prev[r1 + (ret << 1)] += b_y[r1 + 2] * a[ret + 2];
          x_prev[r1 + (ret << 1)] += b_y[r1 + 4] * a[ret + 4];
          x_prev[r1 + (ret << 1)] += b_y[r1 + 6] * a[ret + 6];
        }
      }

      t_prev = x_prev[0];
      t_prev += obj_0->R[0];
      x_prev[0] = t_prev;
      t_prev = x_prev[1];
      t_prev += obj_0->R[1];
      x_prev[1] = t_prev;
      t_prev = x_prev[2];
      t_prev += obj_0->R[2];
      x_prev[2] = t_prev;
      t_prev = x_prev[3];
      t_prev += obj_0->R[3];
      x_prev[3] = t_prev;
      for (r1 = 0; r1 < 16; r1++) {
        W21[r1] = obj_0->P[r1];
      }

      for (r1 = 0; r1 < 4; r1++) {
        for (ret = 0; ret < 2; ret++) {
          b_y[r1 + (ret << 2)] = 0.0;
          b_y[r1 + (ret << 2)] += W21[r1] * a[ret];
          b_y[r1 + (ret << 2)] += W21[r1 + 4] * a[ret + 2];
          b_y[r1 + (ret << 2)] += W21[r1 + 8] * a[ret + 4];
          b_y[r1 + (ret << 2)] += W21[r1 + 12] * a[ret + 6];
        }
      }

      t_prev = x_prev[1];
      p_ball_0 = fabs(t_prev);
      a22 = p_ball_0;
      t_prev = x_prev[0];
      p_ball_0 = fabs(t_prev);
      t_prev = p_ball_0;
      if (a22 > t_prev) {
        r1 = 1;
        ret = 0;
      } else {
        r1 = 0;
        ret = 1;
      }

      t_prev = x_prev[ret] / x_prev[r1];
      a22 = x_prev[ret + 2] - x_prev[r1 + 2] * t_prev;
      K_gain[r1 << 2] = b_y[0] / x_prev[r1];
      K_gain[ret << 2] = (b_y[4] - K_gain[r1 << 2] * x_prev[r1 + 2]) / a22;
      K_gain[r1 << 2] -= K_gain[ret << 2] * t_prev;
      K_gain[(r1 << 2) + 1] = b_y[1] / x_prev[r1];
      K_gain[(ret << 2) + 1] = (b_y[5] - K_gain[(r1 << 2) + 1] * x_prev[r1 + 2])
        / a22;
      K_gain[(r1 << 2) + 1] -= K_gain[(ret << 2) + 1] * t_prev;
      K_gain[(r1 << 2) + 2] = b_y[2] / x_prev[r1];
      K_gain[(ret << 2) + 2] = (b_y[6] - K_gain[(r1 << 2) + 2] * x_prev[r1 + 2])
        / a22;
      K_gain[(r1 << 2) + 2] -= K_gain[(ret << 2) + 2] * t_prev;
      K_gain[(r1 << 2) + 3] = b_y[3] / x_prev[r1];
      K_gain[(ret << 2) + 3] = (b_y[7] - K_gain[(r1 << 2) + 3] * x_prev[r1 + 2])
        / a22;
      K_gain[(r1 << 2) + 3] -= K_gain[(ret << 2) + 3] * t_prev;
      p_ball_0 = y_idx_0;
      p_ball_0 -= theta_0;
      y_idx_0 = p_ball_0;
      p_ball_0 = y_idx_1;
      p_ball_0 -= t2;
      y_idx_1 = p_ball_0;
      for (r1 = 0; r1 <= 2; r1 += 2) {
        tmp_0 = _mm_loadu_pd(&K_gain[r1]);
        tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(y_idx_0));
        tmp_0 = _mm_add_pd(tmp_0, _mm_set1_pd(0.0));
        tmp_1 = _mm_loadu_pd(&K_gain[r1 + 4]);
        tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(y_idx_1));
        tmp_0 = _mm_add_pd(tmp_1, tmp_0);
        _mm_storeu_pd(&x_prev[r1], tmp_0);
      }

      obj_0->xh[0] = x_pred[0] + x_prev[0];
      obj_0->xh[1] = x_pred[1] + x_prev[1];
      obj_0->xh[2] = x_pred[2] + x_prev[2];
      obj_0->xh[3] = x_pred[3] + x_prev[3];
      for (r1 = 0; r1 < 4; r1++) {
        for (ret = 0; ret <= 2; ret += 2) {
          _mm_storeu_pd(&G[ret + (r1 << 2)], _mm_set1_pd(0.0));
          tmp_0 = _mm_loadu_pd(&K_gain[ret]);
          tmp_0 = _mm_mul_pd(_mm_set1_pd(a[r1 << 1]), tmp_0);
          tmp_1 = _mm_loadu_pd(&G[(r1 << 2) + ret]);
          tmp_0 = _mm_add_pd(tmp_0, tmp_1);
          _mm_storeu_pd(&G[ret + (r1 << 2)], tmp_0);
          tmp_0 = _mm_loadu_pd(&K_gain[ret + 4]);
          tmp_0 = _mm_mul_pd(_mm_set1_pd(a[(r1 << 1) + 1]), tmp_0);
          tmp_1 = _mm_loadu_pd(&G[(r1 << 2) + ret]);
          tmp_0 = _mm_add_pd(tmp_0, tmp_1);
          _mm_storeu_pd(&G[ret + (r1 << 2)], tmp_0);
        }
      }

      for (r1 = 0; r1 < 16; r1++) {
        A[r1] = 0.0;
      }

      A[0] = 1.0;
      A[5] = 1.0;
      A[10] = 1.0;
      A[15] = 1.0;
      for (r1 = 0; r1 <= 14; r1 += 2) {
        tmp_0 = _mm_loadu_pd(&A[r1]);
        tmp_1 = _mm_loadu_pd(&G[r1]);
        tmp_0 = _mm_sub_pd(tmp_0, tmp_1);
        _mm_storeu_pd(&A[r1], tmp_0);
      }

      for (r1 = 0; r1 < 16; r1++) {
        W11[r1] = obj_0->P[r1];
      }

      for (r1 = 0; r1 < 4; r1++) {
        for (ret = 0; ret <= 2; ret += 2) {
          _mm_storeu_pd(&G[ret + (r1 << 2)], _mm_set1_pd(0.0));
          tmp_0 = _mm_loadu_pd(&A[ret]);
          tmp_0 = _mm_mul_pd(_mm_set1_pd(W11[r1 << 2]), tmp_0);
          tmp_1 = _mm_loadu_pd(&G[(r1 << 2) + ret]);
          tmp_0 = _mm_add_pd(tmp_0, tmp_1);
          _mm_storeu_pd(&G[ret + (r1 << 2)], tmp_0);
          tmp_0 = _mm_loadu_pd(&A[ret + 4]);
          tmp_0 = _mm_mul_pd(_mm_set1_pd(W11[(r1 << 2) + 1]), tmp_0);
          tmp_1 = _mm_loadu_pd(&G[(r1 << 2) + ret]);
          tmp_0 = _mm_add_pd(tmp_0, tmp_1);
          _mm_storeu_pd(&G[ret + (r1 << 2)], tmp_0);
          tmp_0 = _mm_loadu_pd(&A[ret + 8]);
          tmp_0 = _mm_mul_pd(_mm_set1_pd(W11[(r1 << 2) + 2]), tmp_0);
          tmp_1 = _mm_loadu_pd(&G[(r1 << 2) + ret]);
          tmp_0 = _mm_add_pd(tmp_0, tmp_1);
          _mm_storeu_pd(&G[ret + (r1 << 2)], tmp_0);
          tmp_0 = _mm_loadu_pd(&A[ret + 12]);
          tmp_0 = _mm_mul_pd(_mm_set1_pd(W11[(r1 << 2) + 3]), tmp_0);
          tmp_1 = _mm_loadu_pd(&G[(r1 << 2) + ret]);
          tmp_0 = _mm_add_pd(tmp_0, tmp_1);
          _mm_storeu_pd(&G[ret + (r1 << 2)], tmp_0);
        }
      }

      for (r1 = 0; r1 < 16; r1++) {
        obj_0->P[r1] = G[r1];
      }

      obj_0->xhat_prev[0] = obj_0->xh[0];
      obj_0->xhat_prev[1] = obj_0->xh[1];
      obj_0->xhat_prev[2] = obj_0->xh[2];
      obj_0->xhat_prev[3] = obj_0->xh[3];
      x_pred[0] = obj_0->xh[0];
      x_pred[1] = obj_0->xh[1];
      x_pred[2] = obj_0->xh[2];
      x_pred[3] = obj_0->xh[3];
      obj->x_hat[0] = x_pred[0];
      obj->x_hat[1] = x_pred[1];
      obj->x_hat[2] = x_pred[2];
      obj->x_hat[3] = x_pred[3];
    } else {
      /*  error("invalid observer") */
      printf("%s\n", "invalid observer");
      fflush(stdout);
    }
  }

  x_hat[0] = obj->x_hat[0];
  x_hat[1] = obj->x_hat[1];
  x_hat[2] = obj->x_hat[2];
  x_hat[3] = obj->x_hat[3];

  /*  Compute control */
  /*  V_servo = obj.feedbackLinearizationController(p_ball, v_ball, theta, dtheta, ... */
  /*      p_ball_ref, v_ball_ref, a_ball_ref); */
  for (r1 = 0; r1 < 6; r1++) {
    a_1[r1] = obj->controller[r1];
  }

  for (r1 = 0; r1 < 6; r1++) {
    b_0[r1] = tmp_2[r1];
  }

  ret = memcmp(&a_1[0], &b_0[0], 6);
  equal = (ret == 0);
  if (equal) {
    obj_0 = obj;
    if (t < 5.0) {
      t_prev = 0.0;
      p_ball_0 = 0.0;
      a22 = 0.0;
    } else if (t < 61.85) {
      theta_0 = t - 5.0;
      t_prev = theta_0 / 56.85;
      if (t_prev < 0.5) {
        t6 = t_prev / 0.5 * 0.090000000000000011 + 0.05;
        t_prev = 0.11423973285781065 * theta_0;
        t_prev = sin(t_prev);
        t_prev = 0.83775804095727813 * theta_0 - 0.2094395102393195 * t_prev /
          0.11423973285781065;
        t_prev = sin(t_prev);
        u_prev = 0.11423973285781065 * theta_0;
        u_prev = sin(u_prev);
        u_prev = 0.83775804095727813 * theta_0 - 0.2094395102393195 * u_prev /
          0.11423973285781065;
        u_prev = cos(u_prev);
        t2 = 0.11423973285781065 * theta_0;
        t2 = cos(t2);
        p_ball_0 = (0.83775804095727813 - 0.2094395102393195 * t2) * (t6 *
          u_prev) + 0.00316622691292876 * t_prev;
        t_prev = 6.2831853071795862 * theta_0 / 55.0;
        t_prev = cos(t_prev);
        t2 = 0.83775804095727813 - 3.1415926535897931 * t_prev / 15.0;
        a22 = t2 * t2;
        t_prev = 6.2831853071795862 * theta_0 / 55.0;
        t_prev = sin(t_prev);
        t_prev = 11.0 * t_prev / 6.0 - 12.566370614359172 * theta_0 / 15.0;
        t_prev = cos(t_prev);
        u_prev = 6.2831853071795862 * theta_0 / 55.0;
        u_prev = cos(u_prev);
        t2 = 6.2831853071795862 * theta_0 / 55.0;
        t2 = sin(t2);
        t2 = 11.0 * t2 / 6.0 - 12.566370614359172 * theta_0 / 15.0;
        t2 = sin(t2);
        t10 = 6.2831853071795862 * theta_0 / 55.0;
        t10 = sin(t10);
        t7 = 6.2831853071795862 * theta_0 / 55.0;
        t7 = sin(t7);
        t7 = 11.0 * t7 / 6.0 - 12.566370614359172 * theta_0 / 15.0;
        t7 = cos(t7);
        a22 = ((0.83775804095727813 - 3.1415926535897931 * u_prev / 15.0) *
               (12.0 * t_prev) / 1895.0 + (6.0 * theta_0 / 1895.0 + 0.05) * t2 *
               a22) + (6.0 * theta_0 / 1895.0 + 0.05) * (19.739208802178716 *
          t10 * t7) / 825.0;
      } else {
        t6 = 0.14;
        t_prev = 0.11423973285781065 * theta_0;
        t_prev = sin(t_prev);
        t_prev = 0.83775804095727813 * theta_0 - 0.2094395102393195 * t_prev /
          0.11423973285781065;
        t_prev = cos(t_prev);
        u_prev = 0.11423973285781065 * theta_0;
        u_prev = cos(u_prev);
        p_ball_0 = (0.83775804095727813 - 0.2094395102393195 * u_prev) * (0.14 *
          t_prev);
        t_prev = 6.2831853071795862 * theta_0 / 55.0;
        t_prev = cos(t_prev);
        t2 = 0.83775804095727813 - 3.1415926535897931 * t_prev / 15.0;
        a22 = t2 * t2;
        t_prev = 6.2831853071795862 * theta_0 / 55.0;
        t_prev = sin(t_prev);
        t_prev = 11.0 * t_prev / 6.0 - 12.566370614359172 * theta_0 / 15.0;
        t_prev = sin(t_prev);
        u_prev = 6.2831853071795862 * theta_0 / 55.0;
        u_prev = sin(u_prev);
        t2 = 6.2831853071795862 * theta_0 / 55.0;
        t2 = sin(t2);
        t2 = 11.0 * t2 / 6.0 - 12.566370614359172 * theta_0 / 15.0;
        t2 = cos(t2);
        a22 = 7.0 * t_prev * a22 / 50.0 + 69.0872308076255 * u_prev * t2 /
          20625.0;
      }

      t_prev = 0.11423973285781065 * theta_0;
      t_prev = sin(t_prev);
      t_prev = 0.83775804095727813 * theta_0 - 0.2094395102393195 * t_prev /
        0.11423973285781065;
      t_prev = sin(t_prev);
      t_prev *= t6;
    } else if (t < 65.0) {
      t_prev = 0.0;
      p_ball_0 = 0.0;
      a22 = 0.0;
    } else if (t < 85.0) {
      t_prev = t - 65.0;
      a22 = t_prev / 20.0;
      if (a22 < 0.5) {
        a22 = 0.05;
      } else {
        a22 = 0.1;
      }

      t_prev *= 0.62831853071795862;
      t_prev = sin(t_prev);
      if (t_prev < 0.0) {
        t_prev = -1.0;
      } else {
        t_prev = (t_prev > 0.0);
      }

      t_prev *= a22;
      p_ball_0 = 0.0;
      a22 = 0.0;
    } else {
      t_prev = 0.0;
      p_ball_0 = 0.0;
      a22 = 0.0;
    }

    obj_1[0] = obj_0->x_eq[2];
    obj_1[1] = obj_0->x_eq[3];
    obj_1[2] = obj_0->x_eq[4];
    simulink_experiment_debu_fsolve(t_prev, p_ball_0, a22, obj_1, mt1);
    mt2[0] = t_prev;
    mt2[1] = p_ball_0;
    mt2[2] = mt1[0];
    mt2[3] = mt1[1];
    mt2[4] = mt1[2];
    for (r1 = 0; r1 < 5; r1++) {
      obj_0->x_eq[r1] = mt2[r1];
    }

    x_pred[0] = mt2[0];
    x_pred[2] = mt2[2];
    x_pred[3] = mt2[3];

    /* A_func */
    /*     A_sym = A_func(IN1,U_SYM) */
    /*     This function was generated by the Symbolic Math Toolbox version 9.2. */
    /*     30-Apr-2025 09:26:50 */
    t_prev = x_pred[0];
    a22 = x_pred[2];
    x4 = x_pred[3];
    t2 = a22;
    t2 = cos(t2);
    u_prev = x4 * x4;
    p_ball_0 = t_prev * 0.7142857142857143;
    theta_0 = t2 * t2;
    t6 = p_ball_0 - 0.15196428571428569;
    a22 = sin(a22);
    A[0] = 0.0;
    A[1] = u_prev * theta_0 * 0.00254530756753206;
    A[2] = 0.0;
    A[3] = 0.0;
    A[4] = 1.0;
    A[5] = 0.0;
    A[6] = 0.0;
    A[7] = 0.0;
    A[8] = 0.0;
    A[9] = t2 * 0.41828772872251141 - t2 * u_prev * t6 * a22 *
      0.007126861189089769;
    A[10] = 0.0;
    A[11] = 0.0;
    A[12] = 0.0;
    A[13] = theta_0 * t6 * x4 * 0.007126861189089769;
    A[14] = 1.0;
    A[15] = -40.0;
    x_pred[0] = 0.0;
    x_pred[1] = 0.0;
    x_pred[2] = 0.0;
    x_pred[3] = 60.0;
    for (r1 = 0; r1 < 16; r1++) {
      W11[r1] = obj_0->Q_tvlqr[r1];
    }

    a22 = obj_0->R_tvlqr;

    /* LQR_VIA_MATRIX_SIGN_FUNCTION Computes the LQR gain using matrix sign function method */
    /*  */
    /*    K = LQR_VIA_MATRIX_SIGN_FUNCTION(A, B, Q, R) returns the optimal gain matrix K */
    /*    for the continuous-time LQR problem: */
    /*  */
    /*        minimize J = &#x222B; (x'Qx + u'Ru) dt */
    /*        subject to dx/dt = Ax + Bu */
    /*  */
    /*    Inputs: */
    /*        A - System dynamics matrix (n x n) */
    /*        B - Input matrix (n x m) */
    /*        Q - State cost matrix (n x n), symmetric positive semi-definite */
    /*        R - Input cost matrix (m x m), symmetric positive definite */
    /*  */
    /*    Output: */
    /*        K - Optimal state feedback gain matrix (m x n) */
    /*  Invert R (assumes R is positive definite) */
    t_prev = 1.0 / a22;
    x_prev[0] = 0.0 * t_prev;
    x_prev[1] = 0.0 * t_prev;
    x_prev[2] = 0.0 * t_prev;
    x_prev[3] = 60.0 * t_prev;

    /*  Construct the Hamiltonian matrix Z */
    for (r1 = 0; r1 < 4; r1++) {
      y_idx_0 = x_pred[r1];
      y_idx_1 = x_prev[0] * y_idx_0;
      Z[r1 << 3] = A[r1 << 2];
      Z[(r1 + 4) << 3] = -y_idx_1;
      Z[(r1 << 3) + 4] = -W11[r1 << 2];
      Z[((r1 + 4) << 3) + 4] = -A[r1];
      y_idx_1 = x_prev[1] * y_idx_0;
      Z[(r1 << 3) + 1] = A[(r1 << 2) + 1];
      Z[((r1 + 4) << 3) + 1] = -y_idx_1;
      Z[(r1 << 3) + 5] = -W11[(r1 << 2) + 1];
      Z[((r1 + 4) << 3) + 5] = -A[r1 + 4];
      y_idx_1 = x_prev[2] * y_idx_0;
      Z[(r1 << 3) + 2] = A[(r1 << 2) + 2];
      Z[((r1 + 4) << 3) + 2] = -y_idx_1;
      Z[(r1 << 3) + 6] = -W11[(r1 << 2) + 2];
      Z[((r1 + 4) << 3) + 6] = -A[r1 + 8];
      y_idx_1 = x_prev[3] * y_idx_0;
      Z[(r1 << 3) + 3] = A[(r1 << 2) + 3];
      Z[((r1 + 4) << 3) + 3] = -y_idx_1;
      Z[(r1 << 3) + 7] = -W11[(r1 << 2) + 3];
      Z[((r1 + 4) << 3) + 7] = -A[r1 + 12];
    }

    /*  Initialize matrix W for iteration */
    /*  Newton iteration to compute the matrix sign function */
    for (ret = 0; ret < 1000; ret++) {
      simulink_experiment_debug_t_inv(Z, tmp);
      for (r1 = 0; r1 <= 62; r1 += 2) {
        tmp_0 = _mm_loadu_pd(&Z[r1]);
        tmp_1 = _mm_loadu_pd(&tmp[r1]);
        tmp_1 = _mm_sub_pd(tmp_0, tmp_1);
        tmp_1 = _mm_mul_pd(_mm_set1_pd(0.5), tmp_1);
        tmp_0 = _mm_sub_pd(tmp_0, tmp_1);
        _mm_storeu_pd(&Z[r1], tmp_0);
      }
    }

    /*  Determine the size of the system */
    /*  Partition W into 4 submatrices */
    for (r1 = 0; r1 < 4; r1++) {
      W11[r1 << 2] = Z[r1 << 3];
      G[r1 << 2] = Z[(r1 + 4) << 3];
      W21[r1 << 2] = Z[(r1 << 3) + 4];
      W22[r1 << 2] = Z[((r1 + 4) << 3) + 4];
      W11[(r1 << 2) + 1] = Z[(r1 << 3) + 1];
      G[(r1 << 2) + 1] = Z[((r1 + 4) << 3) + 1];
      W21[(r1 << 2) + 1] = Z[(r1 << 3) + 5];
      W22[(r1 << 2) + 1] = Z[((r1 + 4) << 3) + 5];
      W11[(r1 << 2) + 2] = Z[(r1 << 3) + 2];
      G[(r1 << 2) + 2] = Z[((r1 + 4) << 3) + 2];
      W21[(r1 << 2) + 2] = Z[(r1 << 3) + 6];
      W22[(r1 << 2) + 2] = Z[((r1 + 4) << 3) + 6];
      W11[(r1 << 2) + 3] = Z[(r1 << 3) + 3];
      G[(r1 << 2) + 3] = Z[((r1 + 4) << 3) + 3];
      W21[(r1 << 2) + 3] = Z[(r1 << 3) + 7];
      W22[(r1 << 2) + 3] = Z[((r1 + 4) << 3) + 7];
    }

    /*  Solve for the unique positive semidefinite solution P to the Riccati equation */
    for (r1 = 0; r1 < 16; r1++) {
      A[r1] = 0.0;
    }

    for (ret = 0; ret < 4; ret++) {
      A[ret + (ret << 2)] = 1.0;
      M[ret << 3] = G[ret << 2];
      M[(ret << 3) + 1] = G[(ret << 2) + 1];
      M[(ret << 3) + 2] = G[(ret << 2) + 2];
      M[(ret << 3) + 3] = G[(ret << 2) + 3];
      M[(ret << 3) + 4] = W22[ret << 2] + A[ret << 2];
      M[(ret << 3) + 5] = W22[(ret << 2) + 1] + A[(ret << 2) + 1];
      M[(ret << 3) + 6] = W22[(ret << 2) + 2] + A[(ret << 2) + 2];
      M[(ret << 3) + 7] = W22[(ret << 2) + 3] + A[(ret << 2) + 3];
    }

    for (r1 = 0; r1 < 16; r1++) {
      A[r1] = 0.0;
    }

    A[0] = 1.0;
    A[5] = 1.0;
    A[10] = 1.0;
    A[15] = 1.0;
    for (r1 = 0; r1 < 4; r1++) {
      N[r1 << 3] = W11[r1 << 2] + A[r1 << 2];
      N[(r1 << 3) + 4] = W21[r1 << 2];
      N[(r1 << 3) + 1] = W11[(r1 << 2) + 1] + A[(r1 << 2) + 1];
      N[(r1 << 3) + 5] = W21[(r1 << 2) + 1];
      N[(r1 << 3) + 2] = W11[(r1 << 2) + 2] + A[(r1 << 2) + 2];
      N[(r1 << 3) + 6] = W21[(r1 << 2) + 2];
      N[(r1 << 3) + 3] = W11[(r1 << 2) + 3] + A[(r1 << 2) + 3];
      N[(r1 << 3) + 7] = W21[(r1 << 2) + 3];
    }

    for (r1 = 0; r1 <= 30; r1 += 2) {
      tmp_0 = _mm_loadu_pd(&N[r1]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(-1.0));
      _mm_storeu_pd(&N_0[r1], tmp_0);
    }

    simulink_experiment_de_mldivide(M, N_0, A);

    /*  Compute the optimal LQR gain */
    y_idx_0 = t_prev * 0.0;
    x_pred[0] = y_idx_0;
    y_idx_0 = t_prev * 0.0;
    x_pred[1] = y_idx_0;
    y_idx_0 = t_prev * 0.0;
    x_pred[2] = y_idx_0;
    y_idx_0 = t_prev * 60.0;
    x_pred[3] = y_idx_0;
    for (r1 = 0; r1 < 4; r1++) {
      t_prev = A[r1 << 2] * x_pred[0];
      t_prev += A[(r1 << 2) + 1] * x_pred[1];
      t_prev += A[(r1 << 2) + 2] * x_pred[2];
      t_prev += A[(r1 << 2) + 3] * x_pred[3];
      t_prev = -t_prev;
      x_prev[r1] = t_prev;
    }

    x_pred[0] = x_hat[0] - mt2[0];
    x_pred[1] = x_hat[1] - mt2[1];
    x_pred[2] = x_hat[2] - mt2[2];
    x_pred[3] = x_hat[3] - mt2[3];
    p_ball_0 = x_prev[0] * x_pred[0];
    p_ball_0 += x_prev[1] * x_pred[1];
    p_ball_0 += x_prev[2] * x_pred[2];
    p_ball_0 += x_prev[3] * x_pred[3];
    *V_servo = p_ball_0 + mt2[4];
  } else {
    /*  error("invalid controller") */
    printf("%s\n", "invalid controller");
    fflush(stdout);
  }

  /*  theta_saturation = 56 * pi / 180; */
  /*  V_servo = clip(V_servo, -10, 10); */
  t_prev = -2.0;
  if (*V_servo >= -2.0) {
    t_prev = *V_servo;
  }

  if (t_prev <= 3.0) {
    *V_servo = t_prev;
  } else {
    *V_servo = 3.0;
  }

  /*  disp(V_servo); */
  /*  disp([p_ball, v_ball, theta, dtheta]); */
  /*  disp([p_ball_obs, v_ball_obs, theta_obs, dtheta_obs]); */
  /*  disp([V_servo]); */
  /*  disp("-----"); */
  /*  % (OPTIONAL) Define safe limits for theta */
  /*  theta_min = -3*pi/8;  % Example lower bound */
  /*  theta_max = 3*pi/8;   % Example upper bound */
  /*  gain = 10;  % Scaling factor for control */
  /*  if theta < theta_min */
  /*      V_servo = max(V_servo, gain * (theta_min - theta));  % Proportional correction */
  /*  elseif theta > theta_max */
  /*      V_servo = min(V_servo, gain * (theta_max - theta));  % Proportional correction */
  /*  end */
  obj->u = *V_servo;

  /*  % (DEFAULT) Decide desired servo angle based on simple proportional feedback. */
  /*  k_p = 3; */
  /*  theta_d = - k_p * (p_ball - p_ball_ref); */
  /*   */
  /*  % Make sure that the desired servo angle does not exceed the physical */
  /*  % limit. This part of code is not necessary but highly recommended */
  /*  % because it addresses the actual physical limit of the servo motor. */
  /*  theta_saturation = 56 * pi / 180;     */
  /*  theta_d = min(theta_d, theta_saturation); */
  /*  theta_d = max(theta_d, -theta_saturation); */
  /*   */
  /*  % Simple position control to control servo angle to the desired */
  /*  % position. */
  /*  k_servo = 10; */
  /*  V_servo = k_servo * (theta_d - theta); */
  /*  % Update class properties if necessary. */
  /*  obj.t_prev = t; */
  /*  obj.theta_d = theta_d; */
}

/* Model output function for TID0 */
void simulink_experiment_debug_type1_output0(void) /* Sample time: [0.0s, 0.0s] */
{
  studentControllerInterface_si_T *obj;
  real_T varargin_2[16];
  real_T varargin_2_0[4];
  real_T varargin_2_1[4];
  real_T x1_0[4];
  real_T R_tvlqr;
  real_T V_servo;
  real_T amplitude;
  real_T k;
  real_T u1;
  real_T u2;
  real_T x1;
  real_T *x1_1;
  int32_T b_k;
  boolean_T exitg1;
  boolean_T p;
  boolean_T p_0;
  boolean_T p_1;

  {                                    /* Sample time: [0.0s, 0.0s] */
    rate_monotonic_scheduler();
  }

  /* S-Function (hil_read_encoder_timebase_block): '<S1>/HIL Read Encoder Timebase' */

  /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_read_encoder
      (simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Task, 1,
       &simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Buffer);
    if (result < 0) {
      simulink_experiment_debug_typ_B.HILReadEncoderTimebase = 0;
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
    } else {
      simulink_experiment_debug_typ_B.HILReadEncoderTimebase =
        simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Buffer;
    }
  }

  /* S-Function (hil_read_analog_block): '<S1>/HIL Read Analog' */

  /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Read Analog (hil_read_analog_block) */
  {
    t_error result = hil_read_analog
      (simulink_experiment_debug_ty_DW.HILInitialize_Card,
       &simulink_experiment_debug_typ_P.HILReadAnalog_channels, 1,
       &simulink_experiment_debug_ty_DW.HILReadAnalog_Buffer);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
    }

    simulink_experiment_debug_typ_B.HILReadAnalog =
      simulink_experiment_debug_ty_DW.HILReadAnalog_Buffer;
  }

  /* Gain: '<S1>/BB01 Sensor  Gain (m//V)' */
  simulink_experiment_debug_typ_B.BB01SensorGainmV =
    simulink_experiment_debug_typ_P.BB01SensorGainmV_Gain *
    simulink_experiment_debug_typ_B.HILReadAnalog;

  /* Gain: '<S1>/Encoder Calibration  (rad//count)' */
  simulink_experiment_debug_typ_B.EncoderCalibrationradcount =
    simulink_experiment_debug_typ_P.EncoderCalibrationradcount_Gain *
    simulink_experiment_debug_typ_B.HILReadEncoderTimebase;

  /* Bias: '<S1>/Bias' */
  simulink_experiment_debug_typ_B.Bias =
    simulink_experiment_debug_typ_B.EncoderCalibrationradcount +
    simulink_experiment_debug_typ_P.Bias_Bias;

  /* Clock: '<Root>/Clock' */
  simulink_experiment_debug_typ_B.Clock =
    simulink_experiment_debug_ty_M->Timing.t[0];

  /* MATLABSystem: '<Root>/MATLAB System' */
  memcpy(&varargin_2[0], &simulink_experiment_debug_typ_P.MATLABSystem_Q_tvlqr[0],
         sizeof(real_T) << 4U);
  R_tvlqr = simulink_experiment_debug_typ_P.MATLABSystem_R_tvlqr;
  V_servo = simulink_experiment_debug_typ_P.MATLABSystem_V_servo;
  varargin_2_0[0] = simulink_experiment_debug_typ_P.MATLABSystem_X_prior[0];
  varargin_2_1[0] = simulink_experiment_debug_typ_P.MATLABSystem_initialState[0];
  varargin_2_0[1] = simulink_experiment_debug_typ_P.MATLABSystem_X_prior[1];
  varargin_2_1[1] = simulink_experiment_debug_typ_P.MATLABSystem_initialState[1];
  varargin_2_0[2] = simulink_experiment_debug_typ_P.MATLABSystem_X_prior[2];
  varargin_2_1[2] = simulink_experiment_debug_typ_P.MATLABSystem_initialState[2];
  varargin_2_0[3] = simulink_experiment_debug_typ_P.MATLABSystem_X_prior[3];
  varargin_2_1[3] = simulink_experiment_debug_typ_P.MATLABSystem_initialState[3];
  amplitude = simulink_experiment_debug_typ_B.Clock;
  u1 = simulink_experiment_debug_typ_B.BB01SensorGainmV;
  u2 = simulink_experiment_debug_typ_B.Bias;
  x1_1 = &simulink_experiment_debug_ty_DW.obj.Q_tvlqr[0];
  p = false;
  p_0 = true;
  b_k = 0;
  exitg1 = false;
  while ((!exitg1) && (b_k < 16)) {
    k = (real_T)b_k + 1.0;
    x1 = x1_1[(int32_T)k - 1];
    k = varargin_2[(int32_T)k - 1];
    p_1 = (x1 == k);
    if (!p_1) {
      p_0 = false;
      exitg1 = true;
    } else {
      b_k++;
    }
  }

  if (p_0) {
    p = true;
  }

  if (!p) {
    memcpy(&simulink_experiment_debug_ty_DW.obj.Q_tvlqr[0], &varargin_2[0],
           sizeof(real_T) << 4U);
  }

  if (simulink_experiment_debug_ty_DW.obj.R_tvlqr != R_tvlqr) {
    simulink_experiment_debug_ty_DW.obj.R_tvlqr = R_tvlqr;
  }

  x1_0[0] = simulink_experiment_debug_ty_DW.obj.X_prior[0];
  x1_0[1] = simulink_experiment_debug_ty_DW.obj.X_prior[1];
  x1_0[2] = simulink_experiment_debug_ty_DW.obj.X_prior[2];
  x1_0[3] = simulink_experiment_debug_ty_DW.obj.X_prior[3];
  p = false;
  p_0 = true;
  b_k = 0;
  exitg1 = false;
  while ((!exitg1) && (b_k < 4)) {
    k = (real_T)b_k + 1.0;
    x1 = x1_0[(int32_T)k - 1];
    k = varargin_2_0[(int32_T)k - 1];
    p_1 = (x1 == k);
    if (!p_1) {
      p_0 = false;
      exitg1 = true;
    } else {
      b_k++;
    }
  }

  if (p_0) {
    p = true;
  }

  if (!p) {
    simulink_experiment_debug_ty_DW.obj.X_prior[0] = varargin_2_0[0];
    simulink_experiment_debug_ty_DW.obj.X_prior[1] = varargin_2_0[1];
    simulink_experiment_debug_ty_DW.obj.X_prior[2] = varargin_2_0[2];
    simulink_experiment_debug_ty_DW.obj.X_prior[3] = varargin_2_0[3];
  }

  if (simulink_experiment_debug_ty_DW.obj.V_servo != V_servo) {
    simulink_experiment_debug_ty_DW.obj.V_servo = V_servo;
  }

  x1_0[0] = simulink_experiment_debug_ty_DW.obj.initialState[0];
  x1_0[1] = simulink_experiment_debug_ty_DW.obj.initialState[1];
  x1_0[2] = simulink_experiment_debug_ty_DW.obj.initialState[2];
  x1_0[3] = simulink_experiment_debug_ty_DW.obj.initialState[3];
  p = false;
  p_0 = true;
  b_k = 0;
  exitg1 = false;
  while ((!exitg1) && (b_k < 4)) {
    k = (real_T)b_k + 1.0;
    x1 = x1_0[(int32_T)k - 1];
    k = varargin_2_1[(int32_T)k - 1];
    p_1 = (x1 == k);
    if (!p_1) {
      p_0 = false;
      exitg1 = true;
    } else {
      b_k++;
    }
  }

  if (p_0) {
    p = true;
  }

  if (!p) {
    simulink_experiment_debug_ty_DW.obj.initialState[0] = varargin_2_1[0];
    simulink_experiment_debug_ty_DW.obj.initialState[1] = varargin_2_1[1];
    simulink_experiment_debug_ty_DW.obj.initialState[2] = varargin_2_1[2];
    simulink_experiment_debug_ty_DW.obj.initialState[3] = varargin_2_1[3];
  }

  obj = &simulink_experiment_debug_ty_DW.obj;
  studentControllerInterface_step(obj, amplitude, u1, u2, &R_tvlqr, varargin_2);

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_debug_typ_B.MATLABSystem_o1 = R_tvlqr;

  /* MATLABSystem: '<Root>/MATLAB System' */
  memcpy(&simulink_experiment_debug_typ_B.MATLABSystem_o2[0], &varargin_2[0],
         sizeof(real_T) << 4U);

  /* Saturate: '<Root>/+//-10V' */
  amplitude = simulink_experiment_debug_typ_B.MATLABSystem_o1;
  u1 = simulink_experiment_debug_typ_P.u0V_LowerSat;
  u2 = simulink_experiment_debug_typ_P.u0V_UpperSat;
  if (amplitude > u2) {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_debug_typ_B.u0V = u2;
  } else if (amplitude < u1) {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_debug_typ_B.u0V = u1;
  } else {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_debug_typ_B.u0V = amplitude;
  }

  /* End of Saturate: '<Root>/+//-10V' */

  /* Gain: '<S1>/Motor  Gain (V//V)' */
  simulink_experiment_debug_typ_B.MotorGainVV =
    simulink_experiment_debug_typ_P.MotorGainVV_Gain *
    simulink_experiment_debug_typ_B.u0V;

  /* S-Function (hil_write_analog_block): '<S1>/HIL Write Analog' */

  /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Write Analog (hil_write_analog_block) */
  {
    t_error result;
    result = hil_write_analog(simulink_experiment_debug_ty_DW.HILInitialize_Card,
      &simulink_experiment_debug_typ_P.HILWriteAnalog_channels, 1,
      &simulink_experiment_debug_typ_B.MotorGainVV);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
    }
  }

  /* MATLAB Function: '<Root>/MATLAB Function' */
  /* MATLAB Function 'MATLAB Function': '<S2>:1' */
  /* '<S2>:1:3' */
  if (simulink_experiment_debug_typ_B.Clock < 5.0) {
    simulink_experiment_debug_typ_B.p_ref = 0.0;
    simulink_experiment_debug_typ_B.v_ref = 0.0;
    simulink_experiment_debug_typ_B.a_ref = 0.0;
  } else if (simulink_experiment_debug_typ_B.Clock < 61.85) {
    amplitude = (simulink_experiment_debug_typ_B.Clock - 5.0) / 56.85;
    if (amplitude < 0.5) {
      amplitude = amplitude / 0.5 * 0.090000000000000011 + 0.05;
      simulink_experiment_debug_typ_B.v_ref = cos
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 -
         sin((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065)
         * 0.2094395102393195 / 0.11423973285781065) * amplitude *
        (0.83775804095727813 - cos((simulink_experiment_debug_typ_B.Clock - 5.0)
          * 0.11423973285781065) * 0.2094395102393195) + sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 -
         sin((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065)
         * 0.2094395102393195 / 0.11423973285781065) * 0.00316622691292876;
      u1 = 0.83775804095727813 - cos((simulink_experiment_debug_typ_B.Clock -
        5.0) * 6.2831853071795862 / 55.0) * 3.1415926535897931 / 15.0;
      simulink_experiment_debug_typ_B.a_ref = (cos(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * 12.0 * (0.83775804095727813 - cos
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 3.1415926535897931 / 15.0) / 1895.0 + sin(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * ((simulink_experiment_debug_typ_B.Clock -
        5.0) * 6.0 / 1895.0 + 0.05) * (u1 * u1)) + cos(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * (sin((simulink_experiment_debug_typ_B.Clock
        - 5.0) * 6.2831853071795862 / 55.0) * 19.739208802178716) *
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.0 / 1895.0 + 0.05) /
        825.0;
    } else {
      amplitude = 0.14;
      simulink_experiment_debug_typ_B.v_ref = cos
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 -
         sin((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065)
         * 0.2094395102393195 / 0.11423973285781065) * 0.14 *
        (0.83775804095727813 - cos((simulink_experiment_debug_typ_B.Clock - 5.0)
          * 0.11423973285781065) * 0.2094395102393195);
      u1 = 0.83775804095727813 - cos((simulink_experiment_debug_typ_B.Clock -
        5.0) * 6.2831853071795862 / 55.0) * 3.1415926535897931 / 15.0;
      simulink_experiment_debug_typ_B.a_ref = sin(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * 7.0 * (u1 * u1) / 50.0 + cos(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * (sin((simulink_experiment_debug_typ_B.Clock
        - 5.0) * 6.2831853071795862 / 55.0) * 69.0872308076255) / 20625.0;
    }

    simulink_experiment_debug_typ_B.p_ref = sin
      ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 - sin
       ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065) *
       0.2094395102393195 / 0.11423973285781065) * amplitude;
  } else if (simulink_experiment_debug_typ_B.Clock < 65.0) {
    simulink_experiment_debug_typ_B.p_ref = 0.0;
    simulink_experiment_debug_typ_B.v_ref = 0.0;
    simulink_experiment_debug_typ_B.a_ref = 0.0;
  } else if (simulink_experiment_debug_typ_B.Clock < 85.0) {
    if ((simulink_experiment_debug_typ_B.Clock - 65.0) / 20.0 < 0.5) {
      amplitude = 0.05;
    } else {
      amplitude = 0.1;
    }

    u1 = sin((simulink_experiment_debug_typ_B.Clock - 65.0) *
             0.62831853071795862);
    if (rtIsNaN(u1)) {
      u1 = (rtNaN);
    } else if (u1 < 0.0) {
      u1 = -1.0;
    } else {
      u1 = (u1 > 0.0);
    }

    simulink_experiment_debug_typ_B.p_ref = amplitude * u1;
    simulink_experiment_debug_typ_B.v_ref = 0.0;
    simulink_experiment_debug_typ_B.a_ref = 0.0;
  } else {
    simulink_experiment_debug_typ_B.p_ref = 0.0;
    simulink_experiment_debug_typ_B.v_ref = 0.0;
    simulink_experiment_debug_typ_B.a_ref = 0.0;
  }

  /* End of MATLAB Function: '<Root>/MATLAB Function' */

  /* Gain: '<Root>/m to cm' */
  /* '<S2>:1:3' */
  simulink_experiment_debug_typ_B.mtocm[0] =
    simulink_experiment_debug_typ_P.mtocm_Gain *
    simulink_experiment_debug_typ_B.p_ref;
  simulink_experiment_debug_typ_B.mtocm[1] =
    simulink_experiment_debug_typ_P.mtocm_Gain *
    simulink_experiment_debug_typ_B.BB01SensorGainmV;

  /* Gain: '<S3>/Gain' */
  simulink_experiment_debug_typ_B.Gain =
    simulink_experiment_debug_typ_P.Gain_Gain *
    simulink_experiment_debug_typ_B.Bias;

  /* RateTransition: '<Root>/Rate Transition' */
  if (simulink_experiment_debug_ty_M->Timing.RateInteraction.TID1_2) {
    simulink_experiment_debug_ty_DW.RateTransition_Buffer =
      simulink_experiment_debug_typ_B.Clock;

    /* RateTransition: '<Root>/Rate Transition1' */
    simulink_experiment_debug_ty_DW.RateTransition1_Buffer =
      simulink_experiment_debug_typ_B.p_ref;

    /* RateTransition: '<Root>/Rate Transition2' */
    simulink_experiment_debug_ty_DW.RateTransition2_Buffer =
      simulink_experiment_debug_typ_B.MATLABSystem_o1;

    /* RateTransition: '<Root>/Rate Transition3' */
    simulink_experiment_debug_ty_DW.RateTransition3_Buffer =
      simulink_experiment_debug_typ_B.BB01SensorGainmV;

    /* RateTransition: '<Root>/Rate Transition4' */
    simulink_experiment_debug_ty_DW.RateTransition4_Buffer =
      simulink_experiment_debug_typ_B.Bias;
  }

  /* End of RateTransition: '<Root>/Rate Transition' */
}

/* Model update function for TID0 */
void simulink_experiment_debug_type1_update0(void) /* Sample time: [0.0s, 0.0s] */
{
  /* Update absolute time */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++simulink_experiment_debug_ty_M->Timing.clockTick0)) {
    ++simulink_experiment_debug_ty_M->Timing.clockTickH0;
  }

  simulink_experiment_debug_ty_M->Timing.t[0] =
    simulink_experiment_debug_ty_M->Timing.clockTick0 *
    simulink_experiment_debug_ty_M->Timing.stepSize0 +
    simulink_experiment_debug_ty_M->Timing.clockTickH0 *
    simulink_experiment_debug_ty_M->Timing.stepSize0 * 4294967296.0;

  /* Update absolute time */
  /* The "clockTick1" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick1"
   * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick1 and the high bits
   * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++simulink_experiment_debug_ty_M->Timing.clockTick1)) {
    ++simulink_experiment_debug_ty_M->Timing.clockTickH1;
  }

  simulink_experiment_debug_ty_M->Timing.t[1] =
    simulink_experiment_debug_ty_M->Timing.clockTick1 *
    simulink_experiment_debug_ty_M->Timing.stepSize1 +
    simulink_experiment_debug_ty_M->Timing.clockTickH1 *
    simulink_experiment_debug_ty_M->Timing.stepSize1 * 4294967296.0;
}

/* Model output function for TID2 */
void simulink_experiment_debug_type1_output2(void) /* Sample time: [0.01s, 0.0s] */
{
  /* RateTransition: '<Root>/Rate Transition2' */
  simulink_experiment_debug_typ_B.RateTransition2 =
    simulink_experiment_debug_ty_DW.RateTransition2_Buffer;

  /* RateTransition: '<Root>/Rate Transition1' */
  simulink_experiment_debug_typ_B.RateTransition1 =
    simulink_experiment_debug_ty_DW.RateTransition1_Buffer;

  /* RateTransition: '<Root>/Rate Transition3' */
  simulink_experiment_debug_typ_B.RateTransition3 =
    simulink_experiment_debug_ty_DW.RateTransition3_Buffer;

  /* RateTransition: '<Root>/Rate Transition4' */
  simulink_experiment_debug_typ_B.RateTransition4 =
    simulink_experiment_debug_ty_DW.RateTransition4_Buffer;

  /* RateTransition: '<Root>/Rate Transition' */
  simulink_experiment_debug_typ_B.RateTransition =
    simulink_experiment_debug_ty_DW.RateTransition_Buffer;
}

/* Model update function for TID2 */
void simulink_experiment_debug_type1_update2(void) /* Sample time: [0.01s, 0.0s] */
{
  /* Update absolute time */
  /* The "clockTick2" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick2"
   * and "Timing.stepSize2". Size of "clockTick2" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick2 and the high bits
   * Timing.clockTickH2. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++simulink_experiment_debug_ty_M->Timing.clockTick2)) {
    ++simulink_experiment_debug_ty_M->Timing.clockTickH2;
  }

  simulink_experiment_debug_ty_M->Timing.t[2] =
    simulink_experiment_debug_ty_M->Timing.clockTick2 *
    simulink_experiment_debug_ty_M->Timing.stepSize2 +
    simulink_experiment_debug_ty_M->Timing.clockTickH2 *
    simulink_experiment_debug_ty_M->Timing.stepSize2 * 4294967296.0;
}

/* Use this function only if you need to maintain compatibility with an existing static main program. */
void simulink_experiment_debug_type1_output(int_T tid)
{
  switch (tid) {
   case 0 :
    simulink_experiment_debug_type1_output0();
    break;

   case 2 :
    simulink_experiment_debug_type1_output2();
    break;

   default :
    /* do nothing */
    break;
  }
}

/* Use this function only if you need to maintain compatibility with an existing static main program. */
void simulink_experiment_debug_type1_update(int_T tid)
{
  switch (tid) {
   case 0 :
    simulink_experiment_debug_type1_update0();
    break;

   case 2 :
    simulink_experiment_debug_type1_update2();
    break;

   default :
    /* do nothing */
    break;
  }
}

/* Model initialize function */
void simulink_experiment_debug_type1_initialize(void)
{
  {
    studentControllerInterface_si_T *obj;
    const real_T *varargin_1;
    real_T R_tvlqr;
    real_T V_servo;
    real_T varargin_3_idx_0;
    real_T varargin_3_idx_1;
    real_T varargin_3_idx_2;
    real_T varargin_3_idx_3;
    real_T varargin_5_idx_0;
    real_T varargin_5_idx_1;
    real_T varargin_5_idx_2;
    real_T varargin_5_idx_3;

    /* Start for S-Function (hil_initialize_block): '<S1>/HIL Initialize' */

    /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Initialize (hil_initialize_block) */
    {
      t_int result;
      t_boolean is_switching;
      result = hil_open("q2_usb", "0",
                        &simulink_experiment_debug_ty_DW.HILInitialize_Card);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        return;
      }

      is_switching = false;
      result = hil_set_card_specific_options
        (simulink_experiment_debug_ty_DW.HILInitialize_Card,
         "d0=digital;d1=digital;led=auto;update_rate=normal", 50);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        return;
      }

      result = hil_watchdog_clear
        (simulink_experiment_debug_ty_DW.HILInitialize_Card);
      if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        return;
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_AIPStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_AIPEnter &&
           is_switching)) {
        simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[0] =
          (simulink_experiment_debug_typ_P.HILInitialize_AILow);
        simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[1] =
          (simulink_experiment_debug_typ_P.HILInitialize_AILow);
        simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[0] =
          simulink_experiment_debug_typ_P.HILInitialize_AIHigh;
        simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[1] =
          simulink_experiment_debug_typ_P.HILInitialize_AIHigh;
        result = hil_set_analog_input_ranges
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_AIChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[0],
           &simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_AOPStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_AOPEnter &&
           is_switching)) {
        simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[0] =
          (simulink_experiment_debug_typ_P.HILInitialize_AOLow);
        simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[1] =
          (simulink_experiment_debug_typ_P.HILInitialize_AOLow);
        simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[0] =
          simulink_experiment_debug_typ_P.HILInitialize_AOHigh;
        simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[1] =
          simulink_experiment_debug_typ_P.HILInitialize_AOHigh;
        result = hil_set_analog_output_ranges
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_AOChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[0],
           &simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_AOStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_AOEnter && is_switching))
      {
        simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0] =
          simulink_experiment_debug_typ_P.HILInitialize_AOInitial;
        simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[1] =
          simulink_experiment_debug_typ_P.HILInitialize_AOInitial;
        result = hil_write_analog
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_AOChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if (simulink_experiment_debug_typ_P.HILInitialize_AOReset) {
        simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0] =
          simulink_experiment_debug_typ_P.HILInitialize_AOWatchdog;
        simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[1] =
          simulink_experiment_debug_typ_P.HILInitialize_AOWatchdog;
        result = hil_watchdog_set_analog_expiration_state
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_AOChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      result = hil_set_digital_directions
        (simulink_experiment_debug_ty_DW.HILInitialize_Card, NULL, 0U,
         simulink_experiment_debug_typ_P.HILInitialize_DOChannels, 8U);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        return;
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_DOStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_DOEnter && is_switching))
      {
        {
          int_T i1;
          boolean_T *dw_DOBits =
            &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0];
          for (i1=0; i1 < 8; i1++) {
            dw_DOBits[i1] =
              simulink_experiment_debug_typ_P.HILInitialize_DOInitial;
          }
        }

        result = hil_write_digital
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_DOChannels, 8U,
           (t_boolean *) &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if (simulink_experiment_debug_typ_P.HILInitialize_DOReset) {
        {
          int_T i1;
          int32_T *dw_DOStates =
            &simulink_experiment_debug_ty_DW.HILInitialize_DOStates[0];
          for (i1=0; i1 < 8; i1++) {
            dw_DOStates[i1] =
              simulink_experiment_debug_typ_P.HILInitialize_DOWatchdog;
          }
        }

        result = hil_watchdog_set_digital_expiration_state
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_DOChannels, 8U, (const
            t_digital_state *)
           &simulink_experiment_debug_ty_DW.HILInitialize_DOStates[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_EIPStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_EIPEnter &&
           is_switching)) {
        simulink_experiment_debug_ty_DW.HILInitialize_QuadratureModes[0] =
          simulink_experiment_debug_typ_P.HILInitialize_EIQuadrature;
        simulink_experiment_debug_ty_DW.HILInitialize_QuadratureModes[1] =
          simulink_experiment_debug_typ_P.HILInitialize_EIQuadrature;
        result = hil_set_encoder_quadrature_mode
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_EIChannels, 2U,
           (t_encoder_quadrature_mode *)
           &simulink_experiment_debug_ty_DW.HILInitialize_QuadratureModes[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_EIStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_EIEnter && is_switching))
      {
        simulink_experiment_debug_ty_DW.HILInitialize_InitialEICounts[0] =
          simulink_experiment_debug_typ_P.HILInitialize_EIInitial;
        simulink_experiment_debug_ty_DW.HILInitialize_InitialEICounts[1] =
          simulink_experiment_debug_typ_P.HILInitialize_EIInitial;
        result = hil_set_encoder_counts
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_EIChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_InitialEICounts[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }
    }

    /* Start for S-Function (hil_read_encoder_timebase_block): '<S1>/HIL Read Encoder Timebase' */

    /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_create_encoder_reader
        (simulink_experiment_debug_ty_DW.HILInitialize_Card,
         simulink_experiment_debug_typ_P.HILReadEncoderTimebase_SamplesI,
         &simulink_experiment_debug_typ_P.HILReadEncoderTimebase_Channels, 1,
         &simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Task);
      if (result >= 0) {
        result = hil_task_set_buffer_overflow_mode
          (simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Task,
           (t_buffer_overflow_mode)
           (simulink_experiment_debug_typ_P.HILReadEncoderTimebase_Overflow - 1));
      }

      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
      }
    }

    /* Start for MATLABSystem: '<Root>/MATLAB System' */
    varargin_1 = &simulink_experiment_debug_typ_P.MATLABSystem_Q_tvlqr[0];
    R_tvlqr = simulink_experiment_debug_typ_P.MATLABSystem_R_tvlqr;
    V_servo = simulink_experiment_debug_typ_P.MATLABSystem_V_servo;
    varargin_3_idx_0 = simulink_experiment_debug_typ_P.MATLABSystem_X_prior[0];
    varargin_5_idx_0 =
      simulink_experiment_debug_typ_P.MATLABSystem_initialState[0];
    varargin_3_idx_1 = simulink_experiment_debug_typ_P.MATLABSystem_X_prior[1];
    varargin_5_idx_1 =
      simulink_experiment_debug_typ_P.MATLABSystem_initialState[1];
    varargin_3_idx_2 = simulink_experiment_debug_typ_P.MATLABSystem_X_prior[2];
    varargin_5_idx_2 =
      simulink_experiment_debug_typ_P.MATLABSystem_initialState[2];
    varargin_3_idx_3 = simulink_experiment_debug_typ_P.MATLABSystem_X_prior[3];
    varargin_5_idx_3 =
      simulink_experiment_debug_typ_P.MATLABSystem_initialState[3];
    studentControllerInterface_stud(&simulink_experiment_debug_ty_DW.obj);
    simulink_experiment_debug_ty_DW.objisempty = true;
    memcpy(&simulink_experiment_debug_ty_DW.obj.Q_tvlqr[0], &varargin_1[0],
           sizeof(real_T) << 4U);
    simulink_experiment_debug_ty_DW.obj.R_tvlqr = R_tvlqr;
    simulink_experiment_debug_ty_DW.obj.X_prior[0] = varargin_3_idx_0;
    simulink_experiment_debug_ty_DW.obj.X_prior[1] = varargin_3_idx_1;
    simulink_experiment_debug_ty_DW.obj.X_prior[2] = varargin_3_idx_2;
    simulink_experiment_debug_ty_DW.obj.X_prior[3] = varargin_3_idx_3;
    simulink_experiment_debug_ty_DW.obj.V_servo = V_servo;
    simulink_experiment_debug_ty_DW.obj.initialState[0] = varargin_5_idx_0;
    simulink_experiment_debug_ty_DW.obj.initialState[1] = varargin_5_idx_1;
    simulink_experiment_debug_ty_DW.obj.initialState[2] = varargin_5_idx_2;
    simulink_experiment_debug_ty_DW.obj.initialState[3] = varargin_5_idx_3;
    obj = &simulink_experiment_debug_ty_DW.obj;
    obj->isInitialized = 1;
  }
}

/* Model terminate function */
void simulink_experiment_debug_type1_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<S1>/HIL Initialize' */

  /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_digital_outputs = 0;
    hil_task_stop_all(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    hil_monitor_stop_all(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    is_switching = false;
    if ((simulink_experiment_debug_typ_P.HILInitialize_AOTerminate &&
         !is_switching) || (simulink_experiment_debug_typ_P.HILInitialize_AOExit
         && is_switching)) {
      simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0] =
        simulink_experiment_debug_typ_P.HILInitialize_AOFinal;
      simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[1] =
        simulink_experiment_debug_typ_P.HILInitialize_AOFinal;
      num_final_analog_outputs = 2U;
    } else {
      num_final_analog_outputs = 0;
    }

    if ((simulink_experiment_debug_typ_P.HILInitialize_DOTerminate &&
         !is_switching) || (simulink_experiment_debug_typ_P.HILInitialize_DOExit
         && is_switching)) {
      {
        int_T i1;
        boolean_T *dw_DOBits =
          &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0];
        for (i1=0; i1 < 8; i1++) {
          dw_DOBits[i1] = simulink_experiment_debug_typ_P.HILInitialize_DOFinal;
        }
      }

      num_final_digital_outputs = 8U;
    } else {
      num_final_digital_outputs = 0;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_digital_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(simulink_experiment_debug_ty_DW.HILInitialize_Card
                         ,
                         simulink_experiment_debug_typ_P.HILInitialize_AOChannels,
                         num_final_analog_outputs
                         , NULL, 0
                         ,
                         simulink_experiment_debug_typ_P.HILInitialize_DOChannels,
                         num_final_digital_outputs
                         , NULL, 0
                         ,
                         &simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages
                         [0]
                         , NULL
                         , (t_boolean *)
                         &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0]
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog
            (simulink_experiment_debug_ty_DW.HILInitialize_Card,
             simulink_experiment_debug_typ_P.HILInitialize_AOChannels,
             num_final_analog_outputs,
             &simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_digital_outputs > 0) {
          local_result = hil_write_digital
            (simulink_experiment_debug_ty_DW.HILInitialize_Card,
             simulink_experiment_debug_typ_P.HILInitialize_DOChannels,
             num_final_digital_outputs, (t_boolean *)
             &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    hil_monitor_delete_all(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    hil_close(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    simulink_experiment_debug_ty_DW.HILInitialize_Card = NULL;
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/
void MdlOutputs(int_T tid)
{
  if (tid == 1)
    tid = 0;
  simulink_experiment_debug_type1_output(tid);
}

void MdlUpdate(int_T tid)
{
  if (tid == 1)
    tid = 0;
  simulink_experiment_debug_type1_update(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  simulink_experiment_debug_type1_initialize();
}

void MdlTerminate(void)
{
  simulink_experiment_debug_type1_terminate();
}

/* Registration function */
RT_MODEL_simulink_experiment__T *simulink_experiment_debug_type1(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)simulink_experiment_debug_ty_M, 0,
                sizeof(RT_MODEL_simulink_experiment__T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&simulink_experiment_debug_ty_M->solverInfo,
                          &simulink_experiment_debug_ty_M->Timing.simTimeStep);
    rtsiSetTPtr(&simulink_experiment_debug_ty_M->solverInfo, &rtmGetTPtr
                (simulink_experiment_debug_ty_M));
    rtsiSetStepSizePtr(&simulink_experiment_debug_ty_M->solverInfo,
                       &simulink_experiment_debug_ty_M->Timing.stepSize0);
    rtsiSetErrorStatusPtr(&simulink_experiment_debug_ty_M->solverInfo,
                          (&rtmGetErrorStatus(simulink_experiment_debug_ty_M)));
    rtsiSetRTModelPtr(&simulink_experiment_debug_ty_M->solverInfo,
                      simulink_experiment_debug_ty_M);
  }

  rtsiSetSimTimeStep(&simulink_experiment_debug_ty_M->solverInfo,
                     MAJOR_TIME_STEP);
  rtsiSetIsMinorTimeStepWithModeChange
    (&simulink_experiment_debug_ty_M->solverInfo, false);
  rtsiSetSolverName(&simulink_experiment_debug_ty_M->solverInfo,
                    "FixedStepDiscrete");

  /* Initialize timing info */
  {
    int_T *mdlTsMap =
      simulink_experiment_debug_ty_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    mdlTsMap[2] = 2;

    /* polyspace +2 MISRA2012:D4.1 [Justified:Low] "simulink_experiment_debug_ty_M points to
       static memory which is guaranteed to be non-NULL" */
    simulink_experiment_debug_ty_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    simulink_experiment_debug_ty_M->Timing.sampleTimes =
      (&simulink_experiment_debug_ty_M->Timing.sampleTimesArray[0]);
    simulink_experiment_debug_ty_M->Timing.offsetTimes =
      (&simulink_experiment_debug_ty_M->Timing.offsetTimesArray[0]);

    /* task periods */
    simulink_experiment_debug_ty_M->Timing.sampleTimes[0] = (0.0);
    simulink_experiment_debug_ty_M->Timing.sampleTimes[1] = (0.002);
    simulink_experiment_debug_ty_M->Timing.sampleTimes[2] = (0.01);

    /* task offsets */
    simulink_experiment_debug_ty_M->Timing.offsetTimes[0] = (0.0);
    simulink_experiment_debug_ty_M->Timing.offsetTimes[1] = (0.0);
    simulink_experiment_debug_ty_M->Timing.offsetTimes[2] = (0.0);
  }

  rtmSetTPtr(simulink_experiment_debug_ty_M,
             &simulink_experiment_debug_ty_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = simulink_experiment_debug_ty_M->Timing.sampleHitArray;
    int_T *mdlPerTaskSampleHits =
      simulink_experiment_debug_ty_M->Timing.perTaskSampleHitsArray;
    simulink_experiment_debug_ty_M->Timing.perTaskSampleHits =
      (&mdlPerTaskSampleHits[0]);
    mdlSampleHits[0] = 1;
    simulink_experiment_debug_ty_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(simulink_experiment_debug_ty_M, 25.0);
  simulink_experiment_debug_ty_M->Timing.stepSize0 = 0.002;
  simulink_experiment_debug_ty_M->Timing.stepSize1 = 0.002;
  simulink_experiment_debug_ty_M->Timing.stepSize2 = 0.01;

  /* External mode info */
  simulink_experiment_debug_ty_M->Sizes.checksums[0] = (2026849503U);
  simulink_experiment_debug_ty_M->Sizes.checksums[1] = (3021995234U);
  simulink_experiment_debug_ty_M->Sizes.checksums[2] = (192163238U);
  simulink_experiment_debug_ty_M->Sizes.checksums[3] = (3525489453U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[3];
    simulink_experiment_debug_ty_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = &rtAlwaysEnabled;
    systemRan[2] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(simulink_experiment_debug_ty_M->extModeInfo,
      &simulink_experiment_debug_ty_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(simulink_experiment_debug_ty_M->extModeInfo,
                        simulink_experiment_debug_ty_M->Sizes.checksums);
    rteiSetTPtr(simulink_experiment_debug_ty_M->extModeInfo, rtmGetTPtr
                (simulink_experiment_debug_ty_M));
  }

  simulink_experiment_debug_ty_M->solverInfoPtr =
    (&simulink_experiment_debug_ty_M->solverInfo);
  simulink_experiment_debug_ty_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&simulink_experiment_debug_ty_M->solverInfo, 0.002);
  rtsiSetSolverMode(&simulink_experiment_debug_ty_M->solverInfo,
                    SOLVER_MODE_MULTITASKING);

  /* block I/O */
  simulink_experiment_debug_ty_M->blockIO = ((void *)
    &simulink_experiment_debug_typ_B);

  {
    int32_T i;
    for (i = 0; i < 16; i++) {
      simulink_experiment_debug_typ_B.MATLABSystem_o2[i] = 0.0;
    }

    simulink_experiment_debug_typ_B.HILReadEncoderTimebase = 0.0;
    simulink_experiment_debug_typ_B.HILReadAnalog = 0.0;
    simulink_experiment_debug_typ_B.BB01SensorGainmV = 0.0;
    simulink_experiment_debug_typ_B.EncoderCalibrationradcount = 0.0;
    simulink_experiment_debug_typ_B.Bias = 0.0;
    simulink_experiment_debug_typ_B.Clock = 0.0;
    simulink_experiment_debug_typ_B.u0V = 0.0;
    simulink_experiment_debug_typ_B.MotorGainVV = 0.0;
    simulink_experiment_debug_typ_B.mtocm[0] = 0.0;
    simulink_experiment_debug_typ_B.mtocm[1] = 0.0;
    simulink_experiment_debug_typ_B.Gain = 0.0;
    simulink_experiment_debug_typ_B.RateTransition2 = 0.0;
    simulink_experiment_debug_typ_B.RateTransition1 = 0.0;
    simulink_experiment_debug_typ_B.RateTransition3 = 0.0;
    simulink_experiment_debug_typ_B.RateTransition4 = 0.0;
    simulink_experiment_debug_typ_B.RateTransition = 0.0;
    simulink_experiment_debug_typ_B.MATLABSystem_o1 = 0.0;
    simulink_experiment_debug_typ_B.p_ref = 0.0;
    simulink_experiment_debug_typ_B.v_ref = 0.0;
    simulink_experiment_debug_typ_B.a_ref = 0.0;
  }

  /* parameters */
  simulink_experiment_debug_ty_M->defaultParam = ((real_T *)
    &simulink_experiment_debug_typ_P);

  /* states (dwork) */
  simulink_experiment_debug_ty_M->dwork = ((void *)
    &simulink_experiment_debug_ty_DW);
  (void) memset((void *)&simulink_experiment_debug_ty_DW, 0,
                sizeof(DW_simulink_experiment_debug__T));
  simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_FilterFrequency[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_FilterFrequency[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILReadAnalog_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition1_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition2_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition3_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition4_Buffer = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    simulink_experiment_debug_ty_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 22;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Initialize Sizes */
  simulink_experiment_debug_ty_M->Sizes.numContStates = (0);/* Number of continuous states */
  simulink_experiment_debug_ty_M->Sizes.numY = (0);/* Number of model outputs */
  simulink_experiment_debug_ty_M->Sizes.numU = (0);/* Number of model inputs */
  simulink_experiment_debug_ty_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  simulink_experiment_debug_ty_M->Sizes.numSampTimes = (3);/* Number of sample times */
  simulink_experiment_debug_ty_M->Sizes.numBlocks = (30);/* Number of blocks */
  simulink_experiment_debug_ty_M->Sizes.numBlockIO = (20);/* Number of block outputs */
  simulink_experiment_debug_ty_M->Sizes.numBlockPrms = (114);/* Sum of parameter "widths" */
  return simulink_experiment_debug_ty_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
