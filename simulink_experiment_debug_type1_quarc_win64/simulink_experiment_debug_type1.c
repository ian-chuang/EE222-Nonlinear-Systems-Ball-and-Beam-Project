/*
 * simulink_experiment_debug_type1.c
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

#include "simulink_experiment_debug_type1.h"
#include "simulink_experiment_debug_type1_types.h"
#include "rtwtypes.h"
#include <math.h>
#include <emmintrin.h>
#include "simulink_experiment_debug_type1_private.h"
#include "rt_nonfinite.h"
#include <string.h>
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
  static const real_T tmp[16] = { 0.5, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0,
    0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.01 };

  static const real_T tmp_0[16] = { 0.01, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0,
    0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.01 };

  b_obj = obj;
  b_obj->rg_val = 0.0254;
  b_obj->L_val = 0.4255;
  b_obj->g_val = 9.81;
  b_obj->K_val = 1.5;
  b_obj->tau_val = 0.025;
  b_obj->dt = 0.01;
  b_obj->max_V = 3.0;
  b_obj->K[0] = 100.0;
  b_obj->K[1] = 143.0;
  b_obj->K[2] = 52.1;
  b_obj->K[3] = 10.2;
  b_obj->t_prev = 0.0;
  b_obj->u = 0.0;
  for (i = 0; i < 16; i++) {
    b_obj->P[i] = tmp[i];
  }

  b_obj->xhat_prev[0] = 0.0;
  b_obj->xhat_prev[1] = 0.0;
  b_obj->xhat_prev[2] = -0.99483767363676778;
  b_obj->xhat_prev[3] = 0.0;
  for (i = 0; i < 16; i++) {
    b_obj->Q[i] = tmp_0[i];
  }

  b_obj->R[0] = 0.003;
  b_obj->R[1] = 0.0;
  b_obj->R[2] = 0.0;
  b_obj->R[3] = 0.003;
  b_obj->isInitialized = 0;
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

/* Model output function for TID0 */
void simulink_experiment_debug_type1_output0(void) /* Sample time: [0.0s, 0.0s] */
{
  __m128d tmp;
  __m128d tmp_0;
  studentControllerInterface_si_T *obj;
  studentControllerInterface_si_T *varargin_3;
  real_T A[16];
  real_T b[16];
  real_T b_y[16];
  real_T y[16];
  real_T K_gain[8];
  real_T b_y_0[8];
  real_T x_prev[4];
  real_T L;
  real_T ONE;
  real_T a22;
  real_T a_0;
  real_T amp;
  real_T b_0;
  real_T b_x;
  real_T b_x_0;
  real_T b_x_1;
  real_T b_x_2;
  real_T b_x_3;
  real_T b_x_4;
  real_T c;
  real_T c_0;
  real_T c_1;
  real_T c_2;
  real_T c_3;
  real_T c_4;
  real_T c_5;
  real_T c_6;
  real_T c_7;
  real_T c_8;
  real_T c_9;
  real_T c_a;
  real_T c_b;
  real_T c_c;
  real_T c_d;
  real_T c_e;
  real_T phase_zero2_end;
  real_T rg;
  real_T t_sine;
  real_T tau;
  real_T u0;
  real_T u2;
  real_T u_prev;
  real_T x;
  real_T x_0;
  real_T x_pred_idx_0;
  real_T x_pred_idx_1;
  real_T x_pred_idx_2;
  real_T x_pred_idx_3;
  real_T y_idx_1;
  int32_T TWO;
  int32_T r1;
  int8_T a[8];
  static const int8_T tmp_1[8] = { 1, 0, 0, 0, 0, 1, 0, 0 };

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
  u0 = simulink_experiment_debug_typ_B.Clock;
  t_sine = simulink_experiment_debug_typ_B.BB01SensorGainmV;
  u2 = simulink_experiment_debug_typ_B.Bias;
  obj = &simulink_experiment_debug_ty_DW.obj;

  /*  function setupImpl(obj) */
  /*   */
  /*  end */
  /*  This is the main function called every iteration. You have to implement */
  /*  the controller in this function, bu you are not allowed to */
  /*  change the signature of this function.  */
  /*  Input arguments: */
  /*    t: current time */
  /*    p_ball: position of the ball provided by the ball position sensor (m) */
  /*  */
  /*    theta: servo motor angle provided by the encoder of the motor (rad) */
  /*  Output: */
  /*    V_servo: voltage to the servo input.         */
  if (obj->t_prev != 0.0) {
    obj->dt = u0 - obj->t_prev;
  }

  u_prev = obj->u;
  x_prev[0] = obj->xhat_prev[0];
  x_prev[1] = obj->xhat_prev[1];
  x_prev[2] = obj->xhat_prev[2];
  x_prev[3] = obj->xhat_prev[3];
  varargin_3 = obj;
  a_0 = varargin_3->rg_val / varargin_3->L_val * x_prev[3];
  ONE = a_0 * a_0;
  phase_zero2_end = x_prev[2];
  phase_zero2_end = cos(phase_zero2_end);
  a22 = phase_zero2_end * phase_zero2_end;
  a_0 = varargin_3->dt;
  phase_zero2_end = x_prev[2];
  phase_zero2_end = sin(phase_zero2_end);
  x_pred_idx_0 = x_prev[1];
  x_pred_idx_1 = 5.0 * varargin_3->g_val * varargin_3->rg_val / (7.0 *
    varargin_3->L_val) * phase_zero2_end - (varargin_3->L_val / 2.0 - x_prev[0])
    * 0.7142857142857143 * ONE * a22;
  x_pred_idx_2 = x_prev[3];
  x_pred_idx_3 = varargin_3->K_val / varargin_3->tau_val * u_prev + -x_prev[3] /
    varargin_3->tau_val;
  a22 = x_pred_idx_0;
  a22 *= a_0;
  a22 += x_prev[0];
  x_pred_idx_0 = a22;
  a22 = x_pred_idx_1;
  a22 *= a_0;
  a22 += x_prev[1];
  x_pred_idx_1 = a22;
  a22 = x_pred_idx_2;
  a22 *= a_0;
  a22 += x_prev[2];
  x_pred_idx_2 = a22;
  a22 = x_pred_idx_3;
  a22 *= a_0;
  a22 += x_prev[3];
  x_pred_idx_3 = a22;
  b_0 = obj->dt;
  amp = obj->g_val;
  L = obj->L_val;
  rg = obj->rg_val;
  tau = obj->tau_val;
  for (TWO = 0; TWO < 16; TWO++) {
    A[TWO] = 0.0;
  }

  A[0] = 1.0;
  A[5] = 1.0;
  A[10] = 1.0;
  A[4] = b_0;
  a_0 = rg / L;
  ONE = a_0 * a_0;
  a_0 = x_prev[3];
  a22 = a_0 * a_0;
  phase_zero2_end = x_prev[2];
  phase_zero2_end = cos(phase_zero2_end);
  u_prev = phase_zero2_end * phase_zero2_end;
  A[1] = ONE * a22 * u_prev * (0.7142857142857143 * b_0);
  a_0 = rg / L;
  ONE = a_0 * a_0;
  a_0 = x_prev[3];
  a22 = a_0 * a_0;
  phase_zero2_end = x_prev[2];
  phase_zero2_end = cos(phase_zero2_end);
  x = x_prev[2];
  x = cos(x);
  x_0 = x_prev[2];
  x_0 = sin(x_0);
  A[9] = (L / 2.0 - x_prev[0]) * ONE * a22 * x * x_0 * (1.4285714285714286 * b_0)
    + 5.0 * amp * rg / (7.0 * L) * b_0 * phase_zero2_end;
  a_0 = rg / L;
  ONE = a_0 * a_0;
  phase_zero2_end = x_prev[2];
  phase_zero2_end = cos(phase_zero2_end);
  a22 = phase_zero2_end * phase_zero2_end;
  A[13] = (L / 2.0 - x_prev[0]) * ONE * x_prev[3] * a22 * (-1.4285714285714286 *
    b_0);
  A[14] = b_0;
  A[15] = 1.0 - b_0 / tau;
  for (TWO = 0; TWO < 16; TWO++) {
    b[TWO] = obj->P[TWO];
  }

  for (TWO = 0; TWO < 4; TWO++) {
    for (r1 = 0; r1 < 4; r1++) {
      b_y[TWO + (r1 << 2)] = 0.0;
      a22 = b_y[(r1 << 2) + TWO];
      a22 += b[r1 << 2] * A[TWO];
      b_y[TWO + (r1 << 2)] = a22;
      a22 = b_y[(r1 << 2) + TWO];
      a22 += b[(r1 << 2) + 1] * A[TWO + 4];
      b_y[TWO + (r1 << 2)] = a22;
      a22 = b_y[(r1 << 2) + TWO];
      a22 += b[(r1 << 2) + 2] * A[TWO + 8];
      b_y[TWO + (r1 << 2)] = a22;
      a22 = b_y[(r1 << 2) + TWO];
      a22 += b[(r1 << 2) + 3] * A[TWO + 12];
      b_y[TWO + (r1 << 2)] = a22;
    }

    for (r1 = 0; r1 < 4; r1++) {
      y[TWO + (r1 << 2)] = 0.0;
      a_0 = y[(r1 << 2) + TWO];
      a_0 += b_y[TWO] * A[r1];
      y[TWO + (r1 << 2)] = a_0;
      a_0 = y[(r1 << 2) + TWO];
      a_0 += b_y[TWO + 4] * A[r1 + 4];
      y[TWO + (r1 << 2)] = a_0;
      a_0 = y[(r1 << 2) + TWO];
      a_0 += b_y[TWO + 8] * A[r1 + 8];
      y[TWO + (r1 << 2)] = a_0;
      a_0 = y[(r1 << 2) + TWO];
      a_0 += b_y[TWO + 12] * A[r1 + 12];
      y[TWO + (r1 << 2)] = a_0;
    }
  }

  for (TWO = 0; TWO < 16; TWO++) {
    obj->P[TWO] = y[TWO] + obj->Q[TWO];
  }

  y_idx_1 = u2;
  for (TWO = 0; TWO < 8; TWO++) {
    a[TWO] = tmp_1[TWO];
  }

  for (TWO = 0; TWO < 16; TWO++) {
    b[TWO] = obj->P[TWO];
  }

  for (TWO = 0; TWO < 4; TWO++) {
    for (r1 = 0; r1 < 2; r1++) {
      b_y_0[r1 + (TWO << 1)] = 0.0;
      a22 = b_y_0[(TWO << 1) + r1];
      a22 += b[TWO << 2] * (real_T)a[r1];
      b_y_0[r1 + (TWO << 1)] = a22;
      a22 = b_y_0[(TWO << 1) + r1];
      a22 += b[(TWO << 2) + 1] * (real_T)a[r1 + 2];
      b_y_0[r1 + (TWO << 1)] = a22;
      a22 = b_y_0[(TWO << 1) + r1];
      a22 += b[(TWO << 2) + 2] * (real_T)a[r1 + 4];
      b_y_0[r1 + (TWO << 1)] = a22;
      a22 = b_y_0[(TWO << 1) + r1];
      a22 += b[(TWO << 2) + 3] * (real_T)a[r1 + 6];
      b_y_0[r1 + (TWO << 1)] = a22;
    }
  }

  for (TWO = 0; TWO < 8; TWO++) {
    a[TWO] = tmp_1[TWO];
  }

  for (TWO = 0; TWO < 2; TWO++) {
    for (r1 = 0; r1 < 2; r1++) {
      x_prev[TWO + (r1 << 1)] = 0.0;
      x_prev[TWO + (r1 << 1)] += b_y_0[TWO] * (real_T)a[r1];
      x_prev[TWO + (r1 << 1)] += b_y_0[TWO + 2] * (real_T)a[r1 + 2];
      x_prev[TWO + (r1 << 1)] += b_y_0[TWO + 4] * (real_T)a[r1 + 4];
      x_prev[TWO + (r1 << 1)] += b_y_0[TWO + 6] * (real_T)a[r1 + 6];
    }
  }

  a22 = x_prev[0];
  a22 += obj->R[0];
  x_prev[0] = a22;
  a22 = x_prev[1];
  a22 += obj->R[1];
  x_prev[1] = a22;
  a22 = x_prev[2];
  a22 += obj->R[2];
  x_prev[2] = a22;
  a22 = x_prev[3];
  a22 += obj->R[3];
  x_prev[3] = a22;
  for (TWO = 0; TWO < 16; TWO++) {
    b_y[TWO] = obj->P[TWO];
  }

  for (TWO = 0; TWO < 8; TWO++) {
    a[TWO] = tmp_1[TWO];
  }

  for (TWO = 0; TWO < 4; TWO++) {
    for (r1 = 0; r1 < 2; r1++) {
      b_y_0[TWO + (r1 << 2)] = 0.0;
      b_y_0[TWO + (r1 << 2)] += b_y[TWO] * (real_T)a[r1];
      b_y_0[TWO + (r1 << 2)] += b_y[TWO + 4] * (real_T)a[r1 + 2];
      b_y_0[TWO + (r1 << 2)] += b_y[TWO + 8] * (real_T)a[r1 + 4];
      b_y_0[TWO + (r1 << 2)] += b_y[TWO + 12] * (real_T)a[r1 + 6];
    }
  }

  TWO = 1;
  phase_zero2_end = x_prev[1];
  a22 = fabs(phase_zero2_end);
  ONE = a22;
  phase_zero2_end = x_prev[0];
  a22 = fabs(phase_zero2_end);
  if (ONE > a22) {
    r1 = 1;
    TWO = 0;
  } else {
    r1 = 0;
  }

  ONE = x_prev[TWO] / x_prev[r1];
  a22 = x_prev[TWO + 2] - x_prev[r1 + 2] * ONE;
  K_gain[r1 << 2] = b_y_0[0] / x_prev[r1];
  K_gain[TWO << 2] = (b_y_0[4] - K_gain[r1 << 2] * x_prev[r1 + 2]) / a22;
  K_gain[r1 << 2] -= K_gain[TWO << 2] * ONE;
  K_gain[(r1 << 2) + 1] = b_y_0[1] / x_prev[r1];
  K_gain[(TWO << 2) + 1] = (b_y_0[5] - K_gain[(r1 << 2) + 1] * x_prev[r1 + 2]) /
    a22;
  K_gain[(r1 << 2) + 1] -= K_gain[(TWO << 2) + 1] * ONE;
  K_gain[(r1 << 2) + 2] = b_y_0[2] / x_prev[r1];
  K_gain[(TWO << 2) + 2] = (b_y_0[6] - K_gain[(r1 << 2) + 2] * x_prev[r1 + 2]) /
    a22;
  K_gain[(r1 << 2) + 2] -= K_gain[(TWO << 2) + 2] * ONE;
  K_gain[(r1 << 2) + 3] = b_y_0[3] / x_prev[r1];
  K_gain[(TWO << 2) + 3] = (b_y_0[7] - K_gain[(r1 << 2) + 3] * x_prev[r1 + 2]) /
    a22;
  K_gain[(r1 << 2) + 3] -= K_gain[(TWO << 2) + 3] * ONE;
  a_0 = t_sine;
  a_0 -= x_pred_idx_0;
  t_sine = a_0;
  a_0 = y_idx_1;
  a_0 -= x_pred_idx_2;
  y_idx_1 = a_0;
  for (TWO = 0; TWO <= 2; TWO += 2) {
    /* MATLABSystem: '<Root>/MATLAB System' */
    tmp = _mm_loadu_pd(&K_gain[TWO]);
    tmp = _mm_mul_pd(tmp, _mm_set1_pd(t_sine));
    tmp = _mm_add_pd(tmp, _mm_set1_pd(0.0));
    tmp_0 = _mm_loadu_pd(&K_gain[TWO + 4]);
    tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(y_idx_1));
    tmp = _mm_add_pd(tmp_0, tmp);

    /* MATLABSystem: '<Root>/MATLAB System' */
    _mm_storeu_pd(&x_prev[TWO], tmp);
  }

  /* MATLABSystem: '<Root>/MATLAB System' */
  obj->xh[0] = x_pred_idx_0 + x_prev[0];
  obj->xh[1] = x_pred_idx_1 + x_prev[1];
  obj->xh[2] = x_pred_idx_2 + x_prev[2];
  obj->xh[3] = x_pred_idx_3 + x_prev[3];
  for (TWO = 0; TWO < 8; TWO++) {
    a[TWO] = tmp_1[TWO];
  }

  for (TWO = 0; TWO < 4; TWO++) {
    for (r1 = 0; r1 <= 2; r1 += 2) {
      _mm_storeu_pd(&y[r1 + (TWO << 2)], _mm_set1_pd(0.0));
      tmp = _mm_loadu_pd(&K_gain[r1]);
      tmp = _mm_mul_pd(_mm_set1_pd(a[TWO << 1]), tmp);
      tmp_0 = _mm_loadu_pd(&y[(TWO << 2) + r1]);
      tmp = _mm_add_pd(tmp, tmp_0);
      _mm_storeu_pd(&y[r1 + (TWO << 2)], tmp);
      tmp = _mm_loadu_pd(&K_gain[r1 + 4]);
      tmp = _mm_mul_pd(_mm_set1_pd(a[(TWO << 1) + 1]), tmp);
      tmp_0 = _mm_loadu_pd(&y[(TWO << 2) + r1]);
      tmp = _mm_add_pd(tmp, tmp_0);
      _mm_storeu_pd(&y[r1 + (TWO << 2)], tmp);
    }
  }

  for (TWO = 0; TWO < 16; TWO++) {
    A[TWO] = 0.0;
  }

  A[0] = 1.0;
  A[5] = 1.0;
  A[10] = 1.0;
  A[15] = 1.0;
  for (TWO = 0; TWO <= 14; TWO += 2) {
    /* MATLABSystem: '<Root>/MATLAB System' */
    tmp = _mm_loadu_pd(&A[TWO]);
    tmp_0 = _mm_loadu_pd(&y[TWO]);
    tmp = _mm_sub_pd(tmp, tmp_0);

    /* MATLABSystem: '<Root>/MATLAB System' */
    _mm_storeu_pd(&A[TWO], tmp);
  }

  /* MATLABSystem: '<Root>/MATLAB System' */
  for (TWO = 0; TWO < 16; TWO++) {
    b[TWO] = obj->P[TWO];
  }

  for (TWO = 0; TWO < 4; TWO++) {
    for (r1 = 0; r1 <= 2; r1 += 2) {
      _mm_storeu_pd(&y[r1 + (TWO << 2)], _mm_set1_pd(0.0));
      tmp = _mm_loadu_pd(&A[r1]);
      tmp = _mm_mul_pd(_mm_set1_pd(b[TWO << 2]), tmp);
      tmp_0 = _mm_loadu_pd(&y[(TWO << 2) + r1]);
      tmp = _mm_add_pd(tmp, tmp_0);
      _mm_storeu_pd(&y[r1 + (TWO << 2)], tmp);
      tmp = _mm_loadu_pd(&A[r1 + 4]);
      tmp = _mm_mul_pd(_mm_set1_pd(b[(TWO << 2) + 1]), tmp);
      tmp_0 = _mm_loadu_pd(&y[(TWO << 2) + r1]);
      tmp = _mm_add_pd(tmp, tmp_0);
      _mm_storeu_pd(&y[r1 + (TWO << 2)], tmp);
      tmp = _mm_loadu_pd(&A[r1 + 8]);
      tmp = _mm_mul_pd(_mm_set1_pd(b[(TWO << 2) + 2]), tmp);
      tmp_0 = _mm_loadu_pd(&y[(TWO << 2) + r1]);
      tmp = _mm_add_pd(tmp, tmp_0);
      _mm_storeu_pd(&y[r1 + (TWO << 2)], tmp);
      tmp = _mm_loadu_pd(&A[r1 + 12]);
      tmp = _mm_mul_pd(_mm_set1_pd(b[(TWO << 2) + 3]), tmp);
      tmp_0 = _mm_loadu_pd(&y[(TWO << 2) + r1]);
      tmp = _mm_add_pd(tmp, tmp_0);
      _mm_storeu_pd(&y[r1 + (TWO << 2)], tmp);
    }
  }

  for (TWO = 0; TWO < 16; TWO++) {
    obj->P[TWO] = y[TWO];
  }

  obj->xhat_prev[0] = obj->xh[0];
  obj->xhat_prev[1] = obj->xh[1];
  obj->xhat_prev[2] = obj->xh[2];
  obj->xhat_prev[3] = obj->xh[3];
  if (u0 < 5.0) {
    ONE = 0.0;
    L = 0.0;
    a22 = 0.0;
  } else if (u0 < 61.85) {
    t_sine = u0 - 5.0;
    phase_zero2_end = t_sine / 56.85;
    if (phase_zero2_end < 0.5) {
      amp = phase_zero2_end / 0.5 * 0.090000000000000011 + 0.05;
      phase_zero2_end = 0.11423973285781065 * t_sine;
      phase_zero2_end = sin(phase_zero2_end);
      phase_zero2_end = 0.83775804095727813 * t_sine - 0.2094395102393195 *
        phase_zero2_end / 0.11423973285781065;
      phase_zero2_end = sin(phase_zero2_end);
      x = 0.11423973285781065 * t_sine;
      x = sin(x);
      x = 0.83775804095727813 * t_sine - 0.2094395102393195 * x /
        0.11423973285781065;
      x = cos(x);
      x_0 = 0.11423973285781065 * t_sine;
      x_0 = cos(x_0);
      L = (0.83775804095727813 - 0.2094395102393195 * x_0) * (amp * x) +
        0.00316622691292876 * phase_zero2_end;
      phase_zero2_end = 6.2831853071795862 * t_sine / 55.0;
      phase_zero2_end = cos(phase_zero2_end);
      phase_zero2_end = 0.83775804095727813 - 3.1415926535897931 *
        phase_zero2_end / 15.0;
      a22 = phase_zero2_end * phase_zero2_end;
      phase_zero2_end = 6.2831853071795862 * t_sine / 55.0;
      phase_zero2_end = sin(phase_zero2_end);
      phase_zero2_end = 11.0 * phase_zero2_end / 6.0 - 12.566370614359172 *
        t_sine / 15.0;
      phase_zero2_end = cos(phase_zero2_end);
      x = 6.2831853071795862 * t_sine / 55.0;
      x = cos(x);
      x_0 = 6.2831853071795862 * t_sine / 55.0;
      x_0 = sin(x_0);
      x_0 = 11.0 * x_0 / 6.0 - 12.566370614359172 * t_sine / 15.0;
      x_0 = sin(x_0);
      rg = 6.2831853071795862 * t_sine / 55.0;
      rg = sin(rg);
      tau = 6.2831853071795862 * t_sine / 55.0;
      tau = sin(tau);
      tau = 11.0 * tau / 6.0 - 12.566370614359172 * t_sine / 15.0;
      tau = cos(tau);
      a22 = ((0.83775804095727813 - 3.1415926535897931 * x / 15.0) * (12.0 *
              phase_zero2_end) / 1895.0 + (6.0 * t_sine / 1895.0 + 0.05) * x_0 *
             a22) + (6.0 * t_sine / 1895.0 + 0.05) * (19.739208802178716 * rg *
        tau) / 825.0;
    } else {
      amp = 0.14;
      phase_zero2_end = 0.11423973285781065 * t_sine;
      phase_zero2_end = sin(phase_zero2_end);
      phase_zero2_end = 0.83775804095727813 * t_sine - 0.2094395102393195 *
        phase_zero2_end / 0.11423973285781065;
      phase_zero2_end = cos(phase_zero2_end);
      x = 0.11423973285781065 * t_sine;
      x = cos(x);
      L = (0.83775804095727813 - 0.2094395102393195 * x) * (0.14 *
        phase_zero2_end);
      phase_zero2_end = 6.2831853071795862 * t_sine / 55.0;
      phase_zero2_end = cos(phase_zero2_end);
      phase_zero2_end = 0.83775804095727813 - 3.1415926535897931 *
        phase_zero2_end / 15.0;
      a22 = phase_zero2_end * phase_zero2_end;
      phase_zero2_end = 6.2831853071795862 * t_sine / 55.0;
      phase_zero2_end = sin(phase_zero2_end);
      phase_zero2_end = 11.0 * phase_zero2_end / 6.0 - 12.566370614359172 *
        t_sine / 15.0;
      phase_zero2_end = sin(phase_zero2_end);
      x = 6.2831853071795862 * t_sine / 55.0;
      x = sin(x);
      x_0 = 6.2831853071795862 * t_sine / 55.0;
      x_0 = sin(x_0);
      x_0 = 11.0 * x_0 / 6.0 - 12.566370614359172 * t_sine / 15.0;
      x_0 = cos(x_0);
      a22 = 7.0 * phase_zero2_end * a22 / 50.0 + 69.0872308076255 * x * x_0 /
        20625.0;
    }

    phase_zero2_end = 0.11423973285781065 * t_sine;
    phase_zero2_end = sin(phase_zero2_end);
    phase_zero2_end = 0.83775804095727813 * t_sine - 0.2094395102393195 *
      phase_zero2_end / 0.11423973285781065;
    phase_zero2_end = sin(phase_zero2_end);
    ONE = amp * phase_zero2_end;
  } else if (u0 < 65.0) {
    ONE = 0.0;
    L = 0.0;
    a22 = 0.0;
  } else if (u0 < 85.0) {
    ONE = u0 - 65.0;
    a22 = ONE / 20.0;
    if (a22 < 0.5) {
      a22 = 0.05;
    } else {
      a22 = 0.1;
    }

    phase_zero2_end = 0.62831853071795862 * ONE;
    phase_zero2_end = sin(phase_zero2_end);
    if (phase_zero2_end < 0.0) {
      phase_zero2_end = -1.0;
    } else {
      phase_zero2_end = (phase_zero2_end > 0.0);
    }

    ONE = a22 * phase_zero2_end;
    L = 0.0;
    a22 = 0.0;
  } else {
    ONE = 0.0;
    L = 0.0;
    a22 = 0.0;
  }

  x_prev[0] = ONE;
  x_prev[1] = L;
  x_prev[2] = a22;
  varargin_3 = obj;
  x_pred_idx_0 = obj->xh[0];
  x_pred_idx_1 = obj->xh[1];
  x_pred_idx_2 = obj->xh[2];
  x_pred_idx_3 = obj->xh[3];
  t_sine = x_pred_idx_0;
  b_0 = x_pred_idx_1;
  amp = x_pred_idx_2;
  L = x_pred_idx_3;
  phase_zero2_end = x_prev[0];
  x = x_prev[1];
  x_0 = x_prev[2];
  ONE = L * L;
  a_0 = amp;
  a_0 = cos(a_0);
  a22 = a_0 * a_0;
  u_prev = L * L;
  rg = L * L;
  a_0 = amp;
  a_0 = sin(a_0);
  tau = amp;
  tau = cos(tau);
  b_x = amp;
  b_x = sin(b_x);
  b_x_0 = amp;
  b_x_0 = cos(b_x_0);
  b_x_1 = amp;
  b_x_1 = cos(b_x_1);
  b_x_2 = amp;
  b_x_2 = cos(b_x_2);
  amp = sin(amp);
  x_prev[0] = t_sine - phase_zero2_end;
  x_prev[1] = b_0 - x;
  x_prev[2] = ((5.0 * t_sine / 7.0 - 0.15196428571428572) * (64516.0 * ONE * a22)
               / 1.8105025E+7 + 124587.0 * a_0 / 297850.0) - x_0;
  x_prev[3] = (((((108077.0 * u_prev * b_x / 4.0E+8 + 108077.0 * L * b_x_0 /
                   1.0E+7) - 127.0 * t_sine * L * b_x_1 / 2500.0) + 127.0 * b_0 *
                 L * b_x_2 / 200000.0) - 127.0 * t_sine * rg * amp / 100000.0) +
               0.104353875) * (2.032E+7 * L * tau) / 5.069407E+6;
  y_idx_1 = -varargin_3->K[0];
  t_sine = -varargin_3->K[1];
  a22 = -varargin_3->K[2];
  a_0 = -varargin_3->K[3];
  y_idx_1 *= x_prev[0];
  y_idx_1 += t_sine * x_prev[1];
  y_idx_1 += a22 * x_prev[2];
  y_idx_1 += a_0 * x_prev[3];
  t_sine = x_pred_idx_0;
  amp = x_pred_idx_2;
  ONE = L * L;
  a22 = L * L;
  a_0 = amp;
  a_0 = sin(a_0);
  u_prev = rt_powd_snf(a_0, 3.0);
  rg = L * L;
  a_0 = amp;
  a_0 = cos(a_0);
  x_pred_idx_0 = a_0 * a_0;
  x_pred_idx_1 = rt_powd_snf(L, 4.0);
  a_0 = amp;
  a_0 = cos(a_0);
  x_pred_idx_2 = a_0 * a_0;
  x_pred_idx_3 = rt_powd_snf(L, 4.0);
  a_0 = amp;
  a_0 = cos(a_0);
  c = rt_powd_snf(a_0, 4.0);
  c_0 = rt_powd_snf(L, 3.0);
  c_1 = rt_powd_snf(L, 4.0);
  c_2 = rt_powd_snf(L, 4.0);
  c_3 = L * L;
  a_0 = amp;
  a_0 = cos(a_0);
  c_4 = a_0 * a_0;
  c_5 = L * L;
  a_0 = amp;
  a_0 = cos(a_0);
  c_6 = a_0 * a_0;
  c_7 = rt_powd_snf(L, 4.0);
  a_0 = amp;
  a_0 = cos(a_0);
  c_8 = a_0 * a_0;
  c_9 = rt_powd_snf(L, 4.0);
  a_0 = amp;
  a_0 = cos(a_0);
  c_a = rt_powd_snf(a_0, 4.0);
  c_b = rt_powd_snf(L, 3.0);
  c_c = rt_powd_snf(L, 3.0);
  c_d = L * L;
  c_e = L * L;
  a_0 = amp;
  a_0 = sin(a_0);
  tau = amp;
  tau = sin(tau);
  phase_zero2_end = 2.0 * amp;
  phase_zero2_end = sin(phase_zero2_end);
  b_x = amp;
  b_x = cos(b_x);
  x = 2.0 * amp;
  x = sin(x);
  x_0 = 2.0 * amp;
  x_0 = sin(x_0);
  b_x_0 = amp;
  b_x_0 = cos(b_x_0);
  b_x_1 = amp;
  b_x_1 = sin(b_x_1);
  b_x_2 = amp;
  b_x_2 = cos(b_x_2);
  b_x_3 = amp;
  b_x_3 = cos(b_x_3);
  b_x_4 = amp;
  b_x_4 = cos(b_x_4);
  amp = sin(amp);
  ONE = (((((((((((((((5.37476460632559E+14 * ONE * a_0 / 6.4E+17 +
                       2.5698887331649E+13 * y_idx_1 / 1.28E+16) -
                      1.710053628273E+12 * a22 * (tau - u_prev) / 8.0E+17) +
                     6.9581560143053E+13 * rg * x_pred_idx_0 / 1.0E+16) -
                    6.9581560143053E+13 * x_pred_idx_1 * x_pred_idx_2 / 1.6E+19)
                   + 2.21383089491E+11 * x_pred_idx_3 * c / 8.0E+19) +
                  6.9581560143053E+13 * c_0 * phase_zero2_end / 3.2E+17) -
                 8.1764465503E+10 * t_sine * c_1 / 8.0E+15) +
                5.37476460632559E+14 * L * b_x / 1.6E+16) + 6.9581560143053E+13 *
               c_2 / 3.2E+19) - 8.1764465503E+10 * t_sine * c_3 * c_4 / 2.5E+12)
             + 8.1764465503E+10 * b_0 * c_5 * c_6 / 1.0E+14) + 8.1764465503E+10 *
            t_sine * c_7 * c_8 / 4.0E+15) - 2.60144641E+8 * t_sine * c_9 * c_a /
           2.0E+16) - 8.1764465503E+10 * t_sine * c_b * x / 8.0E+13) +
         8.1764465503E+10 * b_0 * c_c * x_0 / 8.0E+15) * 4.0E+9 / ((((((324231.0
    * c_d * b_x_1 / 4.0E+8 + 108077.0 * L * b_x_2 / 5.0E+6) - 127.0 * t_sine * L
    * b_x_3 / 1250.0) + 127.0 * b_0 * L * b_x_4 / 100000.0) - 381.0 * t_sine *
    c_e * amp / 100000.0) + 0.104353875) * (1.931444067E+9 * b_x_0));
  a_0 = u2;
  a_0 = cos(a_0);
  ONE += a_0 * 0.25;

  /* E(0.1 + 0.5*(-p_ball+0.2)/0.4); */
  if (u2 < -1.1780972450961724) {
    u2 = (-1.1780972450961724 - u2) * 10.0;
    if (!(ONE >= u2)) {
      ONE = u2;
    }
  } else if (u2 > 1.1780972450961724) {
    u2 = (1.1780972450961724 - u2) * 10.0;
    if (!(ONE <= u2)) {
      ONE = u2;
    }
  }

  u2 = obj->max_V;
  if ((!(ONE <= u2)) && (!rtIsNaN(u2))) {
    ONE = u2;
  }

  u2 = -obj->max_V;
  if ((!(ONE >= u2)) && (!rtIsNaN(u2))) {
    ONE = u2;
  }

  u2 = obj->xh[0];
  a22 = obj->xh[1];
  phase_zero2_end = obj->xh[2];
  u_prev = obj->xh[3];
  obj->u = ONE;
  obj->t_prev = u0;

  /* MATLABSystem: '<Root>/MATLAB System' */
  /*  %% Sample Controller: Simple Proportional Controller */
  /*  t_prev = obj.t_prev; */
  /*  % Extract reference trajectory at the current timestep. */
  /*  [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t); */
  /*  % Decide desired servo angle based on simple proportional feedback. */
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
  /*   */
  /*  % Update class properties if necessary. */
  /*  obj.t_prev = t; */
  /*  obj.theta_d = theta_d; */
  simulink_experiment_debug_typ_B.MATLABSystem_o1 = ONE;

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_debug_typ_B.MATLABSystem_o2 = u2;

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_debug_typ_B.MATLABSystem_o3 = a22;

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_debug_typ_B.MATLABSystem_o4 = phase_zero2_end;

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_debug_typ_B.MATLABSystem_o5 = u_prev;

  /* Saturate: '<Root>/+//-10V' */
  u0 = simulink_experiment_debug_typ_B.MATLABSystem_o1;
  t_sine = simulink_experiment_debug_typ_P.u0V_LowerSat;
  u2 = simulink_experiment_debug_typ_P.u0V_UpperSat;
  if (u0 > u2) {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_debug_typ_B.u0V = u2;
  } else if (u0 < t_sine) {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_debug_typ_B.u0V = t_sine;
  } else {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_debug_typ_B.u0V = u0;
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
    phase_zero2_end = (simulink_experiment_debug_typ_B.Clock - 5.0) / 56.85;
    if (phase_zero2_end < 0.5) {
      amp = phase_zero2_end / 0.5 * 0.090000000000000011 + 0.05;
      simulink_experiment_debug_typ_B.v_ref = cos
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 -
         sin((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065)
         * 0.2094395102393195 / 0.11423973285781065) * amp *
        (0.83775804095727813 - cos((simulink_experiment_debug_typ_B.Clock - 5.0)
          * 0.11423973285781065) * 0.2094395102393195) + sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 -
         sin((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065)
         * 0.2094395102393195 / 0.11423973285781065) * 0.00316622691292876;
      u0 = 0.83775804095727813 - cos((simulink_experiment_debug_typ_B.Clock -
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
        5.0) * 6.0 / 1895.0 + 0.05) * (u0 * u0)) + cos(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * (sin((simulink_experiment_debug_typ_B.Clock
        - 5.0) * 6.2831853071795862 / 55.0) * 19.739208802178716) *
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.0 / 1895.0 + 0.05) /
        825.0;
    } else {
      amp = 0.14;
      simulink_experiment_debug_typ_B.v_ref = cos
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 -
         sin((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065)
         * 0.2094395102393195 / 0.11423973285781065) * 0.14 *
        (0.83775804095727813 - cos((simulink_experiment_debug_typ_B.Clock - 5.0)
          * 0.11423973285781065) * 0.2094395102393195);
      a_0 = 0.83775804095727813 - cos((simulink_experiment_debug_typ_B.Clock -
        5.0) * 6.2831853071795862 / 55.0) * 3.1415926535897931 / 15.0;
      simulink_experiment_debug_typ_B.a_ref = sin(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * 7.0 * (a_0 * a_0) / 50.0 + cos(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * (sin((simulink_experiment_debug_typ_B.Clock
        - 5.0) * 6.2831853071795862 / 55.0) * 69.0872308076255) / 20625.0;
    }

    simulink_experiment_debug_typ_B.p_ref = sin
      ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 - sin
       ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065) *
       0.2094395102393195 / 0.11423973285781065) * amp;
  } else if (simulink_experiment_debug_typ_B.Clock < 65.0) {
    simulink_experiment_debug_typ_B.p_ref = 0.0;
    simulink_experiment_debug_typ_B.v_ref = 0.0;
    simulink_experiment_debug_typ_B.a_ref = 0.0;
  } else if (simulink_experiment_debug_typ_B.Clock < 85.0) {
    if ((simulink_experiment_debug_typ_B.Clock - 65.0) / 20.0 < 0.5) {
      a22 = 0.05;
    } else {
      a22 = 0.1;
    }

    u0 = sin((simulink_experiment_debug_typ_B.Clock - 65.0) *
             0.62831853071795862);
    if (rtIsNaN(u0)) {
      a_0 = (rtNaN);
    } else if (u0 < 0.0) {
      a_0 = -1.0;
    } else {
      a_0 = (u0 > 0.0);
    }

    simulink_experiment_debug_typ_B.p_ref = a22 * a_0;
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
    studentControllerInterface_stud(&simulink_experiment_debug_ty_DW.obj);
    simulink_experiment_debug_ty_DW.objisempty = true;
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
  simulink_experiment_debug_ty_M->Sizes.checksums[0] = (3455148126U);
  simulink_experiment_debug_ty_M->Sizes.checksums[1] = (3727089829U);
  simulink_experiment_debug_ty_M->Sizes.checksums[2] = (212033803U);
  simulink_experiment_debug_ty_M->Sizes.checksums[3] = (3872749805U);

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
    simulink_experiment_debug_typ_B.MATLABSystem_o2 = 0.0;
    simulink_experiment_debug_typ_B.MATLABSystem_o3 = 0.0;
    simulink_experiment_debug_typ_B.MATLABSystem_o4 = 0.0;
    simulink_experiment_debug_typ_B.MATLABSystem_o5 = 0.0;
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
  simulink_experiment_debug_ty_M->Sizes.numBlocks = (33);/* Number of blocks */
  simulink_experiment_debug_ty_M->Sizes.numBlockIO = (23);/* Number of block outputs */
  simulink_experiment_debug_ty_M->Sizes.numBlockPrms = (88);/* Sum of parameter "widths" */
  return simulink_experiment_debug_ty_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
