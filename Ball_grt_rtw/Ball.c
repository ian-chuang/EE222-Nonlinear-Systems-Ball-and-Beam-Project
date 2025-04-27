/*
 * Ball.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "Ball".
 *
 * Model version              : 1.0
 * Simulink Coder version : 24.2 (R2024b) 21-Jun-2024
 * C source code generated on : Sat Apr 26 19:58:23 2025
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "Ball.h"
#include <math.h>
#include "Ball_private.h"
#include <string.h>

/* Block signals (default storage) */
B_Ball_T Ball_B;

/* Continuous states */
X_Ball_T Ball_X;

/* Disabled State Vector */
XDis_Ball_T Ball_XDis;

/* External inputs (root inport signals with default storage) */
ExtU_Ball_T Ball_U;

/* External outputs (root outports fed by signals with default storage) */
ExtY_Ball_T Ball_Y;

/* Real-time model */
static RT_MODEL_Ball_T Ball_M_;
RT_MODEL_Ball_T *const Ball_M = &Ball_M_;

/*
 * This function updates continuous states using the ODE4 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE4_IntgData *id = (ODE4_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T *f3 = id->f[3];
  real_T temp;
  int_T i;
  int_T nXc = 4;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  Ball_derivatives();

  /* f1 = f(t + (h/2), y + (h/2)*f0) */
  temp = 0.5 * h;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f0[i]);
  }

  rtsiSetT(si, t + temp);
  rtsiSetdX(si, f1);
  Ball_step();
  Ball_derivatives();

  /* f2 = f(t + (h/2), y + (h/2)*f1) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f1[i]);
  }

  rtsiSetdX(si, f2);
  Ball_step();
  Ball_derivatives();

  /* f3 = f(t + h, y + h*f2) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (h*f2[i]);
  }

  rtsiSetT(si, tnew);
  rtsiSetdX(si, f3);
  Ball_step();
  Ball_derivatives();

  /* tnew = t + h
     ynew = y + (h/6)*(f0 + 2*f1 + 2*f2 + 2*f3) */
  temp = h / 6.0;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + temp*(f0[i] + 2.0*f1[i] + 2.0*f2[i] + f3[i]);
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model step function */
void Ball_step(void)
{
  if (rtmIsMajorTimeStep(Ball_M)) {
    /* set solver stop time */
    if (!(Ball_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&Ball_M->solverInfo, ((Ball_M->Timing.clockTickH0 +
        1) * Ball_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&Ball_M->solverInfo, ((Ball_M->Timing.clockTick0 + 1)
        * Ball_M->Timing.stepSize0 + Ball_M->Timing.clockTickH0 *
        Ball_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(Ball_M)) {
    Ball_M->Timing.t[0] = rtsiGetT(&Ball_M->solverInfo);
  }

  /* Outport: '<Root>/x (m)' incorporates:
   *  Integrator: '<S1>/vel to pos:  x'
   */
  Ball_Y.xm = Ball_X.veltoposx_CSTATE;

  /* Outport: '<Root>/theta_l (rad)' incorporates:
   *  Integrator: '<S1>/SRV02: Vel to Pos'
   */
  Ball_Y.theta_lrad = Ball_X.SRV02VeltoPos_CSTATE;

  /* Gain: '<S1>/Model Gain  (m//s^2//rad)' incorporates:
   *  Integrator: '<S1>/SRV02: Vel to Pos'
   *  Trigonometry: '<S1>/Trigonometric Function'
   */
  Ball_B.ModelGainms2rad = Ball_P.K_bb * sin(Ball_X.SRV02VeltoPos_CSTATE);

  /* TransferFcn: '<S1>/Servo Model' */
  Ball_B.ServoModel = Ball_P.ServoModel_C * Ball_X.ServoModel_CSTATE;

  /* Integrator: '<S1>/acc to vel: x_dot' */
  Ball_B.acctovelx_dot = Ball_X.acctovelx_dot_CSTATE;
  if (rtmIsMajorTimeStep(Ball_M)) {
    /* Matfile logging */
    rt_UpdateTXYLogVars(Ball_M->rtwLogInfo, (Ball_M->Timing.t));
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(Ball_M)) {
    /* signal main to stop simulation */
    {                                  /* Sample time: [0.0s, 0.0s] */
      if ((rtmGetTFinal(Ball_M)!=-1) &&
          !((rtmGetTFinal(Ball_M)-(((Ball_M->Timing.clockTick1+
               Ball_M->Timing.clockTickH1* 4294967296.0)) * 0.002)) >
            (((Ball_M->Timing.clockTick1+Ball_M->Timing.clockTickH1*
               4294967296.0)) * 0.002) * (DBL_EPSILON))) {
        rtmSetErrorStatus(Ball_M, "Simulation finished");
      }
    }

    rt_ertODEUpdateContinuousStates(&Ball_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++Ball_M->Timing.clockTick0)) {
      ++Ball_M->Timing.clockTickH0;
    }

    Ball_M->Timing.t[0] = rtsiGetSolverStopTime(&Ball_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.002s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.002, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       * Timer of this task consists of two 32 bit unsigned integers.
       * The two integers represent the low bits Timing.clockTick1 and the high bits
       * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
       */
      Ball_M->Timing.clockTick1++;
      if (!Ball_M->Timing.clockTick1) {
        Ball_M->Timing.clockTickH1++;
      }
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void Ball_derivatives(void)
{
  XDot_Ball_T *_rtXdot;
  _rtXdot = ((XDot_Ball_T *) Ball_M->derivs);

  /* Derivatives for Integrator: '<S1>/vel to pos:  x' */
  _rtXdot->veltoposx_CSTATE = Ball_B.acctovelx_dot;

  /* Derivatives for Integrator: '<S1>/SRV02: Vel to Pos' */
  _rtXdot->SRV02VeltoPos_CSTATE = Ball_B.ServoModel;

  /* Derivatives for TransferFcn: '<S1>/Servo Model' incorporates:
   *  Inport: '<Root>/Vm (V)'
   */
  _rtXdot->ServoModel_CSTATE = Ball_P.ServoModel_A * Ball_X.ServoModel_CSTATE;
  _rtXdot->ServoModel_CSTATE += Ball_U.VmV;

  /* Derivatives for Integrator: '<S1>/acc to vel: x_dot' */
  _rtXdot->acctovelx_dot_CSTATE = Ball_B.ModelGainms2rad;
}

/* Model initialize function */
void Ball_initialize(void)
{
  /* Registration code */

  /* initialize real-time model */
  (void) memset((void *)Ball_M, 0,
                sizeof(RT_MODEL_Ball_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&Ball_M->solverInfo, &Ball_M->Timing.simTimeStep);
    rtsiSetTPtr(&Ball_M->solverInfo, &rtmGetTPtr(Ball_M));
    rtsiSetStepSizePtr(&Ball_M->solverInfo, &Ball_M->Timing.stepSize0);
    rtsiSetdXPtr(&Ball_M->solverInfo, &Ball_M->derivs);
    rtsiSetContStatesPtr(&Ball_M->solverInfo, (real_T **) &Ball_M->contStates);
    rtsiSetNumContStatesPtr(&Ball_M->solverInfo, &Ball_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&Ball_M->solverInfo,
      &Ball_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&Ball_M->solverInfo,
      &Ball_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&Ball_M->solverInfo,
      &Ball_M->periodicContStateRanges);
    rtsiSetContStateDisabledPtr(&Ball_M->solverInfo, (boolean_T**)
      &Ball_M->contStateDisabled);
    rtsiSetErrorStatusPtr(&Ball_M->solverInfo, (&rtmGetErrorStatus(Ball_M)));
    rtsiSetRTModelPtr(&Ball_M->solverInfo, Ball_M);
  }

  rtsiSetSimTimeStep(&Ball_M->solverInfo, MAJOR_TIME_STEP);
  rtsiSetIsMinorTimeStepWithModeChange(&Ball_M->solverInfo, false);
  rtsiSetIsContModeFrozen(&Ball_M->solverInfo, false);
  Ball_M->intgData.y = Ball_M->odeY;
  Ball_M->intgData.f[0] = Ball_M->odeF[0];
  Ball_M->intgData.f[1] = Ball_M->odeF[1];
  Ball_M->intgData.f[2] = Ball_M->odeF[2];
  Ball_M->intgData.f[3] = Ball_M->odeF[3];
  Ball_M->contStates = ((X_Ball_T *) &Ball_X);
  Ball_M->contStateDisabled = ((XDis_Ball_T *) &Ball_XDis);
  Ball_M->Timing.tStart = (0.0);
  rtsiSetSolverData(&Ball_M->solverInfo, (void *)&Ball_M->intgData);
  rtsiSetSolverName(&Ball_M->solverInfo,"ode4");
  rtmSetTPtr(Ball_M, &Ball_M->Timing.tArray[0]);
  rtmSetTFinal(Ball_M, 4.0);
  Ball_M->Timing.stepSize0 = 0.002;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = (NULL);
    Ball_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(Ball_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(Ball_M->rtwLogInfo, (NULL));
    rtliSetLogT(Ball_M->rtwLogInfo, "tout");
    rtliSetLogX(Ball_M->rtwLogInfo, "");
    rtliSetLogXFinal(Ball_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(Ball_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(Ball_M->rtwLogInfo, 4);
    rtliSetLogMaxRows(Ball_M->rtwLogInfo, 0);
    rtliSetLogDecimation(Ball_M->rtwLogInfo, 1);
    rtliSetLogY(Ball_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(Ball_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(Ball_M->rtwLogInfo, (NULL));
  }

  /* block I/O */
  (void) memset(((void *) &Ball_B), 0,
                sizeof(B_Ball_T));

  /* states (continuous) */
  {
    (void) memset((void *)&Ball_X, 0,
                  sizeof(X_Ball_T));
  }

  /* disabled states */
  {
    (void) memset((void *)&Ball_XDis, 0,
                  sizeof(XDis_Ball_T));
  }

  /* external inputs */
  Ball_U.VmV = 0.0;

  /* external outputs */
  (void)memset(&Ball_Y, 0, sizeof(ExtY_Ball_T));

  /* Matfile logging */
  rt_StartDataLoggingWithStartTime(Ball_M->rtwLogInfo, 0.0, rtmGetTFinal(Ball_M),
    Ball_M->Timing.stepSize0, (&rtmGetErrorStatus(Ball_M)));

  /* InitializeConditions for Integrator: '<S1>/vel to pos:  x' */
  Ball_X.veltoposx_CSTATE = Ball_P.veltoposx_IC;

  /* InitializeConditions for Integrator: '<S1>/SRV02: Vel to Pos' */
  Ball_X.SRV02VeltoPos_CSTATE = Ball_P.SRV02VeltoPos_IC;

  /* InitializeConditions for TransferFcn: '<S1>/Servo Model' */
  Ball_X.ServoModel_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S1>/acc to vel: x_dot' */
  Ball_X.acctovelx_dot_CSTATE = Ball_P.acctovelx_dot_IC;
}

/* Model terminate function */
void Ball_terminate(void)
{
  /* (no terminate code required) */
}
