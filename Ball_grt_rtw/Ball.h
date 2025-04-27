/*
 * Ball.h
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

#ifndef Ball_h_
#define Ball_h_
#ifndef Ball_COMMON_INCLUDES_
#define Ball_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "rt_logging.h"
#include "rt_nonfinite.h"
#include "math.h"
#endif                                 /* Ball_COMMON_INCLUDES_ */

#include "Ball_types.h"
#include <float.h>
#include <string.h>
#include <stddef.h>

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
#define rtmGetContStateDisabled(rtm)   ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
#define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
#define rtmGetContStates(rtm)          ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
#define rtmSetContStates(rtm, val)     ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
#define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
#define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetFinalTime
#define rtmGetFinalTime(rtm)           ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetIntgData
#define rtmGetIntgData(rtm)            ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
#define rtmSetIntgData(rtm, val)       ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
#define rtmGetOdeF(rtm)                ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
#define rtmSetOdeF(rtm, val)           ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
#define rtmGetOdeY(rtm)                ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
#define rtmSetOdeY(rtm, val)           ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
#define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
#define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
#define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
#define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetRTWLogInfo
#define rtmGetRTWLogInfo(rtm)          ((rtm)->rtwLogInfo)
#endif

#ifndef rtmGetZCCacheNeedsReset
#define rtmGetZCCacheNeedsReset(rtm)   ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
#define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
#define rtmGetdX(rtm)                  ((rtm)->derivs)
#endif

#ifndef rtmSetdX
#define rtmSetdX(rtm, val)             ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTFinal
#define rtmGetTFinal(rtm)              ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

#ifndef rtmGetTStart
#define rtmGetTStart(rtm)              ((rtm)->Timing.tStart)
#endif

/* Block signals (default storage) */
typedef struct {
  real_T ModelGainms2rad;              /* '<S1>/Model Gain  (m//s^2//rad)' */
  real_T ServoModel;                   /* '<S1>/Servo Model' */
  real_T acctovelx_dot;                /* '<S1>/acc to vel: x_dot' */
} B_Ball_T;

/* Continuous states (default storage) */
typedef struct {
  real_T veltoposx_CSTATE;             /* '<S1>/vel to pos:  x' */
  real_T SRV02VeltoPos_CSTATE;         /* '<S1>/SRV02: Vel to Pos' */
  real_T ServoModel_CSTATE;            /* '<S1>/Servo Model' */
  real_T acctovelx_dot_CSTATE;         /* '<S1>/acc to vel: x_dot' */
} X_Ball_T;

/* State derivatives (default storage) */
typedef struct {
  real_T veltoposx_CSTATE;             /* '<S1>/vel to pos:  x' */
  real_T SRV02VeltoPos_CSTATE;         /* '<S1>/SRV02: Vel to Pos' */
  real_T ServoModel_CSTATE;            /* '<S1>/Servo Model' */
  real_T acctovelx_dot_CSTATE;         /* '<S1>/acc to vel: x_dot' */
} XDot_Ball_T;

/* State disabled  */
typedef struct {
  boolean_T veltoposx_CSTATE;          /* '<S1>/vel to pos:  x' */
  boolean_T SRV02VeltoPos_CSTATE;      /* '<S1>/SRV02: Vel to Pos' */
  boolean_T ServoModel_CSTATE;         /* '<S1>/Servo Model' */
  boolean_T acctovelx_dot_CSTATE;      /* '<S1>/acc to vel: x_dot' */
} XDis_Ball_T;

#ifndef ODE4_INTG
#define ODE4_INTG

/* ODE4 Integration Data */
typedef struct {
  real_T *y;                           /* output */
  real_T *f[4];                        /* derivatives */
} ODE4_IntgData;

#endif

/* External inputs (root inport signals with default storage) */
typedef struct {
  real_T VmV;                          /* '<Root>/Vm (V)' */
} ExtU_Ball_T;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real_T xm;                           /* '<Root>/x (m)' */
  real_T theta_lrad;                   /* '<Root>/theta_l (rad)' */
} ExtY_Ball_T;

/* Parameters (default storage) */
struct P_Ball_T_ {
  real_T K_bb;                         /* Variable: K_bb
                                        * Referenced by: '<S1>/Model Gain  (m//s^2//rad)'
                                        */
  real_T veltoposx_IC;                 /* Expression: 0
                                        * Referenced by: '<S1>/vel to pos:  x'
                                        */
  real_T SRV02VeltoPos_IC;             /* Expression: 0
                                        * Referenced by: '<S1>/SRV02: Vel to Pos'
                                        */
  real_T ServoModel_A;                 /* Computed Parameter: ServoModel_A
                                        * Referenced by: '<S1>/Servo Model'
                                        */
  real_T ServoModel_C;                 /* Computed Parameter: ServoModel_C
                                        * Referenced by: '<S1>/Servo Model'
                                        */
  real_T acctovelx_dot_IC;             /* Expression: 0
                                        * Referenced by: '<S1>/acc to vel: x_dot'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_Ball_T {
  const char_T *errorStatus;
  RTWLogInfo *rtwLogInfo;
  RTWSolverInfo solverInfo;
  X_Ball_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  XDis_Ball_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[4];
  real_T odeF[4][4];
  ODE4_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    time_T tStart;
    time_T tFinal;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Block parameters (default storage) */
extern P_Ball_T Ball_P;

/* Block signals (default storage) */
extern B_Ball_T Ball_B;

/* Continuous states (default storage) */
extern X_Ball_T Ball_X;

/* Disabled states (default storage) */
extern XDis_Ball_T Ball_XDis;

/* External inputs (root inport signals with default storage) */
extern ExtU_Ball_T Ball_U;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY_Ball_T Ball_Y;

/* Model entry point functions */
extern void Ball_initialize(void);
extern void Ball_step(void);
extern void Ball_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Ball_T *const Ball_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Note that this particular code originates from a subsystem build,
 * and has its own system numbers different from the parent model.
 * Refer to the system hierarchy for this subsystem below, and use the
 * MATLAB hilite_system command to trace the generated code back
 * to the parent model.  For example,
 *
 * hilite_system('simulink_simulation/Ball and Beam Model')    - opens subsystem simulink_simulation/Ball and Beam Model
 * hilite_system('simulink_simulation/Ball and Beam Model/Kp') - opens and selects block Kp
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'simulink_simulation'
 * '<S1>'   : 'simulink_simulation/Ball and Beam Model'
 */
#endif                                 /* Ball_h_ */
