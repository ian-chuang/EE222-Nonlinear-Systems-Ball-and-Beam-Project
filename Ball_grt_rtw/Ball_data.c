/*
 * Ball_data.c
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

/* Block parameters (default storage) */
P_Ball_T Ball_P = {
  /* Variable: K_bb
   * Referenced by: '<S1>/Model Gain  (m//s^2//rad)'
   */
  0.82,

  /* Expression: 0
   * Referenced by: '<S1>/vel to pos:  x'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S1>/SRV02: Vel to Pos'
   */
  0.0,

  /* Computed Parameter: ServoModel_A
   * Referenced by: '<S1>/Servo Model'
   */
  -10.0,

  /* Computed Parameter: ServoModel_C
   * Referenced by: '<S1>/Servo Model'
   */
  100.0,

  /* Expression: 0
   * Referenced by: '<S1>/acc to vel: x_dot'
   */
  0.0
};
