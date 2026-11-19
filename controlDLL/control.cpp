// controlDLL.cpp : Defines the entry point for the DLL application.

#include "servo.h"
#include "param.h"
#include "control.h"
// #include "UiAgent.h"
#include "PrVector.h"
#include "PrMatrix.h"
#include "Utils.h" // misc. utility functions, such as toRad, toDeg, etc.
#include <math.h>
#include <algorithm>
using std::min;
using std::max;

void PrintDebug(GlobalVariables& gv);

// *******************************************************************
// Cubic spline trajectories in joint space (assignment 2, part B)
// *******************************************************************

// Parameters of one cubic spline per joint:
//   q_d(t) = a0 + a1 (t - t0) + a2 (t - t0)^2 + a3 (t - t0)^3,  t0 <= t <= t0 + tf
// t0 is the start time (simulator time), tf the duration of the trajectory.
struct CubicSpline {
  double t0, tf;
  PrVector a0, a1, a2, a3;
};

// Global, so that the init function can compute the spline once and the
// control function can follow it in every servo cycle.
CubicSpline spline;

// Computes the duration tf of a spline from the current joint angles gv.q to
// the desired joint angles gv.qd, starting and ending at rest.
// For such a spline the largest velocity occurs at t = tf/2 and the largest
// acceleration at t = 0:
//   |dq|max  = 1.5 * |qd - q| / tf      <= dqmax   =>  tf >= 1.5 * |qd - q| / dqmax
//   |ddq|max = 6 * |qd - q| / tf^2      <= ddqmax  =>  tf >= sqrt(6 * |qd - q| / ddqmax)
// tf is the largest of these values over all joints, so that every joint
// respects its limits and all joints reach their goal at the same time.
double computeTf(GlobalVariables& gv)
{
  double tf = 0.0;
  for (int i = 0; i < gv.dof; i++) {
    double distance = fabs(gv.qd[i] - gv.q[i]);
    double tf_velocity = 1.5 * distance / gv.dqmax[i];
    double tf_acceleration = sqrt(6.0 * distance / gv.ddqmax[i]);
    tf = max(tf, max(tf_velocity, tf_acceleration));
  }
  return tf;
}

// Computes the spline parameters from gv.q to gv.qd with zero velocity at the
// start and at the end (four conditions for the four unknowns of each joint):
//   a0 = q_start, a1 = 0, a2 = 3 (qd - q_start) / tf^2, a3 = -2 (qd - q_start) / tf^3
void computeSpline(GlobalVariables& gv)
{
  spline.t0 = gv.curTime;
  spline.tf = computeTf(gv);

  PrVector distance = gv.qd - gv.q;
  spline.a0 = gv.q;
  spline.a1 = PrVector(gv.dof);
  if (spline.tf > 0.0) {
    spline.a2 = distance * (3.0 / (spline.tf * spline.tf));
    spline.a3 = distance * (-2.0 / (spline.tf * spline.tf * spline.tf));
  } else {
    // already at the goal: the spline stays constant
    spline.a2 = PrVector(gv.dof);
    spline.a3 = PrVector(gv.dof);
  }
}

// Follows the spline with a joint space PD controller with gravity
// compensation. The desired position and velocity come from the spline:
//   tau = -kp (q - q_d) - kv (dq - dq_d) + G
// After tf the trajectory is finished and the robot switches to float control.
void followSpline(GlobalVariables& gv)
{
  double t = gv.curTime - spline.t0;
  if (t > spline.tf) {
    floatControl(gv);
    return;
  }

  gv.qd = spline.a0 + spline.a1 * t + spline.a2 * (t * t) + spline.a3 * (t * t * t);
  gv.dqd = spline.a1 + spline.a2 * (2.0 * t) + spline.a3 * (3.0 * t * t);

  gv.tau = -gv.kp * (gv.q - gv.qd) - gv.kv * (gv.dq - gv.dqd) + gv.G;
}

// *******************************************************************
// Circle in operational space (assignment 2, part C)
// *******************************************************************

const double PI = 3.14159265358979;

// Circle in the simulator frame (x right, y up): it starts at (0.8, 0.35) m,
// the point reached by proj1Control, and goes around the center (0.6, 0.35) m.
// One full circle takes 5 s.
const double CIRCLE_CENTER_X = 0.6;
const double CIRCLE_CENTER_Y = 0.35;
const double CIRCLE_RADIUS = 0.2;
const double CIRCLE_ANGULAR_VELOCITY = 2.0 * PI / 5.0;

// proj3: three full circles, the angular acceleration along the circle is
// limited to 2 pi / 25 s^-2
const double CIRCLE_ANGULAR_ACCELERATION = 2.0 * PI / 25.0;
const double CIRCLE_TOTAL_ANGLE = 3.0 * 2.0 * PI;

// Simulator time at which the circle was started
double circle_start_time = 0.0;

// Sets the desired end effector position and velocity for the circle angle
// beta and its angular velocity dbeta. beta = 0 is the start point, positive
// beta turns around the z axis of the simulator frame:
//   x_d  = center + r (cos(beta), sin(beta)),        alpha_d = 0 (upright)
//   dx_d = r * dbeta * (-sin(beta), cos(beta)),     dalpha_d = 0
void setCircleTarget(GlobalVariables& gv, double beta, double dbeta)
{
  gv.xd[0] = CIRCLE_CENTER_X + CIRCLE_RADIUS * cos(beta);
  gv.xd[1] = CIRCLE_CENTER_Y + CIRCLE_RADIUS * sin(beta);
  gv.xd[2] = 0.0;

  gv.dxd[0] = -CIRCLE_RADIUS * dbeta * sin(beta);
  gv.dxd[1] = CIRCLE_RADIUS * dbeta * cos(beta);
  gv.dxd[2] = 0.0;
}

// PD controller in operational space. The force on the end effector pulls it
// towards the desired position and velocity:
//   F = -kp (x - x_d) - kv (dx - dx_d)
// The Jacobian transpose of the simulator maps this force to joint torques,
// gravity is compensated with the simulator's gravity vector:
//   tau = J^T F + G
void operationalSpaceControl(GlobalVariables& gv)
{
  PrVector F = -gv.kp * (gv.x - gv.xd) - gv.kv * (gv.dx - gv.dxd);
  gv.tau = gv.Jtranspose * F + gv.G;
}

// *******************************************************************
// Bonus experiments (assignment 2, part C5)
// *******************************************************************

// (a) Other trajectories. They start at the same point as the circle (the end
// point of proj1) and move with the same path speed as the circle,
// r * 2 pi / 5 s = 0.25 m/s, so that the results can be compared.
const double PATH_SPEED = CIRCLE_RADIUS * CIRCLE_ANGULAR_VELOCITY;
const double LINE_LENGTH = 0.4;
const double SQUARE_SIDE = 0.2;

// Simulator time at which a bonus experiment was started
double bonus_start_time = 0.0;

// Line along x from (0.8, 0.35) m to (0.4, 0.35) m and back, repeated.
// The speed is constant, so the velocity jumps from -0.25 m/s to +0.25 m/s at
// both ends of the line.
void setLineTarget(GlobalVariables& gv, double t)
{
  double start_x = CIRCLE_CENTER_X + CIRCLE_RADIUS;
  // distance travelled within the current back-and-forth cycle
  double s = fmod(PATH_SPEED * t, 2.0 * LINE_LENGTH);
  if (s < LINE_LENGTH) {
    gv.xd[0] = start_x - s;
    gv.dxd[0] = -PATH_SPEED;
  } else {
    gv.xd[0] = start_x - (2.0 * LINE_LENGTH - s);
    gv.dxd[0] = PATH_SPEED;
  }
  gv.xd[1] = CIRCLE_CENTER_Y;
  gv.dxd[1] = 0.0;
  gv.xd[2] = 0.0;
  gv.dxd[2] = 0.0;
}

// Square with side 0.2 m, counterclockwise like the circle: up, left, down,
// right. At every corner the direction of the velocity changes by 90 degrees.
void setSquareTarget(GlobalVariables& gv, double t)
{
  const double corner_x[4] = {0.8, 0.8, 0.6, 0.6};
  const double corner_y[4] = {0.35, 0.55, 0.55, 0.35};

  // distance travelled within the current lap, and the side we are on
  double s = fmod(PATH_SPEED * t, 4.0 * SQUARE_SIDE);
  int side = (int)(s / SQUARE_SIDE);
  if (side > 3) {
    side = 3;
  }
  int next = (side + 1) % 4;
  double along = s - side * SQUARE_SIDE;

  // unit vector along the current side
  double dir_x = (corner_x[next] - corner_x[side]) / SQUARE_SIDE;
  double dir_y = (corner_y[next] - corner_y[side]) / SQUARE_SIDE;

  gv.xd[0] = corner_x[side] + along * dir_x;
  gv.xd[1] = corner_y[side] + along * dir_y;
  gv.xd[2] = 0.0;
  gv.dxd[0] = PATH_SPEED * dir_x;
  gv.dxd[1] = PATH_SPEED * dir_y;
  gv.dxd[2] = 0.0;
}

// (b) Errors and noise in the commanded torques, tested on the circle of proj2.
// Noise: every servo cycle a random torque of up to 10 % of tau_max is added
// to every joint.
// Model error: only 80 % of the gravity vector is compensated.
const double NOISE_RATIO = 0.1;
const double GRAVITY_COMPENSATION_RATIO = 0.8;

// Simple pseudo random number generator (linear congruential generator), so no
// additional include is needed. The fixed start value makes every run use the
// same noise.
unsigned int noise_state = 12345;

// returns a pseudo random number between -1 and 1
double randomUnit()
{
  noise_state = noise_state * 1103515245u + 12345u;
  return ((noise_state >> 16) & 0x7FFF) / 16383.5 - 1.0;
}

// circle of proj2 with constant angular velocity, started at bonus_start_time
void bonusCircle(GlobalVariables& gv)
{
  double beta = CIRCLE_ANGULAR_VELOCITY * (gv.curTime - bonus_start_time);
  setCircleTarget(gv, beta, CIRCLE_ANGULAR_VELOCITY);
  operationalSpaceControl(gv);
}

// *******************************************************************
// Initialization functions
// *******************************************************************

void InitControl(GlobalVariables& gv)
{
  // This code runs before the first servo loop
}

void PreprocessControl(GlobalVariables& gv)
{
  // This code runs on every servo loop, just before the control law

  // Assignment 2 uses the gravity vector gv.G computed by the simulator,
  // so it is not overwritten here (our model from assignment 1 is not used).
}

void PostprocessControl(GlobalVariables& gv)
{
  // This code runs on every servo loop, just after the control law
}

void initFloatControl(GlobalVariables& gv)
{
  // Control Initialization Code Here
}

void initOpenControl(GlobalVariables& gv)
{
  // Control Initialization Code Here
}

void initNjholdControl(GlobalVariables& gv)
{
  // Control Initialization Code Here
}

void initJholdControl(GlobalVariables& gv)
{
  // Control Initialization Code Here
}

void initNjmoveControl(GlobalVariables& gv)
{
  // Control Initialization Code Here
}

void initJmoveControl(GlobalVariables& gv)
{
  // Control Initialization Code Here
}

void initNjgotoControl(GlobalVariables& gv)
{
  // Control Initialization Code Here
}

void initJgotoControl(GlobalVariables& gv)
{
  // Control Initialization Code Here
}

void initNjtrackControl(GlobalVariables& gv)
{
  // gv.qd holds the goal configuration entered in the GUI
  computeSpline(gv);
}

void initJtrackControl(GlobalVariables& gv)
{
  // Control Initialization Code Here
}

void initNxtrackControl(GlobalVariables& gv)
{
  // bonus (b): circle with torque noise, same noise sequence in every run
  bonus_start_time = gv.curTime;
  noise_state = 12345;
}

void initXtrackControl(GlobalVariables& gv)
{
  // bonus (a): square
  bonus_start_time = gv.curTime;
}

void initNholdControl(GlobalVariables& gv)
{
  // Control Initialization Code Here
}

void initHoldControl(GlobalVariables& gv)
{
  // Control Initialization Code Here
}

void initNgotoControl(GlobalVariables& gv)
{
  // Control Initialization Code Here
}

void initGotoControl(GlobalVariables& gv)
{
  // Control Initialization Code Here
}

void initNtrackControl(GlobalVariables& gv)
{
  // bonus (b): circle with an error in the gravity compensation
  bonus_start_time = gv.curTime;
}

void initTrackControl(GlobalVariables& gv)
{
  // Control Initialization Code Here
}

void initPfmoveControl(GlobalVariables& gv)
{
  // Control Initialization Code Here
}

void initLineControl(GlobalVariables& gv)
{
  // bonus (a): line back and forth
  bonus_start_time = gv.curTime;
}

void initProj1Control(GlobalVariables& gv)
{
  // fixed goal configuration (5.5, 55.4, -60.8) degrees, then the same cubic
  // spline as in njtrack, so all joints reach the goal at the same time
  gv.qd[0] = toRad(5.5);
  gv.qd[1] = toRad(55.4);
  gv.qd[2] = toRad(-60.8);
  computeSpline(gv);
}

void initProj2Control(GlobalVariables& gv)
{
  // the circle starts at beta = 0 when proj2 is started
  circle_start_time = gv.curTime;
}

void initProj3Control(GlobalVariables& gv)
{
  // the circle starts at beta = 0 when proj3 is started
  circle_start_time = gv.curTime;
}


// *******************************************************************
// Control laws
// *******************************************************************

void noControl(GlobalVariables& gv)
{
}

void floatControl(GlobalVariables& gv)
{
  gv.tau = gv.G;
  // this only works on the real robot unless the function is changed to use cout
  // the handed in solution must not contain any printouts
  // PrintDebug(gv);
}

void openControl(GlobalVariables& gv)
{
  floatControl(gv);  // Remove this line when you implement this controller
}

void njholdControl(GlobalVariables& gv)
{
  floatControl(gv);  // Remove this line when you implement this controller
}

void jholdControl(GlobalVariables& gv)
{
  floatControl(gv);  // Remove this line when you implement this controller
}

void njmoveControl(GlobalVariables& gv)
{
  gv.tau = gv.kp * (gv.qd - gv.q);
}

void jmoveControl(GlobalVariables& gv)
{
  floatControl(gv);  // Remove this line when you implement this controller
}

void njgotoControl(GlobalVariables& gv)
{
  gv.tau = gv.kp * (gv.qd - gv.q) + gv.G;
}

void jgotoControl(GlobalVariables& gv)
{
  gv.tau = gv.kp * (gv.qd - gv.q) - gv.kv * gv.dq + gv.G;
}

void njtrackControl(GlobalVariables& gv)
{
  followSpline(gv);
}

void jtrackControl(GlobalVariables& gv)
{
  floatControl(gv);  // Remove this line when you implement this controller
}

// bonus (b): circle of proj2, a random torque is added to every joint
void nxtrackControl(GlobalVariables& gv)
{
  bonusCircle(gv);
  for (int i = 0; i < gv.dof; i++) {
    gv.tau[i] += NOISE_RATIO * gv.taumax[i] * randomUnit();
  }
}

// bonus (a): square in operational space
void xtrackControl(GlobalVariables& gv)
{
  setSquareTarget(gv, gv.curTime - bonus_start_time);
  operationalSpaceControl(gv);
}

void nholdControl(GlobalVariables& gv)
{
  floatControl(gv);  // Remove this line when you implement this controller
}

void holdControl(GlobalVariables& gv)
{
  floatControl(gv);  // Remove this line when you implement this controller
}

void ngotoControl(GlobalVariables& gv)
{
  floatControl(gv);  // Remove this line when you implement this controller
}

void gotoControl(GlobalVariables& gv)
{
  floatControl(gv);  // Remove this line when you implement this controller
}

// bonus (b): circle of proj2, but only 80 % of the gravity vector is
// compensated (tau = J^T F + 0.8 G)
void ntrackControl(GlobalVariables& gv)
{
  bonusCircle(gv);
  gv.tau = gv.tau - gv.G * (1.0 - GRAVITY_COMPENSATION_RATIO);
}

void trackControl(GlobalVariables& gv)
{
  floatControl(gv);  // Remove this line when you implement this controller
}

void pfmoveControl(GlobalVariables& gv)
{
  floatControl(gv);  // Remove this line when you implement this controller
}

// bonus (a): line back and forth in operational space
void lineControl(GlobalVariables& gv)
{
  setLineTarget(gv, gv.curTime - bonus_start_time);
  operationalSpaceControl(gv);
}

void proj1Control(GlobalVariables& gv)
{
  followSpline(gv);
}

void proj2Control(GlobalVariables& gv)
{
  // constant angular velocity: beta grows linearly with time
  double beta = CIRCLE_ANGULAR_VELOCITY * (gv.curTime - circle_start_time);
  setCircleTarget(gv, beta, CIRCLE_ANGULAR_VELOCITY);
  operationalSpaceControl(gv);
}

// Linear segment with parabolic blends for the circle angle beta:
// - blend 1 (0 <= t < tb): constant acceleration, dbeta grows linearly,
//   beta = 1/2 * acc * t^2
// - linear segment: constant angular velocity 2 pi / 5 s
// - blend 2 (last tb seconds): constant deceleration down to dbeta = 0
// - afterwards the end effector holds the end point (beta = 6 pi)
// The blend time follows from dbeta(tb) = tb * acc (hint in the sheet):
//   tb = (2 pi / 5) / (2 pi / 25) = 5 s
// Each blend covers 1/2 * acc * tb^2 = pi, so the linear segment covers the
// remaining 6 pi - 2 pi = 4 pi. The total time is
//   tf = 2 tb + 4 pi / velocity = total angle / velocity + tb = 20 s
void proj3Control(GlobalVariables& gv)
{
  double t = gv.curTime - circle_start_time;
  double acc = CIRCLE_ANGULAR_ACCELERATION;
  double vel = CIRCLE_ANGULAR_VELOCITY;
  double tb = vel / acc;
  double tf = CIRCLE_TOTAL_ANGLE / vel + tb;

  double beta, dbeta;
  if (t < tb) {
    beta = 0.5 * acc * t * t;
    dbeta = acc * t;
  } else if (t < tf - tb) {
    beta = 0.5 * acc * tb * tb + vel * (t - tb);
    dbeta = vel;
  } else if (t < tf) {
    // mirror image of blend 1, measured backwards from the end
    double t_left = tf - t;
    beta = CIRCLE_TOTAL_ANGLE - 0.5 * acc * t_left * t_left;
    dbeta = acc * t_left;
  } else {
    beta = CIRCLE_TOTAL_ANGLE;
    dbeta = 0.0;
  }

  setCircleTarget(gv, beta, dbeta);
  operationalSpaceControl(gv);
}

// *******************************************************************
// Debug function
// *******************************************************************

void PrintDebug(GlobalVariables& gv)
{
  // Replace this code with any debug information you'd like to get
  // when you type "pdebug" at the prompt.
  printf("This sample code prints the torque and mass\n");
  gv.tau.display("tau");
  gv.A.display("A");
}

#ifdef WIN32
// *******************************************************************
// XPrintf(): Replacement for printf() which calls ui->VDisplay()
// whenever the ui object is available.  See utility/XPrintf.h.
// *******************************************************************

int XPrintf(const char* fmt, ...)
{
  int returnValue;
  va_list argptr;
  va_start(argptr, fmt);

  returnValue = vprintf(fmt, argptr);

  va_end(argptr);
  return returnValue;
}
#endif //#ifdef WIN32

/********************************************************

END OF DEFAULT STUDENT FILE

ADD HERE ALL STUDENT DEFINED AND AUX FUNCTIONS

*******************************************************/
