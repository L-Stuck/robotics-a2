/*

 Please fill in the function stubs provided in this file

 For achieving full points on this exercise, explain the implemented equations
 in the accompanying pdf. Especially:
   - which parameters were used for computation?
   - were equations simplified? What was the original equation?
   - when introducing additional variables, what do they contain?


 * Do not use additional includes
 * Do not include/link additional files
 * Do not change the predefined function signature
 * Do not change the header file

 Tip: use the main() function to test your implementation.

*/

#include "forwardkinematics.hpp" //Have a look at this header file, it declares the class ForwardKinematicsPuma2D
#include <cmath>                 //use sin, cos from cmath for your computation
#include <iostream>
using namespace std;


#ifdef UNITTEST
    #define main STUDENTS_MAIN
#endif


/*
Convenience function to print out a homogenous transform
*/
void print_HTransform(HTransform tf)
{
    std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
    std::cout.precision(3);
    cout << "\n---------------------------\n";

    cout << tf[0][0]<<"  "<<tf[0][1]<<"  "<<tf[0][2]<<"  "<<tf[0][3]<<endl;
    cout << tf[1][0]<<"  "<<tf[1][1]<<"  "<<tf[1][2]<<"  "<<tf[1][3]<<endl;
    cout << tf[2][0]<<"  "<<tf[2][1]<<"  "<<tf[2][2]<<"  "<<tf[2][3]<<endl;
    cout << tf[3][0]<<"  "<<tf[3][1]<<"  "<<tf[3][2]<<"  "<<tf[3][3]<<endl;
    cout <<   "---------------------------\n";
}

/*
Convenience function to print out a 3x3 Jacobian matrix
*/
void print_Jacobian(float Jacobian[3][3])
{
    std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
    std::cout.precision(3);
    cout << "\n+++++++++++++++++++++++++++\n";
    cout << Jacobian[0][0]<<"  "<<Jacobian[0][1]<<"  "<<Jacobian[0][2]<<endl;
    cout << Jacobian[1][0]<<"  "<<Jacobian[1][1]<<"  "<<Jacobian[1][2]<<endl;
    cout << Jacobian[2][0]<<"  "<<Jacobian[2][1]<<"  "<<Jacobian[2][2]<<endl;
    cout <<   "+++++++++++++++++++++++++++\n";
}


/*
Convenience function to print out a position
*/
void print_Position(float F[3])
{
    std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
    std::cout.precision(3);
    cout << "\nooooooooooooooooooooooooooo\n";
    cout << F[0]<<"  "<<F[1]<<"  "<< F[2]<<endl;
    cout <<   "ooooooooooooooooooooooooooo\n";
}


/*
already implemented

Set the robot's joint values and recompute the forward kinematics.

a1, a2 and a3 are assumed to be given in radians!
*/
void ForwardKinematicsPuma2D::setJoints(float a1, float a2, float a3)
{
    //store joint angles
    angles[0] = a1;
    angles[1] = a2;
    angles[2] = a3;
    
    //recompute the dependent variables
    computeDH();
    computeT0_1();
    computeT1_2();
    computeT2_3();
    computeT3_E();
    computeT0_E();
    computeF();
    computeJ();
}




/***********************************************/
/******************EDIT BELOW ******************/
/***********************************************/

/*
pi/2 in radians. Needed because link 1 points to the right in the zero
configuration while links 2 and 3 point upwards (see row 2 of the DH table).
*/
const float HALF_PI = 1.57079633f;

/*
updates the variable T0_1

General DH transform between frame i-1 and frame i (Craig, p. 74):

   [ c(theta_i)              -s(theta_i)               0              a_(i-1)             ]
   [ s(theta_i)c(alpha_(i-1)) c(theta_i)c(alpha_(i-1)) -s(alpha_(i-1)) -s(alpha_(i-1))d_i ]
   [ s(theta_i)s(alpha_(i-1)) c(theta_i)s(alpha_(i-1))  c(alpha_(i-1))  c(alpha_(i-1))d_i ]
   [ 0                        0                         0              1                   ]

All joint axes of the planar arm are parallel, so alpha = 0 and d = 0 in every
row. This simplifies every transform to a rotation about z by theta_i plus a
translation a_(i-1) along x:

   [ c(theta_i) -s(theta_i)  0  a_(i-1) ]
   [ s(theta_i)  c(theta_i)  0  0       ]
   [ 0           0           1  0       ]
   [ 0           0           0  1       ]

Row 1 of the DH table: a0 = 0, theta1 = q1. Frame 0 and frame 1 share their
origin, so T0_1 is a pure rotation by q1.
*/
void ForwardKinematicsPuma2D::computeT0_1()
{
   float c1 = cos(angles[0]);
   float s1 = sin(angles[0]);

   //row vector
   T0_1[0][0] = c1;
   T0_1[0][1] = -s1;
   T0_1[0][2] = 0.0;
   T0_1[0][3] = 0.0;

   //row vector
   T0_1[1][0] = s1;
   T0_1[1][1] = c1;
   T0_1[1][2] = 0.0;
   T0_1[1][3] = 0.0;

   //row vector
   T0_1[2][0] = 0.0;
   T0_1[2][1] = 0.0;
   T0_1[2][2] = 1.0;
   T0_1[2][3] = 0.0;

   //row vector
   T0_1[3][0] = 0.0;
   T0_1[3][1] = 0.0;
   T0_1[3][2] = 0.0;
   T0_1[3][3] = 1.0;
}


/*
updates the variable T1_2

Row 2 of the DH table: a1 = l1, theta2 = q2 - pi/2.
Frame 2 sits at the end of link 1 (distance l1 along x1). In the zero
configuration x2 points up while x1 points to the right, which is a rotation of
-pi/2 about z (z points into the drawing plane, positive angles are clockwise).
With the identities
   cos(q2 - pi/2) = sin(q2) = s2
   sin(q2 - pi/2) = -cos(q2) = -c2
the rotation part becomes [ s2 c2 ; -c2 s2 ].
*/
void ForwardKinematicsPuma2D::computeT1_2()
{
   float c2 = cos(angles[1]);
   float s2 = sin(angles[1]);

   //row vector
   T1_2[0][0] = s2;
   T1_2[0][1] = c2;
   T1_2[0][2] = 0.0;
   T1_2[0][3] = l1;

   //row vector
   T1_2[1][0] = -c2;
   T1_2[1][1] = s2;
   T1_2[1][2] = 0.0;
   T1_2[1][3] = 0.0;

   //row vector
   T1_2[2][0] = 0.0;
   T1_2[2][1] = 0.0;
   T1_2[2][2] = 1.0;
   T1_2[2][3] = 0.0;

   //row vector
   T1_2[3][0] = 0.0;
   T1_2[3][1] = 0.0;
   T1_2[3][2] = 0.0;
   T1_2[3][3] = 1.0;
}

/*
updates the variable T2_3

Row 3 of the DH table: a2 = l2, theta3 = q3.
Frame 3 sits at the end of link 2 (distance l2 along x2). x2 and x3 both point
up in the zero configuration, so there is no constant offset in theta3.
*/
void ForwardKinematicsPuma2D::computeT2_3()
{
   float c3 = cos(angles[2]);
   float s3 = sin(angles[2]);

   //row vector
   T2_3[0][0] = c3;
   T2_3[0][1] = -s3;
   T2_3[0][2] = 0.0;
   T2_3[0][3] = l2;

   //row vector
   T2_3[1][0] = s3;
   T2_3[1][1] = c3;
   T2_3[1][2] = 0.0;
   T2_3[1][3] = 0.0;

   //row vector
   T2_3[2][0] = 0.0;
   T2_3[2][1] = 0.0;
   T2_3[2][2] = 1.0;
   T2_3[2][3] = 0.0;

   //row vector
   T2_3[3][0] = 0.0;
   T2_3[3][1] = 0.0;
   T2_3[3][2] = 0.0;
   T2_3[3][3] = 1.0;
}


/*
updates the variable T3_E

Row 4 (E) of the DH table: a3 = l3, theta4 = 0.
The end effector frame is shifted by l3 along x3 and has the same orientation as
frame 3 (there is no joint between frame 3 and E). The transform is constant.
*/
void ForwardKinematicsPuma2D::computeT3_E()
{
   //row vector
   T3_E[0][0] = 1.0;
   T3_E[0][1] = 0.0;
   T3_E[0][2] = 0.0;
   T3_E[0][3] = l3;

   //row vector
   T3_E[1][0] = 0.0;
   T3_E[1][1] = 1.0;
   T3_E[1][2] = 0.0;
   T3_E[1][3] = 0.0;

   //row vector
   T3_E[2][0] = 0.0;
   T3_E[2][1] = 0.0;
   T3_E[2][2] = 1.0;
   T3_E[2][3] = 0.0;

   //row vector
   T3_E[3][0] = 0.0;
   T3_E[3][1] = 0.0;
   T3_E[3][2] = 0.0;
   T3_E[3][3] = 1.0;
}


/*

This function updates the variable T0_E

T0_E = T0_1 * T1_2 * T2_3 * T3_E, simplified analytically:
- All four transforms rotate about the same z axis, so their rotation angles
  add up: q1 + (q2 - pi/2) + q3 + 0. With cos(q123 - pi/2) = s123 and
  sin(q123 - pi/2) = -c123 the rotation part is [ s123 c123 ; -c123 s123 ].
- The position is the sum of the three link vectors, each one rotated by the
  sum of all angles before it:
     l1 * (c1, s1)
   + l2 * (cos(q12 - pi/2), sin(q12 - pi/2))   = l2 * (s12, -c12)
   + l3 * (cos(q123 - pi/2), sin(q123 - pi/2)) = l3 * (s123, -c123)

   [ s123  c123  0  l1*c1 + l2*s12 + l3*s123 ]
   [-c123  s123  0  l1*s1 - l2*c12 - l3*c123 ]
   [ 0     0     1  0                        ]
   [ 0     0     0  1                        ]

Note that y0 points downwards in Fig. 1, so the end effector above the base has
a negative y coordinate.
*/
void ForwardKinematicsPuma2D::computeT0_E()
{
   float c1 = cos(angles[0]);
   float s1 = sin(angles[0]);
   float c12 = cos(angles[0] + angles[1]);
   float s12 = sin(angles[0] + angles[1]);
   float c123 = cos(angles[0] + angles[1] + angles[2]);
   float s123 = sin(angles[0] + angles[1] + angles[2]);

   //row vector
   T0_E[0][0] = s123;
   T0_E[0][1] = c123;
   T0_E[0][2] = 0.0;
   T0_E[0][3] = l1 * c1 + l2 * s12 + l3 * s123;

   //row vector
   T0_E[1][0] = -c123;
   T0_E[1][1] = s123;
   T0_E[1][2] = 0.0;
   T0_E[1][3] = l1 * s1 - l2 * c12 - l3 * c123;

   //row vector
   T0_E[2][0] = 0.0;
   T0_E[2][1] = 0.0;
   T0_E[2][2] = 1.0;
   T0_E[2][3] = 0.0;

   //row vector
   T0_E[3][0] = 0.0;
   T0_E[3][1] = 0.0;
   T0_E[3][2] = 0.0;
   T0_E[3][3] = 1.0;
}



/*
This function updates the variables ee_x, ee_y, ee_alpha

x and y are the translation part of T0_E (last column).
alpha is the angle of the end effector x axis (xE) measured from x0 about z0.
It is the sum of all DH joint angles theta1 + theta2 + theta3 + theta4,
which is the same as atan2(T0_E[1][0], T0_E[0][0]):
   alpha = q1 + q2 + q3 - pi/2
In the zero configuration the end effector points up, i.e. alpha = -pi/2,
because y0 points down.
*/
void ForwardKinematicsPuma2D::computeF()
{
   F[0] = T0_E[0][3]; //x
   F[1] = T0_E[1][3]; //y
   F[2] = angles[0] + angles[1] + angles[2] - HALF_PI; //alpha
}


/*
This function updates the variable J

J = dF/dq, i.e. J[i][j] is the partial derivative of F[i] with respect to q_j.
Derivatives of the components of F:
   x     = l1*c1 + l2*s12 + l3*s123
   y     = l1*s1 - l2*c12 - l3*c123
   alpha = q1 + q2 + q3 - pi/2
q1 appears in every angle sum, q2 only in s12/c12 and s123/c123, q3 only in
s123/c123. Therefore the columns get shorter from left to right:
   dx/dq1 = -l1*s1 + l2*c12 + l3*c123,  dx/dq2 = l2*c12 + l3*c123,  dx/dq3 = l3*c123
   dy/dq1 =  l1*c1 + l2*s12 + l3*s123,  dy/dq2 = l2*s12 + l3*s123,  dy/dq3 = l3*s123
   dalpha/dq1 = dalpha/dq2 = dalpha/dq3 = 1
*/
void ForwardKinematicsPuma2D::computeJ()
{
   float c1 = cos(angles[0]);
   float s1 = sin(angles[0]);
   float c12 = cos(angles[0] + angles[1]);
   float s12 = sin(angles[0] + angles[1]);
   float c123 = cos(angles[0] + angles[1] + angles[2]);
   float s123 = sin(angles[0] + angles[1] + angles[2]);

   //row vector
   J[0][0] = -l1 * s1 + l2 * c12 + l3 * c123;
   J[0][1] = l2 * c12 + l3 * c123;
   J[0][2] = l3 * c123;

   //row vector
   J[1][0] = l1 * c1 + l2 * s12 + l3 * s123;
   J[1][1] = l2 * s12 + l3 * s123;
   J[1][2] = l3 * s123;

   //row vector
   J[2][0] = 1.0;
   J[2][1] = 1.0;
   J[2][2] = 1.0;
}

/* This function computes the DH params
   and sets the internal member variable DH

Each row i holds (alpha_(i-1), a_(i-1), d_i, theta_i) for the frames 1, 2, 3
and E, read off Fig. 1 (Craig convention). Angles are in radians.
- All z axes point into the drawing plane and are parallel: alpha = 0, d = 0.
- a_(i-1) is the distance between the joint axes along x_(i-1):
  0 (frames 0 and 1 coincide), l1, l2, and l3 to the end effector.
- theta_i is the angle from x_(i-1) to x_i about z_i: the joint angle, plus
  -pi/2 in row 2 because x1 points right and x2 points up in the zero
  configuration. There is no joint between frame 3 and E, so theta4 = 0.
*/
void ForwardKinematicsPuma2D::computeDH()
{
   //row vector: frame 1
   DH[0][0] = 0.0;
   DH[0][1] = 0.0;
   DH[0][2] = 0.0;
   DH[0][3] = angles[0];

   //row vector: frame 2
   DH[1][0] = 0.0;
   DH[1][1] = l1;
   DH[1][2] = 0.0;
   DH[1][3] = angles[1] - HALF_PI;

   //row vector: frame 3
   DH[2][0] = 0.0;
   DH[2][1] = l2;
   DH[2][2] = 0.0;
   DH[2][3] = angles[2];

   //row vector: frame E
   DH[3][0] = 0.0;
   DH[3][1] = l3;
   DH[3][2] = 0.0;
   DH[3][3] = 0.0;
}


/*
Example code to test your functions:

You are free to change main() as you like
*/
int main()
{
 ForwardKinematicsPuma2D* fk = new ForwardKinematicsPuma2D();
 fk->setJoints(0.0,0.0,0.0); //example, try out different values
 cout << "********************Testing Transforms**************"<<endl;
 print_HTransform(fk->T0_1);
 print_HTransform(fk->T1_2);
 print_HTransform(fk->T2_3);
 print_HTransform(fk->T3_E);
 print_HTransform(fk->T0_E);
 cout << "********************Testing F***********************"<<endl;
 print_Position(fk->F);
 cout << "********************Testing J***********************"<<endl;
 print_Jacobian(fk->J);
 return 0;
}
