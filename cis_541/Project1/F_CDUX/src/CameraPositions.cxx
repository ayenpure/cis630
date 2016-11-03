#include <iostream>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include "CameraPositions.h"

using std::sin;
using std::cos;
using std::endl;
using std::cout;
using std::abs;

#define X 'x'
#define Y 'y'
#define Z 'z'

void rotate(double* camera_position, double angle, char axis,double* rotated_camera) {
  double rotation_matrix[3][3];
  if(axis == X) {
    double x_rotation_matrix[3][3] = {
      1,0,0,
      0,cos(angle),-sin(angle),
      0,sin(angle),cos(angle)
    };
    memcpy(rotation_matrix, x_rotation_matrix, 9*sizeof(double));
  } else if (axis == Y) {
    double y_rotation_matrix[3][3] = {
      cos(angle),0,-sin(angle),
      0,1,0,
      sin(angle),0,cos(angle)
    };
    memcpy(rotation_matrix, y_rotation_matrix, 9*sizeof(double));
  } else if (axis == Z) {
    double z_rotation_matrix[3][3] = {
      cos(angle),-sin(angle),0,
      sin(angle),cos(angle),0,
      0,0,1
    };
    memcpy(rotation_matrix, z_rotation_matrix, 9*sizeof(double));
  }
  /*for (int i = 0 ; i < 4 ; i++)
  {
      char str[256];
      sprintf(str, "(%.7f %.7f %.7f)\n", rotation_matrix[i][0], rotation_matrix[i][1], rotation_matrix[i][2]);
      cout << str;
  }*/
  for(int i = 0; i < 3; i++) {
    rotated_camera[i] = 0;
    for(int j = 0; j < 3; j++) {
      rotated_camera[i] = rotated_camera[i] + camera_position[j]*rotation_matrix[j][i];
    }
    if(abs(rotated_camera[i]) < 0.0000001)
      rotated_camera[i] = 0;
  }
  //cout << "rotated position {" << rotated_camera[0] << ", " << rotated_camera[1] << ", " << rotated_camera[2] << " }" << endl;
}

void get_camera_positions(double (*camera_positions)[3]) {
  double calc_camera_positions[8][16][3];
  cout << "PI :" << M_PI << endl;
  double first_rotation[16] = {
    0, (M_PI/6), (M_PI/4), (M_PI/3),
    (M_PI/2), (2*M_PI/3), (3*M_PI/4), (5*M_PI/6),
    (M_PI), (7*M_PI/6), (5*M_PI/4), (4*M_PI/3),
    (3*M_PI/2), (5*M_PI/3),(7*M_PI/4) ,(11*M_PI/6)
  };
  double second_rotation[7] = {
    (M_PI/6), (M_PI/4), (M_PI/3), (M_PI/2), (2*M_PI/3), (3*M_PI/4), (5*M_PI/6)
  };
  double camera_position[3] = {40,0,0};
  for(int i = 0; i < 16; i++) {
    rotate(camera_position, first_rotation[i], Y, calc_camera_positions[0][i]);
  };
  for(int i = 1; i < 8; i++) {
    for(int j = 0; j < 16; j++) {
      rotate(calc_camera_positions[0][j], second_rotation[i-1], X, calc_camera_positions[i][j]);
    }
  };
  int positions_count = 0;
  for(int i = 0; i < 8 ; i++) {
    for (int j = 0; j <16; j++) {
      if(i > 0 && (j==0 || j==8))
        continue;
      camera_positions[positions_count][0] = calc_camera_positions[i][j][0];
      camera_positions[positions_count][1] = calc_camera_positions[i][j][1];
      camera_positions[positions_count][2] = calc_camera_positions[i][j][2];
      positions_count++;
    }
  }
}
