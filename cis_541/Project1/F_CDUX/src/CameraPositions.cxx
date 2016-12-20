#include <iostream>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
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

void get_config_random(double *camera_positions[3]) {
  double calc_camera_positions[8][16][3];
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

void get_config_helix(double *camera_positions[3]) {
  double y = -15;
  int cam_index = 0;
  double camera_position[3] = {40, y,0};
  double rotated_camera[3];
  camera_positions[cam_index][0] = camera_position[0];
  camera_positions[cam_index][1] = camera_position[1];
  camera_positions[cam_index][2] = camera_position[2];
  //cout << camera_position[0] << ", " << camera_position[1] << ", " << camera_position[2] << endl;
  cam_index++;
  do {
    rotate(camera_position, M_PI / 8, Y, rotated_camera);
    camera_position[0] = rotated_camera[0];
    y += .375;
    camera_position[1] = y;
    camera_position[2] = rotated_camera[2];
    //cout << camera_position[0] << ", " << camera_position[1] << ", " << camera_position[2] << endl;
    camera_positions[cam_index][0] = camera_position[0];
    camera_positions[cam_index][1] = camera_position[1];
    camera_positions[cam_index][2] = camera_position[2];
    cam_index++;
  } while ( y < 15);
}

void get_camera_position_and_focus(int config_id, int cam_index, int num_cameras ,double** camera_positions,double* camera_position,double* focus_point) {
  int focus_index = 0;
  camera_position[0] = camera_positions[cam_index][0];
  camera_position[1] = camera_positions[cam_index][1];
  camera_position[2] = camera_positions[cam_index][2];
  switch(config_id) {
    default :
    case 0 :
      focus_point[0] = 0;
      focus_point[1] = 0;
      focus_point[2] = 0;
    break;
    case 1 :
      if(cam_index < 16)
        focus_index = (cam_index + 6) % 16;
      else {
        int quotient = (cam_index - 16) / 14;
        int pseudo_index = (cam_index - 16) % 14;
        int pseudo_focus_index = (pseudo_index + 5) % 14;
        focus_index = ((quotient*14) + 16) + pseudo_focus_index;
      }
      focus_point[0] = camera_positions[focus_index][0];
      focus_point[1] = camera_positions[focus_index][1];
      focus_point[2] = camera_positions[focus_index][2];
    break;
    case 2 :
      if(cam_index < 16)
        focus_index = (cam_index + 7) % 16;
      else {
        int quotient = (cam_index - 16) / 14;
        int pseudo_index = (cam_index - 16) % 14;
        int pseudo_focus_index = (pseudo_index + 6) % 14;
        focus_index = ((quotient*14) + 16) + pseudo_focus_index;
      }
      focus_point[0] = camera_positions[focus_index][0];
      focus_point[1] = camera_positions[focus_index][1];
      focus_point[2] = camera_positions[focus_index][2];
    break;
    case 3 :
      focus_index = cam_index + 7;
      focus_point[0] = camera_positions[focus_index][0];
      focus_point[1] = camera_positions[focus_index][1];
      focus_point[2] = camera_positions[focus_index][2];
    break;
    case 4 :
      focus_index = cam_index + 7;
      focus_point[0] = 0;
      focus_point[1] = camera_positions[cam_index][1];
      focus_point[2] = 0;
    break;
  }
}

double** allocate_memory(int num_cameras) {
  double **camera_positions;
  camera_positions = (double**)malloc(num_cameras*sizeof(double*));
  for(int i = 0; i < num_cameras; i++)
    camera_positions[i] = (double*)malloc(3*sizeof(double));
  return camera_positions;
}

double** get_camera_positions(int config_id, int *num_cameras) {
  double **camera_positions;
  switch(config_id) {
    case 0:
    case 1:
    case 2:
    default :
      *num_cameras = 114;
      camera_positions = allocate_memory(114);
      get_config_random(camera_positions);
    break;
    case 3:
    case 4:
      *num_cameras = 74;
      camera_positions = allocate_memory(81);
      get_config_helix(camera_positions);
    break;
  }
  return camera_positions;
}
