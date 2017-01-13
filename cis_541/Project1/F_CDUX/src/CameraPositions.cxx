#include <iostream>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "CameraPositions.h"
#include "Utilities.h"
#include "TriangleOperations.h"

using std::sin;
using std::cos;
using std::endl;
using std::cout;
using std::abs;

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
    rotate(camera_position, first_rotation[i], 'y', calc_camera_positions[0][i]);
  };
  for(int i = 1; i < 8; i++) {
    for(int j = 0; j < 16; j++) {
      rotate(calc_camera_positions[0][j], second_rotation[i-1], 'x', calc_camera_positions[i][j]);
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

void get_uniform_cinema_distribution(double *camera_positions[3]) {
  std::vector< std::vector<Triangle> > proc_triangles;
  proc_triangles = GetTrianglesForProcs(2,1,0);
  cout << proc_triangles.size() << endl;
  for(int i=0;i<proc_triangles.size();i++) {
    Triangle t = proc_triangles[i][0];
    for(int k = 0; k < 3; k++) {
      double ptMag = sqrt(t.X[k]*t.X[k]+
                          t.Y[k]*t.Y[k]+
                          t.Z[k]*t.Z[k]);
      t.X[k] = (t.X[k] / ptMag)*40;
      t.Y[k] = (t.Y[k] / ptMag)*40;
      t.Z[k] = (t.Z[k] / ptMag)*40;
    }
    camera_positions[i][0] = (t.X[0] + t.X[1] + t.X[2]) / 3;
    camera_positions[i][1] = (t.Y[0] + t.Y[1] + t.Y[2]) / 3;
    camera_positions[i][2] = (t.Z[0] + t.Z[1] + t.Z[2]) / 3;
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
    rotate(camera_position, M_PI / 8, 'y', rotated_camera);
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
    case 5:
      *num_cameras = 64;
      camera_positions = allocate_memory(64);
      get_uniform_cinema_distribution(camera_positions);
    break;
  }
  return camera_positions;
}
