#include <iostream>
#include <cmath>

using std::sin;
using std::cos;
using std::endl;
using std::cout;

void multiply(double* camera_position, double rotation_matrix[][3]) {
  double rotated_camera[3];
  for(int i = 0; i < 3; i++) {
    rotated_camera[i] = 0;
    for(int j = 0; j < 3; j++) {
      rotated_camera[i] = rotated_camera[i] + camera_position[j]*rotation_matrix[j][i];
    }
  }
  cout << "rotated position {" << rotated_camera[0] << ", " << rotated_camera[1] << ", " << rotated_camera[2] << " }" << endl;
}

int main() {
  double camera_position[3] = {20,0,0};
  double rotation_matrix[3][3] = {
    1,0,0,
    0,1,0,
    0,0,1
  };
  multiply(camera_position,rotation_matrix);
  return 0;
}
