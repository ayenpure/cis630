#include <iostream>

using std::cout;
using std::cout;

int main() {
  int Y[3] = {5,7,2};
  cout << ((Y[0] < Y[1]) ? Y[0] : ((Y[1] < Y[2]) ? Y[1] : Y[2])) << endl;
  cout << ((Y[0] > Y[1]) ? Y[0] : ((Y[1] > Y[2]) ? Y[1] : Y[2])) << endl;
}
