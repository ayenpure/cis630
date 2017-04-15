#include <iostream>
#include <vtkm/Types.h>

using std::cout;
using std::endl;

int main() {
        vtkm::Vec<vtkm::Int32,3> vec(10,11,12);
        for(vtkm::IdComponent index = 0; index < vec.NUM_COMPONENTS; index++ ) {
                cout << "value : " << vec[index] << endl;
        }
        return 0;
}
