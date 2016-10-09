#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

using std::cerr;
using std::endl;

vtkImageData *
NewImage(int width, int height) {
	vtkImageData *img = vtkImageData::New();
	img->SetDimensions(width, height, 1);
	img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

	return img;
}

void WriteImage(vtkImageData *img, const char *filename) {
	std::string full_filename = filename;
	full_filename += ".png";
	vtkPNGWriter *writer = vtkPNGWriter::New();
	writer->SetInputData(img);
	writer->SetFileName(full_filename.c_str());
	writer->Write();
	writer->Delete();
}

int main() {
	std::cerr << "In main!" << endl;
	vtkImageData *image = NewImage(1024, 1350);
	unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0, 0, 0);
	int buffer_index = 0;
	for (int strip_num = 0; strip_num < 27; strip_num++) {
		//Calculate RGB values for the current strip
		int red, blue, green;
		red = (strip_num / 9 == 0) ? 0 : ((strip_num / 9 == 1) ? 128 : 255);
		green = ((strip_num / 3) % 3 == 0) ? 0 : (((strip_num / 3) % 3 == 1) ? 128 : 255);
		blue = (strip_num % 3 == 0) ? 0 : ((strip_num % 3 == 1) ? 128 : 255);
		//Paint the strip
		for (int x_co = 0; x_co < 50; x_co++) {
			for (int y_co = 0; y_co < 1024; y_co++) {
				buffer[buffer_index++] = red;
				buffer[buffer_index++] = green;
				buffer[buffer_index++] = blue;
			}
		}
	}
	WriteImage(image, "project1A");
}
