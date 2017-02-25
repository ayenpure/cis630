#include <iostream>
#include <string.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <sstream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetWriter.h>
#include "CameraPositions.h"
#include "Utilities.h"
#include "TriangleOperations.h"
#include "MatrixOperations.h"
#include "Camera.h"
#include "Screen.h"
#include "RenderFunctions.h"

using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;
using std::pow;
using std::tan;
using std::sin;

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

int main(int argc, char *argv[]) {
cerr << "In the right main (camera active pixels)!!" << endl;
	int config_id = 0,
	 		no_of_procs = 0,
			num_cameras = 0,
			read_files = 0;
	std::vector<Triangle> triangles;
	if(argc < 4) {
		cout << "Incorrect number of arguments for execution" << endl;
		exit (EXIT_FAILURE);
	} else {
cerr << "Reading ... ";
		no_of_procs = atoi(argv[1]);
		config_id = atoi(argv[2]);
		read_files = atoi(argv[3]);
		if(!read_files)
			triangles = GetTriangles(2,0,7);
		else {
			triangles = GetTrianglesFromFiles(no_of_procs, argv[4]);
		}
cerr << "done!" << endl;
	}
	double** camera_positions = get_camera_positions(config_id,&num_cameras);
	int active_pixels[num_cameras];
	int no_of_triangles = triangles.size();
	for (int i = 0; i < 1/*num_cameras*/; i++)
		active_pixels[i] = 0;
	for (int cam_index = 0; cam_index < 1; cam_index++) {
cerr << "Working on camera " << cam_index << endl;
		vtkImageData *image = NewImage(1000, 1000);
		unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0, 0,
				0);
		int npixels = 1000 * 1000;
		for (int i = 0; i < npixels * 3; i++)
			buffer[i] = 0;
		double *depth_buffer = (double*)malloc(npixels*sizeof(double));
		for (int i = 0; i < npixels; i++)
			depth_buffer[i] = -1;
		Screen screen;
		screen.buffer = buffer;
		screen.depth_buffer = depth_buffer;
		screen.width = 1000;
		screen.height = 1000;
		double camera_position[3], focus_point[3];
		get_camera_position_and_focus(config_id, cam_index, num_cameras,
			 camera_positions, camera_position, focus_point);
		Camera camera = GetCamera(camera_position, focus_point);
		Matrix camera_transform = camera.CameraTransform();
		cout<<"\nCamera Transform Matrix :\n";camera_transform.Print(std::cout);
		Matrix view_transform = camera.ViewTransform();
		cout<<"\nView Transform Matrix :\n";view_transform.Print(std::cout);
		Matrix device_transform = camera.DeviceTransform(screen);
		cout<<"\nDevice Transform Matrix :\n";device_transform.Print(std::cout);
		Matrix composite = get_total_transform_matrix(camera_transform,
				view_transform, device_transform);
		cout<<"\nComposite Matrix :\n";composite.Print(std::cout);
//cerr << "Working on " << triangles.size() << " triangles" << endl;
		for (int vecIndex = 0; vecIndex < 1/*triangles.size()*/; vecIndex++) {
			//cout << "Triangle #" << vecIndex << endl;
			Triangle t;
			t.X[0] = 4.6122;
			t.Y[0] = -5.00049;
			t.Z[0] = 9.5;
			t.X[1] = 4.6122;
			t.Y[1] = -5.02041;
			t.Z[1] = 9.418;
			t.X[2] = 4.5714;
			t.Y[2] = -4.99509;
			t.Z[2] = 9.5;
			//Triangle t = triangles[vecIndex];
			cerr << "Before transform, triangle is " << endl;
			t.Print(cerr);
			//print_triangle(t);
			transformTriangle(&t, composite, camera);
			//print_triangle(t);
			cerr << "After transform, triangle is " << endl;
			t.Print(cerr);

			if (t.is_flat_bottom_triangle()) {
				/*triangle_pixels[vecIndex] = triangle_pixels[vecIndex] + */scan_line(&t, &screen);
			} else {
				Triangle t1, t2;
				t.split_triangle(&t1, &t2);
				/*triangle_pixels[vecIndex] = triangle_pixels[vecIndex] + */scan_line(&t1, &screen);
				/*triangle_pixels[vecIndex] = triangle_pixels[vecIndex] + */scan_line(&t2, &screen);
			}
		}
		std::ostringstream oss;
		oss << "camera_" << cam_index;
		active_pixels[cam_index] = 0;
cerr << "Counting non-black pixels" << endl;
		for(int pxl =0; pxl < 3*npixels; pxl+=3) {
			if(!(buffer[pxl] == 0 &&
				buffer[pxl+1] == 0 &&
				buffer[pxl+2] == 0 ))
				active_pixels[cam_index]++;
		}
		active_pixels[cam_index] = active_pixels[cam_index] / no_of_procs;
		WriteImage(image, oss.str().c_str());
		oss.str("");
		oss.clear();
	}
	for (int i = 0; i < num_cameras; i++) {
		cout << active_pixels[i] << endl;
	}
}
