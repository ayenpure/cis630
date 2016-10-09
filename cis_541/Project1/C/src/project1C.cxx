#include <iostream>
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

using std::cerr;
using std::endl;

double ceil441(double f) {
	return ceil(f - 0.00001);
}

double floor441(double f) {
	return floor(f + 0.00001);
}

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

class Triangle {
public:
	double X[3];
	double Y[3];
	unsigned char color[3];
	/*
	 * The offset_vertex is the vertex[X,Y] that does not belong to the base
	 * The left_vertex is the vertex[X,Y] that is to the left side of the base
	 * The right_vertex is the vertex[X,Y] that is to the rigth side of the base
	 * slope_on right will hold the slope of the triangle's side on the right and
	 * similarly slope_on_left will hold the slope of the side on the left of the
	 * triangle
	 * left_offset and right_offset are the constants in the line equation for the
	 * sides of the triangle
	 */
	double offset_vertex[2];
	double left_vertex[2], slope_on_left, left_offset;
	double right_vertex[2], slope_on_right, right_offset;

	double getlowestY() {
		return ((Y[0] < Y[1]) ? Y[0] : ((Y[1] < Y[2]) ? Y[1] : Y[2]));
	}

	double gethighestY() {
		return ((Y[0] > Y[1]) ? Y[0] : ((Y[1] > Y[2]) ? Y[1] : Y[2]));
	}

	double getlowestX() {
		return ((X[0] < X[1]) ? X[0] : ((X[1] < X[2]) ? X[1] : X[2]));
	}

	double gethighestX() {
		return ((X[0] > X[1]) ? X[0] : ((X[1] > X[2]) ? X[1] : X[2]));
	}

	void determine_triangle_orientation() {
		/*
		 * We determine the vertices that represent the base and the other
		 * offset vertex here.
		 */
		int top_index = ((Y[0] == Y[1]) ? 2 : ((Y[1] == Y[2]) ? 0 : 1));
		offset_vertex[0] = X[top_index];
		offset_vertex[1] = Y[top_index];
		int left_index = 0, right_index = 0;
		if (X[(top_index + 1) % 3] < X[(top_index + 2) % 3]) {
			left_index = (top_index + 1) % 3;
			right_index = (top_index + 2) % 3;
		} else {
			left_index = (top_index + 2) % 3;
			right_index = (top_index + 1) % 3;
		}
		left_vertex[0] = X[left_index];
		left_vertex[1] = Y[left_index];
		right_vertex[0] = X[right_index];
		right_vertex[1] = Y[right_index];

		/*
		 * We solve for the slope and the offset for the sides of the triangle.
		 * the slope is calculated by two point form slope = (Y2 - Y1) / (X2 - X1)
		 * And the offsets are calculated as c = y - slope * x;
		 */
		double diff;
		diff = (offset_vertex[0] - left_vertex[0]);
		slope_on_left =
				(diff == 0) ? 0 : (offset_vertex[1] - left_vertex[1]) / diff;
		left_offset = left_vertex[1] - (double(slope_on_left * left_vertex[0]));

		diff = (offset_vertex[0] - right_vertex[0]);
		slope_on_right =
				(diff == 0) ? 0 : (offset_vertex[1] - right_vertex[1]) / diff;
		right_offset = right_vertex[1]
				- (double(slope_on_right * right_vertex[0]));

	}

	double get_left_x_intercept(int y) {
		/*
		 * We solve for the intercept for the scan line on the left side of the
		 * triangle
		 */
		double intercept_left =
				(slope_on_left == 0) ?
						getlowestX() :
						(double(y - left_offset)) / slope_on_left;
		return intercept_left;
	}

	double get_right_x_intercept(int y) {
		/*
		 * We solve for the intercept for the scan line on the right side of the
		 * triangle
		 */
		double intercept_right =
				(slope_on_right == 0) ?
						gethighestX() :
						(double(y - right_offset)) / slope_on_right;
		return intercept_right;
	}

	bool is_flat_bottom_triangle() {
		/*
		 * Method to determine if the triangle is a flat bottom or a flat
		 * top triangle
		 */
		if (Y[0] == Y[1] || Y[1] == Y[2] || Y[2] == Y[0])
			return true;
		return false;
	}

	void split_triangle(Triangle *t1, Triangle *t2) {
		/*
		 * Method to split an arbitrary traingle
		 */
		double top_vertex[2];
		double bottom_vertex[2];
		double middle_vertex[2];
		int top_index, bottom_index, middle_index;
		top_index = (Y[0] > Y[1] && Y[0] > Y[2]) ? 0 :
					(Y[1] > Y[0] && Y[1] > Y[2]) ? 1 : 2;
		if (Y[(top_index + 1) % 3] > Y[(top_index + 2) % 3]) {
			middle_index = (top_index + 1) % 3;
			bottom_index = (top_index + 2) % 3;
		} else {
			middle_index = (top_index + 2) % 3;
			bottom_index = (top_index + 1) % 3;
		}
		top_vertex[0] = X[top_index];
		top_vertex[1] = Y[top_index];
		bottom_vertex[0] = X[bottom_index];
		bottom_vertex[1] = Y[bottom_index];
		middle_vertex[0] = X[middle_index];
		middle_vertex[1] = Y[middle_index];

		double diff = top_vertex[0] - bottom_vertex[0];
		double slope =
				(diff == 0) ? 0 : (top_vertex[1] - bottom_vertex[1]) / diff;
		double offset = top_vertex[1] - (double(slope * top_vertex[0]));

		double split_vertex[2];
		split_vertex[1] = middle_vertex[1];
		split_vertex[0] =
				(slope == 0) ?
						top_vertex[0] :
						(double(split_vertex[1] - offset)) / slope;

		t1->X[0] = top_vertex[0];
		t1->Y[0] = top_vertex[1];
		t1->X[1] = middle_vertex[0];
		t1->Y[1] = middle_vertex[1];
		t1->X[2] = split_vertex[0];
		t1->Y[2] = split_vertex[1];
		t1->color[0] = color[0];
		t1->color[1] = color[1];
		t1->color[2] = color[2];

		t2->X[0] = bottom_vertex[0];
		t2->Y[0] = bottom_vertex[1];
		t2->X[1] = middle_vertex[0];
		t2->Y[1] = middle_vertex[1];
		t2->X[2] = split_vertex[0];
		t2->Y[2] = split_vertex[1];
		t2->color[0] = color[0];
		t2->color[1] = color[1];
		t2->color[2] = color[2];
	}
};

class Screen {
public:
	unsigned char *buffer;
	int width, height;
	void findPixelAndColor(int x, int y, unsigned char* color) {
		/*
		 * Ensure the pixels to be painted are in the frame.
		 */
		if ((x >= 0 && x < width) && (y >= 0 && y < height)) {
			int buffer_index = (y * 3 * width) + (x * 3);
			if (buffer_index < width * height * 3) {
				buffer[buffer_index++] = color[0];
				buffer[buffer_index++] = color[1];
				buffer[buffer_index] = color[2];
			}
		}
	}
};

std::vector<Triangle> GetTriangles(void) {
	vtkPolyDataReader *rdr = vtkPolyDataReader::New();
	rdr->SetFileName("proj1c_geometry.vtk");
	cerr << "Reading" << endl;
	rdr->Update();
	if (rdr->GetOutput()->GetNumberOfCells() == 0) {
		cerr << "Unable to open file!!" << endl;
		exit (EXIT_FAILURE);
	}
	vtkPolyData *pd = rdr->GetOutput();
	int numTris = pd->GetNumberOfCells();
	vtkPoints *pts = pd->GetPoints();
	vtkCellArray *cells = pd->GetPolys();
	vtkFloatArray *colors = (vtkFloatArray *) pd->GetPointData()->GetArray(
			"color_nodal");
	float *color_ptr = colors->GetPointer(0);
	std::vector<Triangle> tris(numTris);
	vtkIdType npts;
	vtkIdType *ptIds;
	int idx;
	for (idx = 0, cells->InitTraversal(); cells->GetNextCell(npts, ptIds);
			idx++) {
		if (npts != 3) {
			cerr << "Non-triangles!! ???" << endl;
			exit (EXIT_FAILURE);
		}
		tris[idx].X[0] = pts->GetPoint(ptIds[0])[0];
		tris[idx].X[1] = pts->GetPoint(ptIds[1])[0];
		tris[idx].X[2] = pts->GetPoint(ptIds[2])[0];
		tris[idx].Y[0] = pts->GetPoint(ptIds[0])[1];
		tris[idx].Y[1] = pts->GetPoint(ptIds[1])[1];
		tris[idx].Y[2] = pts->GetPoint(ptIds[2])[1];
		tris[idx].color[0] = (unsigned char) color_ptr[4 * ptIds[0] + 0];
		tris[idx].color[1] = (unsigned char) color_ptr[4 * ptIds[0] + 1];
		tris[idx].color[2] = (unsigned char) color_ptr[4 * ptIds[0] + 2];
	}
	cerr << "Done reading" << endl;
	return tris;
}

void scan_line(Triangle *t, Screen *s) {
	/*
	 * Execute scan line algorithm for the current triangle.
	 */
	double y_min = t->getlowestY();
	double y_max = t->gethighestY();
	double x_min = t->getlowestX();
	double x_max = t->gethighestX();
	// Determine the orientation for the triangle
	t->determine_triangle_orientation();
	// Color the pixels that are inside the triangle
	for (int scan_pos = ceil441(y_min); scan_pos <= floor441(y_max);
			scan_pos++) {
		double left_intercept = t->get_left_x_intercept(scan_pos);
		double right_intercept = t->get_right_x_intercept(scan_pos);
		for (int scan_lim = ceil441(left_intercept);
				scan_lim <= floor441(right_intercept); scan_lim++) {
			s->findPixelAndColor(scan_lim, scan_pos, t->color);
		}
	}
}

int main() {
	vtkImageData *image = NewImage(1786, 1344);
	unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0, 0, 0);
	int npixels = 1786 * 1344;
	for (int i = 0; i < npixels * 3; i++)
		buffer[i] = 0;

	std::vector<Triangle> triangles = GetTriangles();

	Screen screen;
	screen.buffer = buffer;
	screen.width = 1786;
	screen.height = 1344;
	for (int vecIndex = 0; vecIndex < triangles.size(); vecIndex++) {
		Triangle t = triangles[vecIndex];
		if (t.is_flat_bottom_triangle()) {
			scan_line(&t, &screen);
		} else {
			Triangle t1, t2;
			t.split_triangle(&t1, &t2);
			scan_line(&t1, &screen);
			scan_line(&t2, &screen);
		}
	}
	WriteImage(image, "allTriangles");
}
