#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

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

	int getlowestY() {
		return ((Y[0] < Y[1]) ? Y[0] : ((Y[1] < Y[2]) ? Y[1] : Y[2]));
	}

	int gethighestY() {
		return ((Y[0] > Y[1]) ? Y[0] : ((Y[1] > Y[2]) ? Y[1] : Y[2]));
	}

	int getlowestX() {
		return ((X[0] < X[1]) ? X[0] : ((X[1] < X[2]) ? X[1] : X[2]));
	}

	int gethighestX() {
		return ((X[0] > X[1]) ? X[0] : ((X[1] > X[2]) ? X[1] : X[2]));
	}

	void determine_triangle_orientation() {
		/*
		 * We determine the vertexes that represent the base and the other
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
		int diff;
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
		 * We solve for the intercept for the scan line on the left side of the triangle
		 */
		double intercept_left =
				(slope_on_left == 0) ?
						getlowestX() :
						(double(y - left_offset)) / slope_on_left;
		return intercept_left;
	}

	double get_right_x_intercept(int y) {
		/*
		 * We solve for the intercept for the scan line on the right side of the triangle
		 */
		double intercept_right =
				(slope_on_right == 0) ?
						gethighestX() :
						(double(y - right_offset)) / slope_on_right;
		return intercept_right;
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
			int buffer_index = y * 3000 + x * 3;
			if (buffer_index < width * height * 3) {
				buffer[buffer_index++] = color[0];
				buffer[buffer_index++] = color[1];
				buffer[buffer_index] = color[2];
			}
		}
	}
};

std::vector<Triangle> GetTriangles(void) {
	std::vector<Triangle> rv(100);

	unsigned char colors[6][3] = { { 255, 128, 0 }, { 255, 0, 127 }, { 0, 204,
			204 }, { 76, 153, 0 }, { 255, 204, 204 }, { 204, 204, 0 } };
	for (int i = 0; i < 100; i++) {
		int idxI = i % 10;
		int posI = idxI * 100;
		int idxJ = i / 10;
		int posJ = idxJ * 100;
		int firstPt = (i % 3);
		rv[i].X[firstPt] = posI;
		if (i == 50)
			rv[i].X[firstPt] = -10;
		rv[i].Y[firstPt] = posJ;
		rv[i].X[(firstPt + 1) % 3] = posI + 99;
		rv[i].Y[(firstPt + 1) % 3] = posJ;
		rv[i].X[(firstPt + 2) % 3] = posI + i;
		rv[i].Y[(firstPt + 2) % 3] = posJ + 10 * (idxJ + 1);
		if (i == 95)
			rv[i].Y[(firstPt + 2) % 3] = 1050;
		rv[i].color[0] = colors[i % 6][0];
		rv[i].color[1] = colors[i % 6][1];
		rv[i].color[2] = colors[i % 6][2];
	}

	return rv;
}

void scan_line(Triangle *t, Screen *s) {
	/*
	 * Execute scan line algorithm for the current triangle.
	 */
	int y_min = t->getlowestY();
	int y_max = t->gethighestY();
	int x_min = t->getlowestX();
	int x_max = t->gethighestX();
	// Determine the orientation for the triangle
	t->determine_triangle_orientation();
	// Color the pixels that are inside the triangle
	for (int scan_pos = y_min; scan_pos <= y_max; scan_pos++) {
		double left_intercept = t->get_left_x_intercept(scan_pos);
		double right_intercept = t->get_right_x_intercept(scan_pos);
		for (int scan_lim = ceil441(left_intercept);
				scan_lim <= floor441(right_intercept); scan_lim++) {
			s->findPixelAndColor(scan_lim, scan_pos, t->color);
		}
	}
}

int main() {
	vtkImageData *image = NewImage(1000, 1000);
	unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0, 0, 0);
	int npixels = 1000 * 1000;
	for (int i = 0; i < npixels * 3; i++)
		buffer[i] = 0;

	std::vector<Triangle> triangles = GetTriangles();

	Screen screen;
	screen.buffer = buffer;
	screen.width = 1000;
	screen.height = 1000;

	for (int vecIndex = 0; vecIndex < 100; vecIndex++) {
		Triangle t = triangles[vecIndex];
		scan_line(&t, &screen);
	}
	WriteImage(image, "allTriangles");
}
