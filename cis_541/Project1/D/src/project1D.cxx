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
	double Z[3];
	double colors[3][3];
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
	int offset_index, left_index, right_index;

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
		offset_index = ((Y[0] == Y[1]) ? 2 : ((Y[1] == Y[2]) ? 0 : 1));
		offset_vertex[0] = X[offset_index];
		offset_vertex[1] = Y[offset_index];
		if (X[(offset_index + 1) % 3] < X[(offset_index + 2) % 3]) {
			left_index = (offset_index + 1) % 3;
			right_index = (offset_index + 2) % 3;
		} else {
			left_index = (offset_index + 2) % 3;
			right_index = (offset_index + 1) % 3;
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
		//TODO : Interpolate Z and Colors for the new vertex;
		double z_split = interpolate(top_vertex[1], bottom_vertex[1],
				Z[top_index], Z[bottom_index], split_vertex[1]);
		double split_colors[3];
		split_colors[0] = interpolate(top_vertex[1], bottom_vertex[1],
				colors[top_index][0], colors[bottom_index][0], split_vertex[1]);
		split_colors[1] = interpolate(top_vertex[1], bottom_vertex[1],
				colors[top_index][1], colors[bottom_index][1], split_vertex[1]);
		split_colors[2] = interpolate(top_vertex[1], bottom_vertex[1],
				colors[top_index][2], colors[bottom_index][2], split_vertex[1]);

		t1->X[0] = top_vertex[0];
		t1->Y[0] = top_vertex[1];
		t1->X[1] = middle_vertex[0];
		t1->Y[1] = middle_vertex[1];
		t1->X[2] = split_vertex[0];
		t1->Y[2] = split_vertex[1];
		t1->Z[0] = Z[top_index];
		t1->Z[1] = Z[middle_index];
		t1->Z[2] = z_split;

		t1->colors[0][0] = colors[top_index][0];
		t1->colors[0][1] = colors[top_index][1];
		t1->colors[0][2] = colors[top_index][2];
		t1->colors[1][0] = colors[middle_index][0];
		t1->colors[1][1] = colors[middle_index][1];
		t1->colors[1][2] = colors[middle_index][2];
		t1->colors[2][0] = split_colors[0];
		t1->colors[2][1] = split_colors[1];
		t1->colors[2][2] = split_colors[2];

		t2->X[0] = bottom_vertex[0];
		t2->Y[0] = bottom_vertex[1];
		t2->X[1] = middle_vertex[0];
		t2->Y[1] = middle_vertex[1];
		t2->X[2] = split_vertex[0];
		t2->Y[2] = split_vertex[1];
		t2->Z[0] = Z[bottom_index];
		t2->Z[1] = Z[middle_index];
		t2->Z[2] = z_split;

		t2->colors[0][0] = colors[bottom_index][0];
		t2->colors[0][1] = colors[bottom_index][1];
		t2->colors[0][2] = colors[bottom_index][2];
		t2->colors[1][0] = colors[middle_index][0];
		t2->colors[1][1] = colors[middle_index][1];
		t2->colors[1][2] = colors[middle_index][2];
		t2->colors[2][0] = split_colors[0];
		t2->colors[2][1] = split_colors[1];
		t2->colors[2][2] = split_colors[2];
	}

	double interpolate(double point_1, double point_2, double value_1,
			double value_2, double quest_point) {
		double diff = point_2 - point_1;
		double change_prop =
				(diff != 0) ? (quest_point - point_1) / (point_2 - point_1) : 0;
		double value_at_quest = value_1 + (change_prop * (value_2 - value_1));
		return value_at_quest;
	}
};

class Screen {
public:
	unsigned char *buffer;
	double *depth_buffer;
	int width, height;

	//TODO : make colors proportional
	void findPixelAndColor(int x, int y, double *color, double current_depth) {
		/*
		 * Ensure the pixels to be painted are in the frame.
		 */
		if ((x >= 0 && x < width) && (y >= 0 && y < height)) {
			int buffer_index = (y * 3 * width) + (x * 3);
			int depth_buffer_index = y * width + x;
			if (buffer_index < width * height * 3
					&& current_depth >= depth_buffer[depth_buffer_index]) {
				buffer[buffer_index++] = ceil441(color[0] * 255);
				buffer[buffer_index++] = ceil441(color[1] * 255);
				buffer[buffer_index] = ceil441(color[2] * 255);
				depth_buffer[depth_buffer_index] = current_depth;
			}
		}
	}
};

std::vector<Triangle> GetTriangles(void) {
	vtkPolyDataReader *rdr = vtkPolyDataReader::New();
	rdr->SetFileName("proj1d_geometry.vtk");
	cerr << "Reading" << endl;
	rdr->Update();
	cerr << "Done reading" << endl;
	if (rdr->GetOutput()->GetNumberOfCells() == 0) {
		cerr << "Unable to open file!!" << endl;
		exit (EXIT_FAILURE);
	}
	vtkPolyData *pd = rdr->GetOutput();
	int numTris = pd->GetNumberOfCells();
	vtkPoints *pts = pd->GetPoints();
	vtkCellArray *cells = pd->GetPolys();
	vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray(
			"hardyglobal");
	float *color_ptr = var->GetPointer(0);
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
		tris[idx].Z[0] = pts->GetPoint(ptIds[0])[2];
		tris[idx].Z[1] = pts->GetPoint(ptIds[1])[2];
		tris[idx].Z[2] = pts->GetPoint(ptIds[2])[2];
		// 1->2 interpolate between light blue, dark blue
		// 2->2.5 interpolate between dark blue, cyan
		// 2.5->3 interpolate between cyan, green
		// 3->3.5 interpolate between green, yellow
		// 3.5->4 interpolate between yellow, orange
		// 4->5 interpolate between orange, brick
		// 5->6 interpolate between brick, salmon
		double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
		double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
		unsigned char RGB[8][3] = { { 71, 71, 219 }, { 0, 0, 91 },
				{ 0, 255, 255 }, { 0, 128, 0 }, { 255, 255, 0 }, { 255, 96, 0 },
				{ 107, 0, 0 }, { 224, 76, 76 } };
		for (int j = 0; j < 3; j++) {
			float val = color_ptr[ptIds[j]];
			int r;
			for (r = 0; r < 7; r++) {
				if (mins[r] <= val && val < maxs[r])
					break;
			}
			if (r == 7) {
				cerr << "Could not interpolate color for " << val << endl;
				exit (EXIT_FAILURE);
			}
			double proportion = (val - mins[r]) / (maxs[r] - mins[r]);
			tris[idx].colors[j][0] = (RGB[r][0]
					+ proportion * (RGB[r + 1][0] - RGB[r][0])) / 255.0;
			tris[idx].colors[j][1] = (RGB[r][1]
					+ proportion * (RGB[r + 1][1] - RGB[r][1])) / 255.0;
			tris[idx].colors[j][2] = (RGB[r][2]
					+ proportion * (RGB[r + 1][2] - RGB[r][2])) / 255.0;
		}
	}
	return tris;
}

void calculate_color_for_pixel(Triangle *t, double left_intercept,
		double right_intercept, double current_x_pos, double current_y_pos, double *current_color) {
	double color_at_left_intercept[3];
	double color_at_right_intercept[3];

	color_at_left_intercept[0] = t->interpolate(t->offset_vertex[1],
			t->left_vertex[1], t->colors[t->offset_index][0],
			t->colors[t->left_index][0], current_y_pos);
	color_at_left_intercept[1] = t->interpolate(t->offset_vertex[1],
			t->left_vertex[1], t->colors[t->offset_index][1],
			t->colors[t->left_index][1], current_y_pos);
	color_at_left_intercept[2] = t->interpolate(t->offset_vertex[1],
			t->left_vertex[1], t->colors[t->offset_index][2],
			t->colors[t->left_index][2], current_y_pos);

	color_at_right_intercept[0] = t->interpolate(t->offset_vertex[1],
			t->right_vertex[1], t->colors[t->offset_index][0],
			t->colors[t->right_index][0], current_y_pos);
	color_at_right_intercept[1] = t->interpolate(t->offset_vertex[1],
			t->right_vertex[1], t->colors[t->offset_index][1],
			t->colors[t->right_index][1], current_y_pos);
	color_at_right_intercept[2] = t->interpolate(t->offset_vertex[1],
			t->right_vertex[1], t->colors[t->offset_index][2],
			t->colors[t->right_index][2], current_y_pos);

	/*if (current_x_pos == ceil441(left_intercept)) {
		current_color[0] = color_at_left_intercept[0];
		current_color[1] = color_at_left_intercept[1];
		current_color[2] = color_at_left_intercept[2];
	} else if (current_x_pos == floor441(right_intercept)) {
		current_color[0] = color_at_right_intercept[0];
		current_color[1] = color_at_right_intercept[1];
		current_color[2] = color_at_right_intercept[2];
	} else {*/
		current_color[0] = t->interpolate(left_intercept, right_intercept,
				color_at_left_intercept[0], color_at_right_intercept[0],
				current_x_pos);
		current_color[1] = t->interpolate(left_intercept, right_intercept,
				color_at_left_intercept[1], color_at_right_intercept[1],
				current_x_pos);
		current_color[2] = t->interpolate(left_intercept, right_intercept,
				color_at_left_intercept[2], color_at_right_intercept[2],
				current_x_pos);
	//}
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
		double z_left_intercept = t->interpolate(t->offset_vertex[1],t->left_vertex[1], t->Z[t->offset_index], t->Z[t->left_index], scan_pos);
		double z_right_intercept = t->interpolate(t->offset_vertex[1],t->right_vertex[1], t->Z[t->offset_index], t->Z[t->right_index], scan_pos);
		for (int scan_lim = ceil441(left_intercept);
				scan_lim <= floor441(right_intercept); scan_lim++) {
			double z_current;
			/*if (scan_lim == ceil441(left_intercept))
				z_current = z_left_intercept;
			else if (scan_lim == floor441(right_intercept))
				z_current = z_right_intercept;
			else*/
				z_current = t->interpolate(left_intercept, right_intercept,
						z_left_intercept, z_right_intercept, scan_lim);

			double color[3] = { 0, 0, 0 };
			calculate_color_for_pixel(t, left_intercept, right_intercept,
					scan_lim, scan_pos, color);
			s->findPixelAndColor(scan_lim, scan_pos, color, z_current);
		}
	}
}

int main() {
	vtkImageData *image = NewImage(1000, 1000);
	unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0, 0, 0);
	int npixels = 1000 * 1000;
	for (int i = 0; i < npixels * 3; i++)
		buffer[i] = 0;
	double depth_buffer[npixels];
	for (int i = 0; i < npixels; i++)
		depth_buffer[i] = -1;
	std::vector<Triangle> triangles = GetTriangles();

	Screen screen;
	screen.buffer = buffer;
	screen.depth_buffer = depth_buffer;
	screen.width = 1000;
	screen.height = 1000;
	for (int vecIndex = 0; vecIndex < triangles.size(); vecIndex++) {
		Triangle t = triangles[vecIndex];
		/*cout << "\nFor Triangle : " << vecIndex + 1 << "\n" << t.X[0] << " " << t.X[1] << " " << t.X[2];
		 cout << "\n" << t.Y[0] << " " << t.Y[1] << " " << t.Y[2];
		 cout << "\n" << t.Z[0] << " " << t.Z[1] << " " << t.Z[2];
		 cout << "\n" << t.colors[0][0] << " " << t.colors[0][1] << " " << t.colors[0][2];
		 cout << "\n" << t.colors[1][0] << " " << t.colors[1][1] << " " << t.colors[1][2];
		 cout << "\n" << t.colors[2][0] << " " << t.colors[2][1] << " " << t.colors[2][2];*/
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
