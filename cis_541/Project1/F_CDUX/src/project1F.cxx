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
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetWriter.h>
#include <string.h>
#include <cmath>
#include <algorithm>
#include <sstream>

using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;
using std::pow;
using std::tan;

double ceil441(double f) {
	return ceil(f - 0.00001);
}

double floor441(double f) {
	return floor(f + 0.00001);
}

/**
 * This function interpolates the value at quest_point based on
 * values value_1 and value_2 at points point_1 and point_2 respectively
 */
double interpolate(double point_1, double point_2, double value_1,
		double value_2, double quest_point) {
	double diff = point_2 - point_1;
	double change_prop =
			(diff != 0) ? (quest_point - point_1) / (point_2 - point_1) : 0;
	double value_at_quest = value_1 + (change_prop * (value_2 - value_1));
	return value_at_quest;
}

void normalize_vector(double* quest_normal) {
	double norm = sqrt(
			(quest_normal[0] * quest_normal[0])
					+ (quest_normal[1] * quest_normal[1])
					+ (quest_normal[2] * quest_normal[2]));
	quest_normal[0] = quest_normal[0] / norm;
	quest_normal[1] = quest_normal[1] / norm;
	quest_normal[2] = quest_normal[2] / norm;
}

double dot_product(double* vector_1, double* vector_2) {
	return (vector_1[0] * vector_2[0]) + (vector_1[1] * vector_2[1])
			+ (vector_1[2] * vector_2[2]);
}

void cross_product(double* vector_1, double* vector_2, double *cross_vec) {
 	double product[3] = { (vector_1[1]*vector_2[2]) - (vector_1[2]*vector_2[1]),
												(vector_1[2]*vector_2[0]) - (vector_1[0]*vector_2[2]),
												(vector_1[0]*vector_2[1]) - (vector_1[1]*vector_2[0])};
  memcpy(cross_vec, product, 3 * sizeof(double));
}

void print_vector(double *print_vector) {
	cout << "\n{" << print_vector[0] << ", " << print_vector[1] << ", " << print_vector[2] << "}" << endl;
}

double cot(double angle) {
	return (1/tan(angle));
}

/**
 * This function interpolates over the vectors normal_1 and normal_2 to yield
 * the vector quest_normal for the point (quest_point[0],quest_point[1]) in space.
 */
void interpolate_vector(double point_1, double point_2, double* normal_1,
		double* normal_2, double quest_point, double* quest_normal) {
	double diff = point_2 - point_1;
	double proportion =
			(diff != 0) ? (quest_point - point_1) / (point_2 - point_1) : 0;
	double diff_vector[3] = { normal_2[0] - normal_1[0], normal_2[1]
			- normal_1[1], normal_2[2] - normal_1[2] };
	diff_vector[0] = proportion * diff_vector[0];
	diff_vector[1] = proportion * diff_vector[1];
	diff_vector[2] = proportion * diff_vector[2];
	quest_normal[0] = normal_1[0] + diff_vector[0];
	quest_normal[1] = normal_1[1] + diff_vector[1];
	quest_normal[2] = normal_1[2] + diff_vector[2];
	normalize_vector(quest_normal);
}

class Screen {
public:
	unsigned char *buffer;
	double *depth_buffer;
	int width, height;

	void find_pixel_and_color(int x, int y, double *color, double current_depth/*,
			double shading_amount*/) {
		/*
		 * Ensure the pixels to be painted are in the frame.
		 */
		double shading_amount = 1;
		if ((x >= 0 && x < width) && (y >= 0 && y < height)) {
			int buffer_index = (y * 3 * width) + (x * 3);
			int depth_buffer_index = y * width + x;
			if (buffer_index < width * height * 3
          && (current_depth > -1 && current_depth < 1)
          && current_depth >= depth_buffer[depth_buffer_index]) {
				/*double shading_amount = calculate_phong_shading(lp,
						view_direction, current_normal);*/
				buffer[buffer_index++] = min(
						ceil441((shading_amount * color[0]) * 255),
						(double) 255);
				buffer[buffer_index++] = min(
						ceil441((shading_amount * color[1]) * 255),
						(double) 255);
				buffer[buffer_index] = min(
						ceil441((shading_amount * color[2]) * 255),
						(double) 255);
				depth_buffer[depth_buffer_index] = current_depth;
			}
		}
	}

	/*
	 * This method calculates the color for a certain pixel in a trangle given
	 * its intercepts for a scanline and the colors at the intercepts.
	 */
	void calculate_color_for_pixel(double left_intercept,
			double right_intercept, double current_x_pos,
			double *color_at_left_intercept, double *color_at_right_intercept,
			double *current_color) {

		current_color[0] = interpolate(left_intercept, right_intercept,
				color_at_left_intercept[0], color_at_right_intercept[0],
				current_x_pos);
		current_color[1] = interpolate(left_intercept, right_intercept,
				color_at_left_intercept[1], color_at_right_intercept[1],
				current_x_pos);
		current_color[2] = interpolate(left_intercept, right_intercept,
				color_at_left_intercept[2], color_at_right_intercept[2],
				current_x_pos);
	}
};

//Defining all data structures
class Matrix
{
  public:
    double          A[4][4];

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

void
Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix
Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}

void
Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
             + ptIn[1]*A[1][0]
             + ptIn[2]*A[2][0]
             + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
             + ptIn[1]*A[1][1]
             + ptIn[2]*A[2][1]
             + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
             + ptIn[1]*A[1][2]
             + ptIn[2]*A[2][2]
             + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
             + ptIn[1]*A[1][3]
             + ptIn[2]*A[2][3]
             + ptIn[3]*A[3][3];
}

struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 5.3;
         alpha = 7.5;
    };


    double lightDir[3];  // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};

LightingParameters lp;

class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];

    Matrix CameraTransform(void) {
    	double z_vector[3] = {position[0] - focus[0],
    												position[1] - focus[1],
    												position[2] - focus[2]};
    	normalize_vector(z_vector);
    	double x_vector[3],y_vector[3];
  		cross_product(up, z_vector, x_vector);
  		normalize_vector(x_vector);
  		cross_product(z_vector, x_vector, y_vector);
  		normalize_vector(y_vector);
    	double frame_vec[3] = {
    		0 - position[0],
    		0 - position[1],
    		0 - position[2]
    	};
    	//normalize_vector(frame_vec);
    	print_vector(position);
    	print_vector(x_vector);
    	print_vector(y_vector);
    	print_vector(z_vector);
    	print_vector(frame_vec);
      double cartesian_x[3] = {1,0,0};
    	double cartesian_y[3] = {0,1,0};
    	double cartesian_z[3] = {0,0,1};
    	double camera_transform[4][4] = {
    		dot_product(x_vector,cartesian_x), dot_product(y_vector, cartesian_x), dot_product(z_vector, cartesian_x), 0,
    		dot_product(x_vector,cartesian_y), dot_product(y_vector, cartesian_y), dot_product(z_vector, cartesian_y), 0,
    		dot_product(x_vector,cartesian_z), dot_product(y_vector, cartesian_z), dot_product(z_vector, cartesian_z), 0,
    		dot_product(x_vector,frame_vec), dot_product(y_vector, frame_vec), dot_product(z_vector, frame_vec), 1
    	};
    	Matrix m;
    	memcpy(m.A, camera_transform, 16*sizeof(double));
    	return m;
    }

    Matrix ViewTransform(void) {
    	double view_transform[4][4] = {
    		cot(angle/2), 0, 0 , 0,
    		0, cot(angle/2), 0 , 0,
    		0, 0 , (far + near)/(far - near), -1,
    		0, 0 , (2*far*near)/(far - near), 0
    	};
    	Matrix m;
    	memcpy(m.A, view_transform, 16*sizeof(double));
    	return m;
    }

    Matrix DeviceTransform(Screen screen) {
      double device_transform[4][4] = {
        screen.width/2, 0, 0, 0,
        0, screen.height/2, 0, 0,
        0, 0, 1, 0,
        screen.width/2, screen.height/2, 0, 1
      };
      Matrix m;
      memcpy(m.A, device_transform, 16*sizeof(double));
      return m;
    };

};


double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
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
	double normals[3][3];
  double shading[3];
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

		/*
		 * Calculate for the value of Z buffer and colors for the new vertex
		 */
		double z_split = interpolate(top_vertex[1], bottom_vertex[1],
				Z[top_index], Z[bottom_index], split_vertex[1]);
    /*double shade_split = interpolate(top_vertex[1], bottom_vertex[1],
				shading[top_index], shading[bottom_index], split_vertex[1]);*/

		double split_colors[3];
		split_colors[0] = interpolate(top_vertex[1], bottom_vertex[1],
				colors[top_index][0], colors[bottom_index][0], split_vertex[1]);
		split_colors[1] = interpolate(top_vertex[1], bottom_vertex[1],
				colors[top_index][1], colors[bottom_index][1], split_vertex[1]);
		split_colors[2] = interpolate(top_vertex[1], bottom_vertex[1],
				colors[top_index][2], colors[bottom_index][2], split_vertex[1]);

		/*double split_normal[3];
		interpolate_vector(top_vertex[1], bottom_vertex[1], normals[top_index],
				normals[bottom_index], split_vertex[1], split_normal);*/

		t1->X[0] = top_vertex[0];
		t1->Y[0] = top_vertex[1];
		t1->X[1] = middle_vertex[0];
		t1->Y[1] = middle_vertex[1];
		t1->X[2] = split_vertex[0];
		t1->Y[2] = split_vertex[1];
    t1->Z[0] = Z[top_index];
		t1->Z[1] = Z[middle_index];
		t1->Z[2] = z_split;
    /*t1->shading[0] = shading[top_index];
		t1->shading[1] = shading[middle_index];
		t1->shading[2] = shade_split;*/
		/*memcpy(t1->normals[0], this->normals[top_index], 3 * sizeof(double));
		memcpy(t1->normals[1], this->normals[middle_index], 3 * sizeof(double));
		memcpy(t1->normals[2], split_normal, 3 * sizeof(double));*/

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
    /*t2->shading[0] = shading[bottom_index];
		t2->shading[1] = shading[middle_index];
		t2->shading[2] = shade_split;*/
		/*memcpy(t2->normals[0], this->normals[bottom_index], 3 * sizeof(double));
		memcpy(t2->normals[1], this->normals[middle_index], 3 * sizeof(double));
		memcpy(t2->normals[2], split_normal, 3 * sizeof(double));*/

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

	/*
	 * This method helps to calculate the colors at extremes of the scanline
	 * corresponding to quest_point
	 */
	void calculate_color_for_scanline_extremes(double quest_point,
			double *color_at_left_intercept, double *color_at_right_intercept) {

		color_at_left_intercept[0] = interpolate(offset_vertex[1],
				left_vertex[1], colors[offset_index][0], colors[left_index][0],
				quest_point);
		color_at_left_intercept[1] = interpolate(offset_vertex[1],
				left_vertex[1], colors[offset_index][1], colors[left_index][1],
				quest_point);
		color_at_left_intercept[2] = interpolate(offset_vertex[1],
				left_vertex[1], colors[offset_index][2], colors[left_index][2],
				quest_point);

		color_at_right_intercept[0] = interpolate(offset_vertex[1],
				right_vertex[1], colors[offset_index][0],
				colors[right_index][0], quest_point);
		color_at_right_intercept[1] = interpolate(offset_vertex[1],
				right_vertex[1], colors[offset_index][1],
				colors[right_index][1], quest_point);
		color_at_right_intercept[2] = interpolate(offset_vertex[1],
				right_vertex[1], colors[offset_index][2],
				colors[right_index][2], quest_point);
	}

};

double calculate_for_specular_lighting(LightingParameters lp, double *V,
		double* N) {
	/* R = 2*(L . N)*N - L
	 * return ==> V . R
	 */
	double two_L_dot_N = 2 * dot_product(lp.lightDir, N);
	double two_L_dot_N_N[3] = { two_L_dot_N * N[0], two_L_dot_N * N[1],
			two_L_dot_N * N[2] };
	double R[3] = { two_L_dot_N_N[0] - lp.lightDir[0], two_L_dot_N_N[1] - lp.lightDir[1],
			two_L_dot_N_N[2] - lp.lightDir[2] };
	return dot_product(R, V);
}

double calculate_phong_shading(LightingParameters lp, double *view_direction,
		double *normal) {
	double diffuse_component = dot_product(lp.lightDir, normal);
  double specular_component = pow(calculate_for_specular_lighting(lp, view_direction, normal),lp.alpha);
  if(specular_component < 0 || specular_component != specular_component) {
      specular_component = 0;
  }
	double shading_amount = lp.Ka + lp.Kd * abs(diffuse_component)
			+ lp.Ks * specular_component;
}

<<<<<<< HEAD:cis_541/Project1/F_CDUX/src/project1F.cxx
std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("hardyglobal.0.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    /*vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);*/
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
        /*tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];*/
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
        /*tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];*/
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
        /*tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];*/

        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 },
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 }
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
=======
class Screen {
public:
	unsigned char *buffer;
	double *depth_buffer;
	int width, height;

	void find_pixel_and_color(int x, int y, double *color, double current_depth/*,
			double shading_amount*/) {
		/*
		 * Ensure the pixels to be painted are in the frame.
		 */
		if ((x >= 0 && x < width) && (y >= 0 && y < height)) {
			int buffer_index = (y * 3 * width) + (x * 3);
			int depth_buffer_index = y * width + x;
			if (buffer_index < width * height * 3
          && (current_depth > -1 && current_depth < 1)
          && current_depth >= depth_buffer[depth_buffer_index]) {
				/*double shading_amount = calculate_phong_shading(lp,
						view_direction, current_normal);*/
				buffer[buffer_index++] = min(
						ceil441((/*shading_amount * */color[0]) * 255),
						(double) 255);
				buffer[buffer_index++] = min(
						ceil441((/*shading_amount * */color[1]) * 255),
						(double) 255);
				buffer[buffer_index] = min(
						ceil441((/*shading_amount * */color[2]) * 255),
						(double) 255);
				depth_buffer[depth_buffer_index] = current_depth;
			}
		}
	}

	/*
	 * This method calculates the color for a certain pixel in a trangle given
	 * its intercepts for a scanline and the colors at the intercepts.
	 */
	void calculate_color_for_pixel(double left_intercept,
			double right_intercept, double current_x_pos,
			double *color_at_left_intercept, double *color_at_right_intercept,
			double *current_color) {

		current_color[0] = interpolate(left_intercept, right_intercept,
				color_at_left_intercept[0], color_at_right_intercept[0],
				current_x_pos);
		current_color[1] = interpolate(left_intercept, right_intercept,
				color_at_left_intercept[1], color_at_right_intercept[1],
				current_x_pos);
		current_color[2] = interpolate(left_intercept, right_intercept,
				color_at_left_intercept[2], color_at_right_intercept[2],
				current_x_pos);
	}
};

std::vector<Triangle>
GetTriangles(void)
{
    std::vector<Triangle> triangles(0);
    int index = 0;
    for(int i=0;i<63;i++) {
      vtkPolyDataReader *rdr = vtkPolyDataReader::New();
      /*std::stringstream vtk_file;
      vtk_file << "hardyglobal." << i << ".vtk";
      rdr->SetFileName(vtk_file.c_str());*/
      std::ostringstream stream;
      stream << "hardyglobal." << i << ".vtk";
      std::string str =  stream.str();
      const char* vtk_file = str.c_str();
      rdr->SetFileName(vtk_file);

      cerr << "Reading" << endl;
      rdr->Update();
      cerr << "Done reading" << endl;
      if (rdr->GetOutput()->GetNumberOfCells() == 0)
      {
          cerr << "Unable to open file!!" << endl;
          exit(EXIT_FAILURE);
      }
      vtkPolyData *pd = rdr->GetOutput();
      /*
      vtkDataSetWriter *writer = vtkDataSetWriter::New();
      writer->SetInput(pd);
      writer->SetFileName("hrc.vtk");
      writer->Write();
       */
      int numTris = pd->GetNumberOfCells();
      vtkPoints *pts = pd->GetPoints();
      vtkCellArray *cells = pd->GetPolys();
      vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
      double *color_ptr = var->GetPointer(0);
      //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
      //float *color_ptr = var->GetPointer(0);
      /*vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
      float *normals = n->GetPointer(0);*/
      triangles.resize(triangles.size() + numTris);
      vtkIdType npts;
      vtkIdType *ptIds;
      int idx;
      for (idx = index, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++, index++)
      {
          if (npts != 3)
          {
              cerr << "Non-triangles!! ???" << endl;
              exit(EXIT_FAILURE);
          }
          double *pt = NULL;
          pt = pts->GetPoint(ptIds[0]);
          triangles[idx].X[0] = pt[0];
          triangles[idx].Y[0] = pt[1];
          triangles[idx].Z[0] = pt[2];
          /*tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
          tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
          tris[idx].normals[0][2] = normals[3*ptIds[0]+2];*/
          pt = pts->GetPoint(ptIds[1]);
          triangles[idx].X[1] = pt[0];
          triangles[idx].Y[1] = pt[1];
          triangles[idx].Z[1] = pt[2];
          /*tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
          tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
          tris[idx].normals[1][2] = normals[3*ptIds[1]+2];*/
          pt = pts->GetPoint(ptIds[2]);
          triangles[idx].X[2] = pt[0];
          triangles[idx].Y[2] = pt[1];
          triangles[idx].Z[2] = pt[2];
          /*tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
          tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
          tris[idx].normals[2][2] = normals[3*ptIds[2]+2];*/

          // 1->2 interpolate between light blue, dark blue
          // 2->2.5 interpolate between dark blue, cyan
          // 2.5->3 interpolate between cyan, green
          // 3->3.5 interpolate between green, yellow
          // 3.5->4 interpolate between yellow, orange
          // 4->5 interpolate between orange, brick
          // 5->6 interpolate between brick, salmon
          double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
          double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
          unsigned char RGB[8][3] = { { 71, 71, 219 },
                                      { 0, 0, 91 },
                                      { 0, 255, 255 },
                                      { 0, 128, 0 },
                                      { 255, 255, 0 },
                                      { 255, 96, 0 },
                                      { 107, 0, 0 },
                                      { 224, 76, 76 }
                                    };
          for (int j = 0 ; j < 3 ; j++)
          {
              float val = color_ptr[ptIds[j]];
              int r;
              for (r = 0 ; r < 7 ; r++)
              {
                  if (mins[r] <= val && val < maxs[r])
                      break;
              }
              if (r == 7)
              {
                  /*cerr << "Could not interpolate color for " << val << endl;
                  exit(EXIT_FAILURE);*/
                  r = ((int)val)%7;
              }
              double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
              triangles[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
              triangles[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
              triangles[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
          }
      }
>>>>>>> 55cca6a5ab9ec6227a3200879db076a7779d7040:cis_541/Project1/F-CDUX/src/project1Fcdux.cxx
    }
    return triangles;
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
	for (int current_y = ceil441(y_min); current_y <= floor441(y_max);
			current_y++) {
		double left_intercept = t->get_left_x_intercept(current_y);
		double right_intercept = t->get_right_x_intercept(current_y);
		double z_left_intercept = interpolate(t->offset_vertex[1],
				t->left_vertex[1], t->Z[t->offset_index], t->Z[t->left_index],
				current_y);
		double z_right_intercept = interpolate(t->offset_vertex[1],
				t->right_vertex[1], t->Z[t->offset_index], t->Z[t->right_index],
				current_y);
		double color_at_left_intercept[3], color_at_right_intercept[3];
		t->calculate_color_for_scanline_extremes(current_y,
				color_at_left_intercept, color_at_right_intercept);

    /*double shading_left_intercept = interpolate(t->offset_vertex[1],
				t->left_vertex[1], t->shading[t->offset_index], t->shading[t->left_index],
				current_y);
		double shading_right_intercept = interpolate(t->offset_vertex[1],
				t->right_vertex[1], t->shading[t->offset_index], t->shading[t->right_index],
				current_y);*/

		for (int current_x = ceil441(left_intercept);
				current_x <= floor441(right_intercept); current_x++) {
			double current_z = interpolate(left_intercept, right_intercept,
					z_left_intercept, z_right_intercept, current_x);
      /*double current_shading = interpolate(left_intercept, right_intercept,
					shading_left_intercept, shading_right_intercept, current_x);*/
			double color_for_current_pixel[3] = { 0, 0, 0 };
			s->calculate_color_for_pixel(left_intercept, right_intercept,
					current_x, color_at_left_intercept,
					color_at_right_intercept, color_for_current_pixel);

			/*double current_normal[3];
			interpolate_vector(left_intercept, right_intercept, normal_on_left,
					normal_on_right, current_x, current_normal);*/
			s->find_pixel_and_color(current_x, current_y,
					color_for_current_pixel, current_z/*, current_shading*/);
		}
	}
}

void transformTriangle(Triangle *t,Matrix composite, Camera camera) {
  for(int i =0;i < 3;i++) {
    double view_dir[3] = {
      t->X[i] - camera.position[0],
      t->Y[i] - camera.position[1],
      t->Z[i] - camera.position[2]
    };
    normalize_vector(view_dir);
    t->shading[i] = calculate_phong_shading(lp, view_dir, t->normals[i]);
    double current_quadro[4] = {
      t->X[i],t->Y[i],t->Z[i], 1
    };
    double transformed_vertex[4];
    composite.TransformPoint(current_quadro, transformed_vertex);
    if(transformed_vertex[3] != 1) {
      for(int j=0;j < 4;j++)
        transformed_vertex[j] = transformed_vertex[j] / transformed_vertex[3];
    }
    t->X[i] = transformed_vertex[0];
    t->Y[i] = transformed_vertex[1];
    t->Z[i] = transformed_vertex[2];
  }
}

void print_triangle(Triangle t) {
  cout << "\nPrinting Triangle" << endl;
  cout << t.X[0] << ", " <<  t.Y[0] << ", " <<  t.Z[0] << endl;
  cout << t.X[1] << ", " <<  t.Y[1] << ", " <<  t.Z[1] << endl;
  cout << t.X[2] << ", " <<  t.Y[2] << ", " <<  t.Z[2] << endl;
}

Matrix get_total_transform_matrix(Matrix camera_transform, Matrix view_transform, Matrix device_transform) {
  Matrix composite;
  composite = Matrix::ComposeMatrices(camera_transform, view_transform);
  composite = Matrix::ComposeMatrices(composite, device_transform);
  return composite;
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

  Camera camera = GetCamera(0,1000);

  Matrix camera_transform = camera.CameraTransform();
	camera_transform.Print(std::cout);
	Matrix view_transform = camera.ViewTransform();
  view_transform.Print(std::cout);
  Matrix device_transform = camera.DeviceTransform(screen);
  device_transform.Print(std::cout);
  Matrix composite = get_total_transform_matrix(camera_transform,view_transform,device_transform);
  composite.Print(std::cout);

	for (int vecIndex = 0; vecIndex < triangles.size(); vecIndex++) {
		Triangle t = triangles[vecIndex];
    //print_triangle(t);
    transformTriangle(&t, composite, camera);
    //print_triangle(t);
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
