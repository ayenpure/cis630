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
#include <cmath>

using std::cerr;
using std::endl;

double camera_position[3] = {5,5,-20};

double ceil441(double f) {
	return ceil(f - 0.00001);
}

double floor441(double f) {
	return floor(f + 0.00001);
}

double interpolate(double point_1, double point_2, double value_1,
		double value_2, double quest_point) {
	double diff = point_2 - point_1;
	double change_prop =
			(diff != 0) ? (quest_point - point_1) / (point_2 - point_1) : 0;
	double value_at_quest = value_1 + (change_prop * (value_2 - value_1));
	return value_at_quest;
}

double vector_magnitude(double* quest_vec) {
	return sqrt( (quest_vec[0] * quest_vec[0])
						 + (quest_vec[1] * quest_vec[1])
						 + (quest_vec[2] * quest_vec[2]));
}

void normalize_vector(double* quest_normal) {
	double norm = vector_magnitude(quest_normal);
	quest_normal[0] = quest_normal[0] / norm;
	quest_normal[1] = quest_normal[1] / norm;
	quest_normal[2] = quest_normal[2] / norm;
}

double dot_product(double* vector_1, double* vector_2) {
	return (vector_1[0] * vector_2[0]) + (vector_1[1] * vector_2[1])
			+ (vector_1[2] * vector_2[2]);
}

void cross_product(double* vector_1, double* vector_2, double *cross_vec) {
	double product[3] = { (vector_1[1] * vector_2[2])
			- (vector_1[2] * vector_2[1]), (vector_1[2] * vector_2[0])
			- (vector_1[0] * vector_2[2]), (vector_1[0] * vector_2[1])
			- (vector_1[1] * vector_2[0]) };
	memcpy(cross_vec, product, 3 * sizeof(double));
}

class Screen {
public:
	unsigned char *buffer;
	int width, height;

	void find_pixel_and_color(int x, int y, double *color) {
		/*
		 * Ensure the pixels to be painted are in the frame.
		 */
		if ((x >= 0 && x < width) && (y >= 0 && y < height)) {
			int buffer_index = (y * 3 * width) + (x * 3);
			if (buffer_index < width * height * 3) {
				/*double shading_amount = calculate_phong_shading(lp,
						view_direction, current_normal);*/
				buffer[buffer_index++] = color[0];
				buffer[buffer_index++] = color[1];
				buffer[buffer_index] = color[2];
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
	double normal[3];


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

	double getlowestZ() {
		return ((Z[0] < Z[1]) ? Z[0] : ((Z[1] < Z[2]) ? Z[1] : Z[2]));
	}

	double gethighestZ() {
		return ((Z[0] > Z[1]) ? Z[0] : ((Z[1] > Z[2]) ? Z[1] : Z[2]));
	}

	void calculate_normals() {
			if((X[0] == X[1] && X[1] == X[2] && X[2] == X[0])
			&& (Y[0] == Y[1] && Y[1] == Y[2] && Y[2] == Y[0])
			&& (Z[0] == Z[1] && Z[1] == Z[2] && Z[2] == Z[0])){
				normal[0] = X[0];
				normal[1] = Y[0];
				normal[2] = Z[0];
			} else if (X[0] == X[1] && Y[0] == Y[1] && Z[0] == Z[1] ) {
				normal[0] = (X[0] + X[1]) /2;
				normal[1] = (Y[0] + Y[1]) /2;
				normal[2] = (Z[0] + Z[1]) /2;
			} else if (X[1] == X[2] && Y[1] == Y[2] && Z[1] == Z[2] ) {
				normal[0] = (X[1] + X[2]) /2;
				normal[1] = (Y[1] + Y[2]) /2;
				normal[2] = (Z[1] + Z[2]) /2;
			} else if (X[2] == X[0] && Y[2] == Y[0] && Z[2] == Z[0] ) {
				normal[0] = (X[2] + X[0]) /2;
				normal[1] = (Y[2] + Y[0]) /2;
				normal[2] = (Z[2] + Z[0]) /2;
			} else {
				int i = 0,adj_1 = 1,adj_2 = 2;
				double adj_1_vector[3] = { X[adj_1] - X[i], Y[adj_1] - Y[i],
					Z[adj_1] - Z[i] };
				normalize_vector(adj_1_vector);
				double adj_2_vector[3] = { X[adj_2] - X[i], Y[adj_2] - Y[i],
					Z[adj_2] - Z[i] };
				normalize_vector(adj_2_vector);
				cross_product(adj_1_vector, adj_2_vector, normal);
			}
			normalize_vector(normal);
			cout << normal[0] << ", " << normal[1] << ", " << normal[2] << endl;
		cout << endl;
	}
};

std::vector<Triangle> GetTriangles(const char *filename) {

	/*vtkPolyDataReader *rdr = vtkPolyDataReader::New();
	rdr->SetFileName(filename);
	//cerr << "Reading" << endl;
	rdr->Update();
	//cerr << "Done reading" << endl;
	if (rdr->GetOutput()->GetNumberOfCells() == 0) {
		cerr << "Unable to open file!!" << endl;
		exit (EXIT_FAILURE);
	}
	vtkPolyData *pd = rdr->GetOutput();
	int numTris = pd->GetNumberOfCells();
	cout << "Reading " << numTris << " triangles" << endl;
	std::vector<Triangle> tris(numTris);
	vtkPoints *pts = pd->GetPoints();
	vtkCellArray *cells = pd->GetPolys();
	vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
	float *color_ptr = var->GetPointer(0);
	vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
	vtkIdType npts;
	vtkIdType *ptIds;
	int idx;
	for (idx = 0, cells->InitTraversal();
			cells->GetNextCell(npts, ptIds); idx++) {
		if (npts != 3) {
			cerr << "Non-triangles!! ???" << endl;
			exit (EXIT_FAILURE);
		}

		double *pt = NULL;
		pt = pts->GetPoint(ptIds[0]);
		tris[idx].X[0] = pt[0];
		tris[idx].Y[0] = pt[1];
		tris[idx].Z[0] = pt[2];

		pt = pts->GetPoint(ptIds[1]);
		tris[idx].X[1] = pt[0];
		tris[idx].Y[1] = pt[1];
		tris[idx].Z[1] = pt[2];

		pt = pts->GetPoint(ptIds[2]);
		tris[idx].X[2] = pt[0];
		tris[idx].Y[2] = pt[1];
		tris[idx].Z[2] = pt[2];

		double mins[4] = { 1, 2.25, 3.5, 4.75};
		double maxs[4] = { 2.25, 3.5, 4.75, 6};
		unsigned char RGB[5][3] = {
			{0, 0, 255},
			{0, 204, 255},
			{0, 153, 0},
			{255, 204, 0},
			{255, 0, 0},
		 };
		for (int j = 0; j < 3; j++) {
			float val = color_ptr[ptIds[j]];
			int r;
			 for (r = 0 ; r < 4 ; r++)
			 {
			 if (mins[r] <= val && val < maxs[r])
			 break;
			 }
			 if (r == 4)
			 {
			 cerr << "Could not interpolate color for " << val << endl;
			 exit (EXIT_FAILURE);
			 }
			double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
			tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
			tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
			tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
		}
		tris[idx].calculate_normals();
	}
	return tris;*/
	std::vector<Triangle> tris(5);
	Triangle t1;
	t1.X[0] = -5;t1.X[1] = -5;t1.X[2] = 5;
	t1.Y[0] = -5;t1.Y[1] = 5;t1.Y[2] = -5;
	t1.Z[0] = 0;t1.Z[1] = 0;t1.Z[2] = 0;
	t1.colors[0][0] = 0;t1.colors[0][1] = 0;t1.colors[0][2] = 255;
	t1.colors[1][0] = 0;t1.colors[1][1] = 255;t1.colors[1][2] = 0;
	t1.colors[2][0] = 255;t1.colors[2][1] = 0;t1.colors[2][2] = 0;
	t1.calculate_normals();
	tris[0] = t1;
	Triangle t2;
	t2.X[0] = -2.5;t2.X[1] = -2.5;t2.X[2] = 2.5;
	t2.Y[0] = -2.5;t2.Y[1] = 2.5;t2.Y[2] = -2.5;
	t2.Z[0] = 10;t2.Z[1] = 10;t2.Z[2] = 10;
	t2.colors[0][0] = 0;t2.colors[0][1] = 127;t2.colors[0][2] = 255;
	t2.colors[1][0] = 0;t2.colors[1][1] = 127;t2.colors[1][2] = 255;
	t2.colors[2][0] = 0;t2.colors[2][1] = 127;t2.colors[2][2] = 255;
	t2.calculate_normals();
	tris[1] = t2;
	Triangle t3;
	t3.X[0] = -5;t3.X[1] = -5;t3.X[2] = 5;
	t3.Y[0] = -5;t3.Y[1] = 5;t3.Y[2] = -5;
	t3.Z[0] = -5;t3.Z[1] = -5;t3.Z[2] = -5;
	t3.colors[0][0] = 255;t3.colors[0][1] = 127;t3.colors[0][2] = 0;
	t3.colors[1][0] = 255;t3.colors[1][1] = 127;t3.colors[1][2] = 0;
	t3.colors[2][0] = 255;t3.colors[2][1] = 127;t3.colors[2][2] = 0;
	t3.calculate_normals();
	tris[2] = t3;
	Triangle t4;
	t4.X[0] = -20;t4.X[1] = -20;t4.X[2] = 20;
	t4.Y[0] = 10;t4.Y[1] = 10;t4.Y[2] = 10;
	t4.Z[0] = -20;t4.Z[1] = 20;t4.Z[2] = 20;
	t4.colors[0][0] = 255;t4.colors[0][1] = 255;t4.colors[0][2] = 255;
	t4.colors[1][0] = 255;t4.colors[1][1] = 255;t4.colors[1][2] = 255;
	t4.colors[2][0] = 255;t4.colors[2][1] = 255;t4.colors[2][2] = 255;
	t4.calculate_normals();
	tris[3] = t4;
	Triangle t5;
	t5.X[0] = -20;t5.X[1] = 20;t5.X[2] = 20;
	t5.Y[0] = 10;t5.Y[1] = 10;t5.Y[2] = 10;
	t5.Z[0] = -20;t5.Z[1] = -20;t5.Z[2] = 20;
	t5.colors[0][0] = 255;t5.colors[0][1] = 255;t5.colors[0][2] = 255;
	t5.colors[1][0] = 255;t5.colors[1][1] = 255;t5.colors[1][2] = 255;
	t5.colors[2][0] = 255;t5.colors[2][1] = 255;t5.colors[2][2] = 255;
	t5.calculate_normals();
	tris[4] = t5;
	return tris;
}

double calculate_area(double *edge_1, double *edge_2, double *area_vec) {
	cross_product(edge_1,edge_2,area_vec);
	double area = vector_magnitude(area_vec) / 2;
}

bool is_point_inside_triangle(double *intersect_point, Triangle triangle, double *barycentric) {
	/*double temp[3] = {0,0,0};
	for(int i = 0; i < 3; i++) {
		int adj_1 = (i + 1) % 3;
		double vector_1[3] = {
			triangle.X[adj_1] -  triangle.X[i],
			triangle.Y[adj_1] -  triangle.Y[i],
			triangle.Z[adj_1] -  triangle.Z[i]
		};
		double vector_2[3] = {
			intersect_point[0] - triangle.X[i],
			intersect_point[1] - triangle.Y[i],
			intersect_point[2] - triangle.Z[i]
		};
		cross_product(vector_1, vector_2, temp);
		double dot = dot_product(triangle.normal, temp);
		if(dot < 0)
			return false;
	}
	return true;*/
	double edge_1[3] = {
		triangle.X[1] - triangle.X[0],
		triangle.Y[1] - triangle.Y[0],
		triangle.Z[1] - triangle.Z[0]
	};
	double edge_2[3] = {
		triangle.X[2] - triangle.X[0],
		triangle.Y[2] - triangle.Y[0],
		triangle.Z[2] - triangle.Z[0]
	};
	double area_vec[3];
	double t_area = calculate_area(edge_1, edge_2, area_vec);
	//calculate barycentric coordinates
	for(int i = 0; i < 3; i++) {
		int adj = (i + 1) % 3;
		edge_1[0] = triangle.X[adj] - triangle.X[i];
		edge_1[1] = triangle.Y[adj] - triangle.Y[i];
		edge_1[2] = triangle.Z[adj] - triangle.Z[i];

		edge_2[0] = intersect_point[0] - triangle.X[i];
		edge_2[1] = intersect_point[1] - triangle.Y[i];
		edge_2[2] = intersect_point[2] - triangle.Z[i];
		double area = calculate_area(edge_1, edge_2, area_vec);
		if(dot_product(triangle.normal, area_vec) < 0)
			return false;
		barycentric[i] = area / t_area;
	}
	return true;
}

double calculate_shading(double (*colors)[3], double *barycentric, int component) {
		if(colors[0][component] == colors[1][component]
		&& colors[1][component] == colors[2][component]
		&& colors[2][component] == colors[0][component])
			return colors[0][component];
		else
			return barycentric[0]*colors[0][component]
						+ barycentric[1]*colors[1][component]
						+ barycentric[2]*colors[2][component];
}

void get_color_for_pixel(double *ray, std::vector<Triangle> triangles, double * color) {
	double distance = std::numeric_limits<double>::max();
	for(int i = 0; i < triangles.size(); i++) {
		//Check if triangle is a single pixel triangle
		//if not get point of intersection of the triangle with the triangle plane
		if(dot_product(ray, triangles[i].normal) != 0) {
			double triangle_vertex[3] = {
				triangles[i].X[0],
				triangles[i].Y[0],
				triangles[i].Z[0]
			};
			double o_distance = dot_product(triangles[i].normal,triangle_vertex);
			//float t = (dot(N, orig) + D) / dot(N, dir);
			double c_distance = - ((dot_product(triangles[i].normal, camera_position) + o_distance) / dot_product(triangles[i].normal, ray));
			//Vec3f Phit = orig + t * dirtrina
			double intersect_point[3] = {
				camera_position[0] + c_distance*ray[0],
				camera_position[1] + c_distance*ray[1],
				camera_position[2] + c_distance*ray[2]
			};
			double barycentric[3] = {0,0,0};
			if(is_point_inside_triangle(intersect_point, triangles[i], barycentric)) {
				if(c_distance < distance){
					distance = c_distance;
					color[0] = calculate_shading(triangles[i].colors, barycentric, 0);
					color[1] = calculate_shading(triangles[i].colors, barycentric, 1);
					color[2] = calculate_shading(triangles[i].colors, barycentric, 2);
				}
			}
		}
	}
};

int main() {
	int height = 1000,width = 1000;
	std::vector<Triangle> triangles = GetTriangles("hardyglobal.0.vtk");
	cout << "Number of triangles " << triangles.size() << endl;
	vtkImageData *image = NewImage(height, width);
	unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0, 0, 0);
	int npixels = height * width;
	for (int i = 0; i < npixels * 3; i++)
		buffer[i] = 0;
	Screen screen;
	screen.buffer = buffer;
	screen.width = height;
	screen.height = width;
	for(int x = 0; x < width; x++) {
		for(int y = 0; y < height; y++) {
			double translated_x = (20*(double)x)/(double)width - 10;
			double translated_y = (20*(double)y)/(double)height - 10;
			double look_at[3] = {translated_x, translated_y, -10};
			double ray[3] = {
				look_at[0] - camera_position[0],
				look_at[1] - camera_position[1],
				look_at[2] - camera_position[2]
			};
			normalize_vector(ray);
			double color[3] = {0,0,0};
			get_color_for_pixel(ray, triangles, color);
			screen.find_pixel_and_color(x,y, color);
			if(color[0] != 0 || color[1] != 0 || color[2] != 0)
				cout << "coloring pixel " << x << ", " << y << endl;
		}
	}
	WriteImage(image, "raytracer");
	free(buffer);
}
