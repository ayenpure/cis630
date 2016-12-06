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

#define MAX_DEPTH 2

using std::cerr;
using std::endl;
using std::min;

double camera_position[3] = {5,5,-20};
//double light_position[3] = {-10,15,-20};

double ceil441(double f) {
	return ceil(f - 0.00001);
}

double floor441(double f) {
	return floor(f + 0.00001);
}

struct LightingParameters
{
    LightingParameters(void)
    {
         light_position[0] = -10;
         light_position[1] = 15;
         light_position[2] = -20;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 5.3;
         alpha = 7.5;
    };
    double light_position[3];  // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};

LightingParameters lp;

void vector_copy(double *destination, double *source) {
	destination[0] = source[0];
	destination[1] = source[1];
	destination[2] = source[2];
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
	vector_copy(cross_vec, product);
}

double calculate_for_specular_lighting(LightingParameters lp, double *V,
		double* N) {
	/* R = 2*(L . N)*N - L
	 * return ==> V . R
	 */
	double two_L_dot_N = 2 * dot_product(lp.light_position, N);
	double two_L_dot_N_N[3] = { two_L_dot_N * N[0], two_L_dot_N * N[1],
			two_L_dot_N * N[2] };
	double R[3] = { two_L_dot_N_N[0] - lp.light_position[0], two_L_dot_N_N[1] - lp.light_position[1],
			two_L_dot_N_N[2] - lp.light_position[2] };
	return dot_product(R, V);
}

double calculate_phong_shading(LightingParameters lp, double *view_direction,
		double *normal) {
	double diffuse_component = dot_product(lp.light_position, normal);
  double specular_component = pow(calculate_for_specular_lighting(lp, view_direction, normal),lp.alpha);
  if(specular_component < 0 || specular_component != specular_component) {
      specular_component = 0;
  }
	double shading_amount = lp.Ka + lp.Kd * abs(diffuse_component)
			+ lp.Ks * specular_component;
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
	double reflection;

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
	t1.reflection = 0;
	t1.X[0] = -5;t1.X[1] = -5;t1.X[2] = 5;
	t1.Y[0] = -10;t1.Y[1] = 0;t1.Y[2] = -10;
	t1.Z[0] = 0;t1.Z[1] = 0;t1.Z[2] = 0;
	t1.colors[0][0] = 255;t1.colors[0][1] = 0;t1.colors[0][2] = 0;
	t1.colors[1][0] = 0;t1.colors[1][1] = 255;t1.colors[1][2] = 0;
	t1.colors[2][0] = 0;t1.colors[2][1] = 0;t1.colors[2][2] = 255;
	t1.calculate_normals();
	tris[0] = t1;
	Triangle t2;
	t2.reflection = 0;
	t2.X[0] = -2.5;t2.X[1] = -2.5;t2.X[2] = 2.5;
	t2.Y[0] = -2.5;t2.Y[1] = 2.5;t2.Y[2] = -2.5;
	t2.Z[0] = 10;t2.Z[1] = 10;t2.Z[2] = 10;
	t2.colors[0][0] = 0;t2.colors[0][1] = 127;t2.colors[0][2] = 255;
	t2.colors[1][0] = 0;t2.colors[1][1] = 127;t2.colors[1][2] = 255;
	t2.colors[2][0] = 0;t2.colors[2][1] = 127;t2.colors[2][2] = 255;
	t2.calculate_normals();
	tris[1] = t2;
	Triangle t3;
	t3.reflection = 0;
	t3.X[0] = -5;t3.X[1] = -5;t3.X[2] = 5;
	t3.Y[0] = -5;t3.Y[1] = 5;t3.Y[2] = -5;
	t3.Z[0] = -5;t3.Z[1] = -5;t3.Z[2] = -5;
	t3.colors[0][0] = 255;t3.colors[0][1] = 127;t3.colors[0][2] = 0;
	t3.colors[1][0] = 255;t3.colors[1][1] = 127;t3.colors[1][2] = 0;
	t3.colors[2][0] = 255;t3.colors[2][1] = 127;t3.colors[2][2] = 0;
	t3.calculate_normals();
	tris[2] = t3;
	Triangle t4;
	t4.reflection = 0.5;
	t4.X[0] = -20;t4.X[1] = -20;t4.X[2] = 20;
	t4.Y[0] = 10;t4.Y[1] = 10;t4.Y[2] = 10;
	t4.Z[0] = 20;t4.Z[1] = -20;t4.Z[2] = 20;
	t4.colors[0][0] = 127;t4.colors[0][1] = 127;t4.colors[0][2] = 127;
	t4.colors[1][0] = 127;t4.colors[1][1] = 127;t4.colors[1][2] = 127;
	t4.colors[2][0] = 127;t4.colors[2][1] = 127;t4.colors[2][2] = 127;
	t4.calculate_normals();
	tris[3] = t4;
	Triangle t5;
	t5.reflection = 0.5;
	t5.X[0] = 20;t5.X[1] = -20;t5.X[2] = 20;
	t5.Y[0] = 10;t5.Y[1] = 10;t5.Y[2] = 10;
	t5.Z[0] = 20;t5.Z[1] = -20;t5.Z[2] = -20;
	t5.colors[0][0] = 127;t5.colors[0][1] = 127;t5.colors[0][2] = 127;
	t5.colors[1][0] = 127;t5.colors[1][1] = 127;t5.colors[1][2] = 127;
	t5.colors[2][0] = 127;t5.colors[2][1] = 127;t5.colors[2][2] = 127;
	t5.calculate_normals();
	tris[4] = t5;
	return tris;
}

double calculate_area(double *edge_1, double *edge_2, double *area_vec) {
	cross_product(edge_1,edge_2,area_vec);
	double area = vector_magnitude(area_vec) / 2;
}

bool is_point_inside_triangle(double *intersect_point, Triangle triangle, double *barycentric) {
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

void get_color_for_pixel(double *ray, double * ray_origin, std::vector<Triangle> triangles, double * color, int depth,int skip_index) {
	double distance = std::numeric_limits<double>::max();
	int object_index = -1;
	double barycentric[3];
	double intersect_point[3];
	for(int i = 0; i < triangles.size(); i++) {
		if(i == skip_index)
			continue;
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
			double c_distance = - ((dot_product(triangles[i].normal, ray_origin) + o_distance) / dot_product(triangles[i].normal, ray));
			//Vec3f Phit = orig + t * dirtrina
			if(c_distance < 0)
				continue;
			double c_intersect_point[3] = {
				ray_origin[0] + c_distance*ray[0],
				ray_origin[1] + c_distance*ray[1],
				ray_origin[2] + c_distance*ray[2]
			};
			double c_barycentric[3] = {0,0,0};
			if(is_point_inside_triangle(c_intersect_point, triangles[i], c_barycentric)) {
				if(c_distance < distance){
					distance = c_distance;
					object_index = i;
					vector_copy(barycentric, c_barycentric);
					vector_copy(intersect_point, c_intersect_point);
				}
			}
		}
	}
	if(object_index != -1 && triangles[object_index].reflection > 0 && depth < MAX_DEPTH) {
		//Reflection
		double angle = dot_product(triangles[object_index].normal, ray);
		double reflection_ray[3] = {
			ray[0] - 2*angle*triangles[object_index].normal[0],
			ray[1] - 2*angle*triangles[object_index].normal[1],
			ray[2] - 2*angle*triangles[object_index].normal[2],
		};
		normalize_vector(reflection_ray);
		get_color_for_pixel(reflection_ray, intersect_point, triangles, color, ++depth, object_index);
		--depth;
		color[0] = min(255.,calculate_shading(triangles[object_index].colors,barycentric,0) + triangles[object_index].reflection*color[0]);
		color[1] = min(255.,calculate_shading(triangles[object_index].colors,barycentric,1) + triangles[object_index].reflection*color[1]);
		color[2] = min(255.,calculate_shading(triangles[object_index].colors,barycentric,2) + triangles[object_index].reflection*color[2]);
	} else {
		if(object_index == -1)
			return;
		color[0] = calculate_shading(triangles[object_index].colors,barycentric,0);
		color[1] = calculate_shading(triangles[object_index].colors,barycentric,1);
		color[2] = calculate_shading(triangles[object_index].colors,barycentric,2);
	}
	if(depth == 0) {
		double dummy_color[3] = {0,0,0};
		double shadow_ray[3] =  {
			lp.light_position[0] - intersect_point[0],
			lp.light_position[1] - intersect_point[1],
			lp.light_position[2] - intersect_point[2]
		};
		normalize_vector(shadow_ray);
		get_color_for_pixel(shadow_ray, intersect_point, triangles, dummy_color, MAX_DEPTH, object_index);
		if(dummy_color[0] != 0 || dummy_color[1] != 0 || dummy_color[2] != 0) {
			color[0] /=2;
			color[1] /=2;
			color[2] /=2;
		}
		double shading_amount = calculate_phong_shading(lp,camera_position,triangles[object_index].normal);
		/*color[0] = min(255., shading_amount*color[0]);
		color[1] = min(255., shading_amount*color[1]);
		color[2] = min(255., shading_amount*color[2]);*/
		/*color[0] = shading_amount*color[0];
		color[1] = shading_amount*color[1];
		color[2] = shading_amount*color[2];*/
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
			double translated_x = (20*(double)(x+.5))/(double)width - 10;
			double translated_y = (20*(double)(y+.5))/(double)height - 10;
			double look_at[3] = {translated_x, translated_y, -10};
			double ray[3] = {
				look_at[0] - camera_position[0],
				look_at[1] - camera_position[1],
				look_at[2] - camera_position[2]
			};
			normalize_vector(ray);
			double color[3] = {0,0,0};
			get_color_for_pixel(ray, camera_position, triangles, color,0,-1);
			cout << "Coloring pixel " << x << ", " << y << endl;
			screen.find_pixel_and_color(x,y, color);
		}
	}
	WriteImage(image, "raytracer");
	free(buffer);
}
