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
#include <sstream>
#include "VectorOperations.h"
#include "Triangle.h"
#include "MatrixOperations.h"

#define MAX_DEPTH 2

using std::cerr;
using std::endl;
using std::min;
using std::tan;

double ceil441(double f) {
	return ceil(f - 0.00001);
}

double floor441(double f) {
	return floor(f + 0.00001);
}

double cot(double angle) {
	return (1/tan(angle));
}

class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];

    Matrix CameraTransform(void) {
			double x_vector[3],y_vector[3];
    	double z_vector[3] = {position[0] - focus[0],
    												position[1] - focus[1],
    												position[2] - focus[2]};
    	normalize_vector(z_vector);
  		cross_product(up, z_vector, x_vector);
  		normalize_vector(x_vector);
  		cross_product(z_vector, x_vector, y_vector);
  		normalize_vector(y_vector);
    	double frame_vec[3] = {
    		0 - position[0],
    		0 - position[1],
    		0 - position[2]
    	};

    	Matrix m;
			m.A[0][0] = x_vector[0];
			m.A[0][1] = y_vector[0];
			m.A[0][2] = z_vector[0];
			m.A[0][3] = 0;
			m.A[1][0] = x_vector[1];
			m.A[1][1] = y_vector[1];
			m.A[1][2] = z_vector[1];
			m.A[1][3] = 0;
			m.A[2][0] = x_vector[2];
			m.A[2][1] = y_vector[2];
			m.A[2][2] = z_vector[2];
			m.A[2][3] = 0;
			m.A[3][0] = dot_product(x_vector,frame_vec);
			m.A[3][1] = dot_product(y_vector,frame_vec);
			m.A[3][2] = dot_product(z_vector,frame_vec);
			m.A[3][3] = 1;
    	//memcpy(m.A, camera_transform, 16*sizeof(double));
    	return m;
    }

		Matrix ViewTransform(void) {
    	/*double view_transform[4][4] = {
    		cot(angle/2), 0, 0 , 0,
    		0, cot(angle/2), 0 , 0,
    		0, 0 , (far + near)/(far - near), -1,
    		0, 0 , (2*far*near)/(far - near), 0
    	};*/
			Matrix m;
			m.A[0][0] = cot(angle/2);
			m.A[0][1] = 0;
			m.A[0][2] = 0;
			m.A[0][3] = 0;
			m.A[1][0] = 0;
			m.A[1][1] = cot(angle/2);
			m.A[1][2] = 0;
			m.A[1][3] = 0;
			m.A[2][0] = 0;
			m.A[2][1] = 0;
			m.A[2][2] = (far + near)/(far - near);
			m.A[2][3] = -1;
			m.A[3][0] = 0;
			m.A[3][1] = 0;
			m.A[3][2] = (2*far*near)/(far - near);
			m.A[3][3] = 0;
			//memcpy(m.A, camera_transform, 16*sizeof(double));
			return m;
    }
};

Camera
GetCamera(/*double x, double y, double z*/)
{
    Camera c;
    c.near = -80;
    c.far = 80;
    c.angle = 2*M_PI/3;
    c.position[0] = 0;
    c.position[1] = 20;
    c.position[2] = 30;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}

Camera
GetCamera(double x, double y, double z)
{
    Camera c;
    c.near = -80;
    c.far = 80;
    c.angle = 2*M_PI/3;
    c.position[0] = x;
    c.position[1] = y;
    c.position[2] = z;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}

struct LightingParameters
{
    LightingParameters(void)
    {
         light_position[0] = -10;
         light_position[1] = 15;
         light_position[2] = 20;
         Ka = 0.3;
         Kd = 0.05;
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

double calculate_area(double *edge_1, double *edge_2, double *area_vec) {
	cross_product(edge_1,edge_2,area_vec);
	double area = vector_magnitude(area_vec) / 2;
}

void recalculate_barycentric(double *intersect_point, Triangle triangle, double *barycentric) {
	double dummy_intersect[3];
	vector_copy(dummy_intersect, intersect_point);
	if(triangle.X[0] == triangle.X[1] && triangle.X[1] == triangle.X[2] && triangle.X[1] == triangle.X[2])
		dummy_intersect[0] = triangle.X[0];
	if(triangle.Y[0] == triangle.Y[1] && triangle.Y[1] == triangle.Y[2] && triangle.Y[1] == triangle.Y[2])
		dummy_intersect[1] = triangle.Y[0];
	if(triangle.Z[0] == triangle.Z[1] && triangle.Z[1] == triangle.Z[2] && triangle.Z[1] == triangle.Z[2])
		dummy_intersect[2] = triangle.Z[0];
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

		edge_2[0] = dummy_intersect[0] - triangle.X[i];
		edge_2[1] = dummy_intersect[1] - triangle.Y[i];
		edge_2[2] = dummy_intersect[2] - triangle.Z[i];
		double area = calculate_area(edge_1, edge_2, area_vec);
		barycentric[i] = area / t_area;
	}
	barycentric[0] = 1 - (barycentric[1] + barycentric[2]);
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
	}
	recalculate_barycentric(intersect_point, triangle, barycentric);
	return true;
}

double calculate_color(double (*colors)[3], double *barycentric, int component) {
		if(colors[0][component] == colors[1][component]
		&& colors[1][component] == colors[2][component]
		&& colors[2][component] == colors[0][component])
			return colors[0][component];
		else
			return barycentric[1]*colors[0][component]
						+ barycentric[2]*colors[1][component]
						+ barycentric[0]*colors[2][component];
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
		color[0] = min(255.,calculate_color(triangles[object_index].colors,barycentric,0) + triangles[object_index].reflection*color[0]);
		color[1] = min(255.,calculate_color(triangles[object_index].colors,barycentric,1) + triangles[object_index].reflection*color[1]);
		color[2] = min(255.,calculate_color(triangles[object_index].colors,barycentric,2) + triangles[object_index].reflection*color[2]);
	} else {
		if(object_index == -1)
			return;
		color[0] = calculate_color(triangles[object_index].colors,barycentric,0);
		color[1] = calculate_color(triangles[object_index].colors,barycentric,1);
		color[2] = calculate_color(triangles[object_index].colors,barycentric,2);
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
		//if(object_index == 2 || object_index == 3) {
		/*Camera camera = GetCamera();
		double shading_amount = calculate_phong_shading(lp,camera.position,triangles[object_index].normal);*/
		//cout << "shading :" << shading_amount << " triangle " << object_index <<  endl;
		//}
		/*color[0] = min(255., shading_amount*color[0]);
		color[1] = min(255., shading_amount*color[1]);
		color[2] = min(255., shading_amount*color[2]);*/
		/*color[0] = shading_amount*color[0];
		color[1] = shading_amount*color[1];
		color[2] = shading_amount*color[2];*/
	}
};

void get_ray_and_origin(double x, double y, double *ray, double *ray_origin,Screen screen, Matrix world_transform) {
	double aspect_ratio = 1;
	double translated_x = 2*(x+.5)/screen.width - 1;
	double translated_y = 2*(y+.5)/screen.height - 1;
	if(screen.width > screen.height) {
			aspect_ratio = screen.width / screen.height;
			translated_x = translated_x * aspect_ratio;
	} else if (screen.width < screen.height) {
		aspect_ratio = screen.height / screen.width;
		translated_y = translated_y * aspect_ratio;
	}
	double h_focus[4] = {translated_x, translated_y, -1, 1};
	double h_origin[4] = {0,0,0,1};
	double t_h_focus[4], t_h_origin[4];
	world_transform.TransformPoint(h_focus, t_h_focus);
	world_transform.TransformPoint(h_origin, t_h_origin);
	if(t_h_focus[3] != 1) {
		for(int i=0;i < 3;i++)
			t_h_focus[i] = t_h_focus[i] / t_h_focus[3];
	}
	if(t_h_origin[3] != 1) {
		for(int i=0;i < 3;i++)
			t_h_origin[i] = t_h_origin[i] / t_h_origin[3];
	}
	vector_copy(ray_origin, t_h_origin);
	ray[0] = t_h_focus[0] - t_h_origin[0];
	ray[1] = t_h_focus[1] - t_h_origin[1];
	ray[2] = t_h_focus[2] - t_h_origin[2];
	normalize_vector(ray);
}

int main() {
	int height = 1000,width = 1000;
	std::vector<Triangle> triangles = GetTriangles("hardyglobal.vtk");
	cout << "Number of triangles " << triangles.size() << endl;

	double start_x=-30, start_y=0, start_z=30, steps = 500;
	double end_x=-30, end_y=20, end_z=-20;

	/*Camera camera = GetCamera();
	Matrix camera_transform = camera.CameraTransform();
	Matrix world_transform = camera_transform.inverse();*/
	/*Matrix view_transform = camera.ViewTransform();
	Matrix forInverse = Matrix::ComposeMatrices(camera_transform, view_transform);
	Matrix world_transform = forInverse.inverse();*/

	double x_diff = (start_x - end_x) / steps;
	double y_diff = (start_y - end_y) / steps;
	double z_diff = (start_z - end_z) / steps;
	double curr_x = start_x, curr_y = start_y, curr_z = start_z;
	for(int step = 0; step <= 500; step++) {
		vtkImageData *image = NewImage(height, width);
		unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0, 0, 0);
		int npixels = height * width;
		for (int i = 0; i < npixels * 3; i++)
			buffer[i] = 0;

		Screen screen;
		screen.buffer = buffer;
		screen.width = width;
		screen.height = height;

		cout << "Step : " << step << ", x :" << curr_x << ", y :" << curr_y << ", z :" << curr_z << endl;
		Camera camera = GetCamera(curr_x, curr_y, curr_z);
		Matrix camera_transform = camera.CameraTransform();
		Matrix world_transform = camera_transform.inverse();
		for(double x = 0; x < width; x++) {
			for(double y = 0; y < height; y++) {
				double ray[3], ray_origin[3];
				get_ray_and_origin(x, y, ray, ray_origin, screen, world_transform);
				double color[3] = {0,0,0};
				get_color_for_pixel(ray, ray_origin, triangles, color,0,-1);
				screen.find_pixel_and_color(x,y, color);
			}
		}
		std::ostringstream oss;
		/*if((step) / 10 == 0)
			oss << "raytracer00" << step + 500;
		else if (step / 100 == 0)
			oss << "raytracer0" << step + 500;
		else*/
		oss << "raytracer" << step + 500;
		WriteImage(image, oss.str().c_str());
		oss.str("");
		oss.clear();
		free(buffer);
		curr_x = curr_x - x_diff;
		curr_y = curr_y - y_diff;
		curr_z = curr_z - z_diff;
	}
}
