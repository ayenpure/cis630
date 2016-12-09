#ifndef CAMERAPOSITIONS_H_
#define CAMERAPOSITIONS_H_

#include <iostream>
#include "VectorOperations.h"

using std::cout;
using std::endl;

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
			tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]));///255.0;
			tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]));//255.0;
			tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]));//255.0;
		}
		tris[idx].calculate_normals();
	}*/
	//return tris;
	std::vector<Triangle> tris(20);
	Triangle t1;
	t1.reflection = 0.5;
	t1.X[0] = -20;t1.X[1] = -20;t1.X[2] = 20;
	t1.Y[0] = 10;t1.Y[1] = 10;t1.Y[2] = 10;
	t1.Z[0] = 20;t1.Z[1] = -20;t1.Z[2] = 20;
	t1.colors[0][0] = 100;t1.colors[0][1] = 100;t1.colors[0][2] = 100;
	t1.colors[1][0] = 100;t1.colors[1][1] = 100;t1.colors[1][2] = 100;
	t1.colors[2][0] = 100;t1.colors[2][1] = 100;t1.colors[2][2] = 100;
	t1.calculate_normals();
	tris[0] = t1;

	Triangle t2;
	t2.reflection = 0.5;
	t2.X[0] = 20;t2.X[1] = -20;t2.X[2] = 20;
	t2.Y[0] = 10;t2.Y[1] = 10;t2.Y[2] = 10;
	t2.Z[0] = 20;t2.Z[1] = -20;t2.Z[2] = -20;
	t2.colors[0][0] = 100;t2.colors[0][1] = 100;t2.colors[0][2] = 100;
	t2.colors[1][0] = 100;t2.colors[1][1] = 100;t2.colors[1][2] = 100;
	t2.colors[2][0] = 100;t2.colors[2][1] = 100;t2.colors[2][2] = 100;
	t2.calculate_normals();
	tris[1] = t2;

	Triangle t3;
	t3.reflection = 0;
	t3.X[0] = -15;t3.X[1] = -15; t3.X[2] = -5;
	t3.Y[0] = 10;t3.Y[1] = -10; t3.Y[2] = -10;
	t3.Z[0] = -5;t3.Z[1] = -5; t3.Z[2] = -5;
	t3.colors[0][0] = 255;t3.colors[0][1] = 0;t3.colors[0][2] = 0;
	t3.colors[1][0] = 0;t3.colors[1][1] = 255;t3.colors[1][2] = 0;
	t3.colors[2][0] = 0;t3.colors[2][1] = 0;t3.colors[2][2] = 255;
	t3.calculate_normals();
	tris[2] = t3;

	Triangle t4;
	t4.reflection = 0;
	t4.X[0] = -15;t4.X[1] = -5;t4.X[2] = -5;
	t4.Y[0] = 10;t4.Y[1] = -10;t4.Y[2] = 10;
	t4.Z[0] = -5;t4.Z[1] = -5;t4.Z[2] = -5;
	t4.colors[0][0] = 255;t4.colors[0][1] = 0;t4.colors[0][2] = 0;
	t4.colors[1][0] = 0;t4.colors[1][1] = 0;t4.colors[1][2] = 255;
	t4.colors[2][0] = 0;t4.colors[2][1] = 255;t4.colors[2][2] = 0;
	t4.calculate_normals();
	tris[3] = t4;

	Triangle t5;
	t5.reflection = 0;
	t5.X[0] = 5;t5.X[1] =  5;t5.X[2] = 5;
	t5.Y[0] = 10;t5.Y[1] = -10;t5.Y[2] = -10;
	t5.Z[0] = 5;t5.Z[1] = 5;t5.Z[2] = -5;
	t5.colors[0][0] = 0;t5.colors[0][1] = 255;t5.colors[0][2] = 0;
	t5.colors[1][0] = 0;t5.colors[1][1] = 0;t5.colors[1][2] = 255;
	t5.colors[2][0] = 255;t5.colors[2][1] = 0;t5.colors[2][2] = 0;
	t5.calculate_normals();
	tris[4] = t5;

	Triangle t6;
	t6.reflection = 0;
	t6.X[0] = 5;t6.X[1] = 5;t6.X[2] = 5;
	t6.Y[0] = 10;t6.Y[1] = -10;t6.Y[2] = 10;
	t6.Z[0] = 5;t6.Z[1] = -5;t6.Z[2] = -5;
	t6.colors[0][0] = 0;t6.colors[0][1] = 255;t6.colors[0][2] = 0;
	t6.colors[1][0] = 255;t6.colors[1][1] = 0;t6.colors[1][2] = 0;
	t6.colors[2][0] = 0;t6.colors[2][1] = 0;t6.colors[2][2] = 255;
	t6.calculate_normals();
	tris[5] = t6;

	Triangle t7;
	t7.reflection = 0;
	t7.X[0] = -15;t7.X[1] = -15; t7.X[2] = -5;
	t7.Y[0] = 10;t7.Y[1] = -10; t7.Y[2] = -10;
	t7.Z[0] = 5;t7.Z[1] = 5; t7.Z[2] = 5;
	t7.colors[0][0] = 255;t7.colors[0][1] = 0;t7.colors[0][2] = 0;
	t7.colors[1][0] = 0;t7.colors[1][1] = 255;t7.colors[1][2] = 0;
	t7.colors[2][0] = 0;t7.colors[2][1] = 0;t7.colors[2][2] = 255;
	t7.calculate_normals();
	tris[6] = t7;

	Triangle t8;
	t8.reflection = 0;
	t8.X[0] = -15;t8.X[1] = -5;t8.X[2] = -5;
	t8.Y[0] = 10;t8.Y[1] = -10;t8.Y[2] = 10;
	t8.Z[0] = 5;t8.Z[1] = 5;t8.Z[2] = 5;
	t8.colors[0][0] = 255;t8.colors[0][1] = 0;t8.colors[0][2] = 0;
	t8.colors[1][0] = 0;t8.colors[1][1] = 0;t8.colors[1][2] = 255;
	t8.colors[2][0] = 0;t8.colors[2][1] = 255;t8.colors[2][2] = 0;
	t8.calculate_normals();
	tris[7] = t8;

	Triangle t9;
	t9.reflection = 0;
	t9.X[0] = 15;t9.X[1] =  15;t9.X[2] = 15;
	t9.Y[0] = 10;t9.Y[1] = -10;t9.Y[2] = -10;
	t9.Z[0] = 5;t9.Z[1] = 5;t9.Z[2] = -5;
	t9.colors[0][0] = 0;t9.colors[0][1] = 255;t9.colors[0][2] = 0;
	t9.colors[1][0] = 0;t9.colors[1][1] = 0;t9.colors[1][2] = 255;
	t9.colors[2][0] = 255;t9.colors[2][1] = 0;t9.colors[2][2] = 0;
	t9.calculate_normals();
	tris[8] = t9;

	Triangle t10;
	t10.reflection = 0;
	t10.X[0] = 15;t10.X[1] = 15;t10.X[2] = 15;
	t10.Y[0] = 10;t10.Y[1] = -10;t10.Y[2] = 10;
	t10.Z[0] = 5;t10.Z[1] = -5;t10.Z[2] = -5;
	t10.colors[0][0] = 0;t10.colors[0][1] = 255;t10.colors[0][2] = 0;
	t10.colors[1][0] = 255;t10.colors[1][1] = 0;t10.colors[1][2] = 0;
	t10.colors[2][0] = 0;t10.colors[2][1] = 0;t10.colors[2][2] = 255;
	t10.calculate_normals();
	tris[9] = t10;

	Triangle t11;
	t11.reflection = 0;
	t11.X[0] = -15;t11.X[1] = -5;t11.X[2] = -5;
	t11.Y[0] = -10;t11.Y[1] = -10;t11.Y[2] = -10;
	t11.Z[0] = 5;t11.Z[1] = 5;t11.Z[2] = -5;
	t11.colors[0][0] = 255;t11.colors[0][1] = 0;t11.colors[0][2] = 0;
	t11.colors[1][0] = 0;t11.colors[1][1] = 255;t11.colors[1][2] = 0;
	t11.colors[2][0] = 0;t11.colors[2][1] = 0;t11.colors[2][2] = 255;
	t11.calculate_normals();
	tris[10] = t11;

	Triangle t12;
	t12.reflection = 0;
	t12.X[0] = -15;t12.X[1] = -15;t12.X[2] = -5;
	t12.Y[0] = -10;t12.Y[1] = -10;t12.Y[2] = -10;
	t12.Z[0] = -5;t12.Z[1] = 5;t12.Z[2] = -5;
	t12.colors[0][0] = 0;t12.colors[0][1] = 255;t12.colors[0][2] = 0;
	t12.colors[1][0] = 0;t12.colors[1][1] = 0;t12.colors[1][2] = 255;
	t12.colors[2][0] = 255;t12.colors[2][1] = 0;t12.colors[2][2] = 0;
	t12.calculate_normals();
	tris[11] = t12;

	Triangle t13;
	t13.reflection = 0;
	t13.X[0] = 10;t13.X[1] = 15; t13.X[2] = 5;
	t13.Y[0] = -10;t13.Y[1] = 0; t13.Y[2] = 0;
	t13.Z[0] = 5;t13.Z[1] = 0; t13.Z[2] = 0;
	t13.colors[0][0] = 255;t13.colors[0][1] = 255;t13.colors[0][2] = 0;
	t13.colors[1][0] = 0;t13.colors[1][1] = 255;t13.colors[1][2] = 255;
	t13.colors[2][0] = 255;t13.colors[2][1] = 0;t13.colors[2][2] = 255;
	t13.calculate_normals();
	tris[12] = t13;

	Triangle t14;
	t14.reflection = 0;
	t14.X[0] = 10;t14.X[1] = 0;t14.X[2] = 5;
	t14.Y[0] = -10;t14.Y[1] = -10;t14.Y[2] = 0;
	t14.Z[0] = 5;t14.Z[1] = 5;t14.Z[2] = 0;
	t14.colors[0][0] = 255;t14.colors[0][1] = 255;t14.colors[0][2] = 0;
	t14.colors[1][0] = 0;t14.colors[1][1] = 255;t14.colors[1][2] = 255;
	t14.colors[2][0] = 255;t14.colors[2][1] = 0;t14.colors[2][2] = 255;
	t14.calculate_normals();
	tris[13] = t14;

	Triangle t15;
	t15.reflection = 0;
	t15.X[0] = 10;t15.X[1] =  20;t15.X[2] = 15;
	t15.Y[0] = -10;t15.Y[1] = -10;t15.Y[2] = 0;
	t15.Z[0] = -5;t15.Z[1] = -5;t15.Z[2] = 0;
	t15.colors[0][0] = 0;t15.colors[0][1] = 255;t15.colors[0][2] = 255;
	t15.colors[1][0] = 255;t15.colors[1][1] = 255;t15.colors[1][2] = 0;
	t15.colors[2][0] = 255;t15.colors[2][1] = 0;t15.colors[2][2] = 255;
	t15.calculate_normals();
	tris[14] = t15;

	Triangle t16;
	t16.reflection = 0;
	t16.X[0] = 5;t16.X[1] =  10;t16.X[2] = 15;
	t16.Y[0] = 0;t16.Y[1] = -10;t16.Y[2] = 0;
	t16.Z[0] = 0;t16.Z[1] = -5;t16.Z[2] = 0;
	t16.colors[0][0] = 255;t16.colors[0][1] = 255;t16.colors[0][2] = 0;
	t16.colors[1][0] = 0;t16.colors[1][1] = 255;t16.colors[1][2] = 255;
	t16.colors[2][0] = 255;t16.colors[2][1] = 0;t16.colors[2][2] = 255;
	t16.calculate_normals();
	tris[15] = t16;

  Triangle t17;
	t17.reflection = 0.5;
	t17.X[0] = -20;t17.X[1] = -20;t17.X[2] = -20;
	t17.Y[0] = -10;t17.Y[1] = -10;t17.Y[2] = 20;
	t17.Z[0] = -20;t17.Z[1] = 20;t17.Z[2] = 20;
	t17.colors[0][0] = 127;t17.colors[0][1] = 127;t17.colors[0][2] = 127;
	t17.colors[1][0] = 127;t17.colors[1][1] = 127;t17.colors[1][2] = 127;
	t17.colors[2][0] = 127;t17.colors[2][1] = 127;t17.colors[2][2] = 127;
	t17.calculate_normals();
	tris[16] = t17;

	Triangle t18;
	t18.reflection = 0.5;
	t18.X[0] = -20;t18.X[1] = -20;t18.X[2] = -20;
	t18.Y[0] = 20;t18.Y[1] = -10;t18.Y[2] = 20;
	t18.Z[0] = 20;t18.Z[1] = -20;t18.Z[2] = -20;
	t18.colors[0][0] = 127;t18.colors[0][1] = 127;t18.colors[0][2] = 127;
	t18.colors[1][0] = 127;t18.colors[1][1] = 127;t18.colors[1][2] = 127;
	t18.colors[2][0] = 127;t18.colors[2][1] = 127;t18.colors[2][2] = 127;
	t18.calculate_normals();
	tris[17] = t18;

  Triangle t19;
	t19.reflection = 0;
	t19.X[0] = -20;t19.X[1] = -20;t19.X[2] = 20;
	t19.Y[0] = 20;t19.Y[1] = -10;t19.Y[2] = -10;
	t19.Z[0] = 20;t19.Z[1] = 20;t19.Z[2] = 20;
	t19.colors[0][0] = 127;t19.colors[0][1] = 127;t19.colors[0][2] = 127;
	t19.colors[1][0] = 127;t19.colors[1][1] = 127;t19.colors[1][2] = 127;
	t19.colors[2][0] = 127;t19.colors[2][1] = 127;t19.colors[2][2] = 127;
	t19.calculate_normals();
	tris[18] = t19;

  Triangle t20;
	t20.reflection = 0;
	t20.X[0] = -20;t20.X[1] = 20;t20.X[2] = 20;
	t20.Y[0] = 20;t20.Y[1] = -10;t20.Y[2] = 20;
	t20.Z[0] = 20;t20.Z[1] = 20;t20.Z[2] = 20;
	t20.colors[0][0] = 127;t20.colors[0][1] = 127;t20.colors[0][2] = 127;
	t20.colors[1][0] = 127;t20.colors[1][1] = 127;t20.colors[1][2] = 127;
	t20.colors[2][0] = 127;t20.colors[2][1] = 127;t20.colors[2][2] = 127;
	t20.calculate_normals();
	tris[19] = t20;

	return tris;
}

#endif
