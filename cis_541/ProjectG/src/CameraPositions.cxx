/*=========================================================================

 Program:   Visualization Toolkit
 Module:    SpecularSpheres.cxx

 Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 All rights reserved.
 See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

 =========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <iostream>
#include <cmath>
#include <string.h>
#include <stdio.h>

#define X 'x'
#define Y 'y'
#define Z 'z'

void rotate(double* camera_position, double angle, char axis,double* rotated_camera) {
	double rotation_matrix[3][3];
	if(axis == X) {
		double x_rotation_matrix[3][3] = {
			1,0,0,
			0,cos(angle),-sin(angle),
			0,sin(angle),cos(angle)
		};
		memcpy(rotation_matrix, x_rotation_matrix, 9*sizeof(double));
	} else if (axis == Y) {
		double y_rotation_matrix[3][3] = {
			cos(angle),0,-sin(angle),
			0,1,0,
			sin(angle),0,cos(angle)
		};
		memcpy(rotation_matrix, y_rotation_matrix, 9*sizeof(double));
	} else if (axis == Z) {
		double z_rotation_matrix[3][3] = {
			cos(angle),-sin(angle),0,
			sin(angle),cos(angle),0,
			0,0,1
		};
		memcpy(rotation_matrix, z_rotation_matrix, 9*sizeof(double));
	}
	/*for (int i = 0 ; i < 4 ; i++)
	{
			char str[256];
			sprintf(str, "(%.7f %.7f %.7f)\n", rotation_matrix[i][0], rotation_matrix[i][1], rotation_matrix[i][2]);
			cout << str;
	}*/
	for(int i = 0; i < 3; i++) {
		rotated_camera[i] = 0;
		for(int j = 0; j < 3; j++) {
			rotated_camera[i] = rotated_camera[i] + camera_position[j]*rotation_matrix[j][i];
		}
		if(abs(rotated_camera[i]) < 0.0000001)
			rotated_camera[i] = 0;
	}
	//cout << "rotated position {" << rotated_camera[0] << ", " << rotated_camera[1] << ", " << rotated_camera[2] << " }" << endl;
}

void get_camera_positions(double (*camera_positions)[3]) {
	double calc_camera_positions[8][16][3];
	double first_rotation[16] = {
		0, (M_PI/6), (M_PI/4), (M_PI/3),
		(M_PI/2), (2*M_PI/3), (3*M_PI/4), (5*M_PI/6),
		(M_PI), (7*M_PI/6), (5*M_PI/4), (4*M_PI/3),
		(3*M_PI/2), (5*M_PI/3),(7*M_PI/4) ,(11*M_PI/6)
	};
	double second_rotation[7] = {
		(M_PI/6), (M_PI/4), (M_PI/3), (M_PI/2), (2*M_PI/3), (3*M_PI/4), (5*M_PI/6)
	};
	double camera_position[3] = {40,0,0};
	for(int i = 0; i < 16; i++) {
		rotate(camera_position, first_rotation[i], Y, calc_camera_positions[0][i]);
	};
	for(int i = 1; i < 8; i++) {
		for(int j = 0; j < 16; j++) {
			rotate(calc_camera_positions[0][j], second_rotation[i-1], X, calc_camera_positions[i][j]);
		}
	};
	int positions_count = 0;
	for(int i = 0; i < 8 ; i++) {
		for (int j = 0; j <16; j++) {
			if(i > 0 && (j==0 || j==8))
				continue;
			camera_positions[positions_count][0] = calc_camera_positions[i][j][0];
			camera_positions[positions_count][1] = calc_camera_positions[i][j][1];
			camera_positions[positions_count][2] = calc_camera_positions[i][j][2];
			positions_count++;
		}
	}
}


class vtk441Mapper: public vtkOpenGLPolyDataMapper {
protected:
	GLuint displayList;
	bool initialized;

public:
	vtk441Mapper() {
		initialized = false;
	}

	void RemoveVTKOpenGLStateSideEffects() {
		float Info[4] = { 0, 0, 0, 1 };
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Info);
		float ambient[4] = { 1, 1, 1, 1.0 };
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
		float diffuse[4] = { 1, 1, 1, 1.0 };
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
		float specular[4] = { 1, 1, 1, 1.0 };
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
		glDisable (GL_TEXTURE_1D);
		glDisable (GL_COLOR_MATERIAL);
	}

	void SetupLight(void) {
		glEnable (GL_LIGHTING);
		glEnable (GL_LIGHT0);
		GLfloat diffuse0[4] = { 0.8, 0.8, 0.8, 1 };
		GLfloat ambient0[4] = { 0.2, 0.2, 0.2, 1 };
		GLfloat specular0[4] = { 0.0, 0.0, 0.0, 1 };
		GLfloat pos0[4] = { 1, 2, 3, 0 };
		glLightfv(GL_LIGHT0, GL_POSITION, pos0);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse0);
		glLightfv(GL_LIGHT0, GL_AMBIENT, ambient0);
		glLightfv(GL_LIGHT0, GL_SPECULAR, specular0);
		glDisable (GL_LIGHT1);
		glDisable (GL_LIGHT2);
		glDisable (GL_LIGHT3);
		glDisable (GL_LIGHT5);
		glDisable (GL_LIGHT6);
		glDisable (GL_LIGHT7);
	}
};

GLfloat ctrlpoints[4][4][3] = { { { -1.5, -1.5, 1.0 }, { -0.5, -1.5, 1.0 }, {
		0.5, -1.5, 1.0 }, { 1.5, -1.5, 1.0 } }, { { -1.5, -0.5, 1.0 }, { -0.5,
		-0.5, 1.0 }, { 0.5, -0.5, 1.0 }, { 1.5, -0.5, 1.0 } }, { { -1.5, 0.5,
		1.0 }, { -0.5, 0.5, 1.0 }, { 0.5, 0.5, 1.0 }, { 1.5, 0.5, 1.0 } }, { {
		-1.5, 1.5, 1.0 }, { -0.5, 1.5, 1.0 }, { 0.5, 1.5, 1.0 },
		{ 1.5, 1.5, 1.0 } } };

class vtk441MapperPart1: public vtk441Mapper {
public:
	static vtk441MapperPart1 *New();

	virtual void RenderPiece(vtkRenderer *ren, vtkActor *act) {
		RemoveVTKOpenGLStateSideEffects();
		//glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );
		glEnable ( GL_COLOR_MATERIAL );
		SetupLight();
		double camera_positions[114][3];
		get_camera_positions(camera_positions);
		for(int cam_index = 0; cam_index < 114; cam_index++) {
			double focus[3];
			int focus_index = 0;
			if (0) {
				focus[0] = 0;
				focus[1] = 0;
				focus[2] = 0;
			} else if(1){
				//double focus[3] = {0,0,0};
				//camera = GetCamera(camera_positions[cam_index], focus);
				if (cam_index < 16)
					focus_index = (cam_index + 7) % 16;
				else if (0) {
					int quotient = (cam_index - 16) / 14;
					int pseudo_index = (cam_index - 16) % 14;
					int pseudo_focus_index = (pseudo_index + 6) % 14;
					focus_index = ((quotient * 14) + 16) + pseudo_focus_index;
				} else {
					int pseudo_index = cam_index -16 + 70;
					if (pseudo_index >= 98) {
						focus_index = (pseudo_index % 98) + 7;
					} else {
						focus_index = pseudo_index;
					}
					focus_index+=16;
				}
				focus[0] = camera_positions[focus_index][0];
				focus[1] = camera_positions[focus_index][1];
				focus[2] = camera_positions[focus_index][2];
			}
			glBegin(GL_LINES);
			glColor3f(1,1,1);
			glVertex3f(camera_positions[cam_index][0],camera_positions[cam_index][1],camera_positions[cam_index][2]);
			glColor3f(1,0,0);
			glVertex3f(focus[0],focus[1],focus[2]);
			glEnd();
		}
	}
};

vtkStandardNewMacro (vtk441MapperPart1);

int main() {
	// Dummy input so VTK pipeline mojo is happy.
	//
	vtkSmartPointer < vtkSphereSource > sphere = vtkSmartPointer
			< vtkSphereSource > ::New();
	sphere->SetThetaResolution(100);
	sphere->SetPhiResolution(50);

	// The mapper is responsible for pushing the geometry into the graphics
	// library. It may also do color mapping, if scalars or other attributes
	// are defined.
	//
	vtkSmartPointer < vtk441MapperPart1 > win1Mapper = vtkSmartPointer
			< vtk441MapperPart1 > ::New();
	win1Mapper->SetInputConnection(sphere->GetOutputPort());

	vtkSmartPointer < vtkActor > win1Actor = vtkSmartPointer < vtkActor
			> ::New();
	win1Actor->SetMapper(win1Mapper);

	vtkSmartPointer < vtkRenderer > ren1 = vtkSmartPointer < vtkRenderer
			> ::New();

	vtkSmartPointer < vtkRenderWindow > renWin = vtkSmartPointer
			< vtkRenderWindow > ::New();
	renWin->AddRenderer(ren1);
	ren1->SetViewport(0, 0, 1, 1);

	vtkSmartPointer < vtkRenderWindowInteractor > iren = vtkSmartPointer
			< vtkRenderWindowInteractor > ::New();
	iren->SetRenderWindow(renWin);

	// Add the actors to the renderer, set the background and size.
	//
	bool doWindow1 = true;
	if (doWindow1)
		ren1->AddActor(win1Actor);
	ren1->SetBackground(0.0, 0.0, 0.0);

	// Set up the lighting.
	//
	vtkRenderer *rens[1] = { ren1 };
	for (int i = 0; i < 1; i++) {
		rens[i]->GetActiveCamera()->SetFocalPoint(0, 0, 0);
		rens[i]->GetActiveCamera()->SetPosition(0, 0, 70);
		rens[i]->GetActiveCamera()->SetViewUp(0, 1, 0);
		rens[i]->GetActiveCamera()->SetClippingRange(-120, 120);
		rens[i]->GetActiveCamera()->SetDistance(70);
	}

	// This starts the event loop and invokes an initial render.
	//
	((vtkInteractorStyle *) iren->GetInteractorStyle())->SetAutoAdjustCameraClippingRange(
			0);
	iren->Initialize();
	iren->Start();

	return EXIT_SUCCESS;
}
