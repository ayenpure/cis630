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

// Purpose: returns a 256x3 array of colors
//
/*unsigned char *
GetColorMap(void) {
	unsigned char controlPts[8][3] = { { 71, 71, 219 }, { 0, 0, 91 }, { 0, 255,
			255 }, { 0, 127, 0 }, { 255, 255, 0 }, { 255, 96, 0 },
			{ 107, 0, 0 }, { 224, 76, 76 }, };
	int textureSize = 256;
	unsigned char *ptr = new unsigned char[textureSize * 3];
	int nControlPts = 8;
	double amountPerPair = ((double) textureSize - 1.0) / (nControlPts - 1.0);
	for (int i = 0; i < textureSize; i++) {
		int lowerControlPt = (int) (i / amountPerPair);
		int upperControlPt = lowerControlPt + 1;
		if (upperControlPt >= nControlPts)
			upperControlPt = lowerControlPt; // happens for i == textureSize-1

		double proportion = (i / amountPerPair) - lowerControlPt;
		for (int j = 0; j < 3; j++)
			ptr[3 * i + j] = controlPts[lowerControlPt][j]
					+ proportion
							* (controlPts[upperControlPt][j]
									- controlPts[lowerControlPt][j]);
	}
	return ptr;
}*/

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

GLfloat ctrlpoints[8][3] = {
        {1,1,0}, {2,6,0}, {3,1,0}, {4,6,0},{5,6,0}, {6,1,0}, {7,6,0}, {8,1,0}};

class vtk441MapperPart1: public vtk441Mapper {
public:
	static vtk441MapperPart1 *New();

	virtual void RenderPiece(vtkRenderer *ren, vtkActor *act) {
		RemoveVTKOpenGLStateSideEffects();
		SetupLight();
		glClearColor(0.0, 0.0, 0.0, 0.0);
   	glShadeModel(GL_FLAT);
   	glMap1f(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 8, &ctrlpoints[0][0]);
   	glEnable(GL_MAP1_VERTEX_3);
		int i;
   	glClear(GL_COLOR_BUFFER_BIT);
		glEnable(GL_COLOR_MATERIAL);
   	glColor3f(1.0, 1.0, 1.0);
   	glBegin(GL_LINE_STRIP);
    	for (i = 0; i <= 30; i++)
       	glEvalCoord1f((GLfloat) i/30.0);
   	glEnd();
   	/* The following code displays the control points as dots. */
   	glPointSize(5.0);
   	glColor3f(1.0, 1.0, 0.0);
   	glBegin(GL_POINTS);
      for (i = 0; i < 8; i++)
         glVertex3fv(&ctrlpoints[i][0]);
   	glEnd();
   	glFlush();
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
