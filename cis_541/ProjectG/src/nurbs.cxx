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
#include <GL/glu.h>

/*void nurbsError(GLenum errorCode)
{
   const GLubyte *estring;
   estring = gluErrorString(errorCode);
   fprintf (stderr, "Nurbs Error: %s\n", estring);
   exit (0);
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

GLfloat ctlpoints[4][4][3] = {
   {{ -1.5, -1.5, 4.0}, { -0.5, -1.5, 2.0},
    {0.5, -1.5, -1.0}, {1.5, -1.5, 2.0}},
   {{ -1.5, -0.5, 1.0}, { -0.5, -0.5, 3.0},
    {0.5, -0.5, 0.0}, {1.5, -0.5, -1.0}},
   {{ -1.5, 0.5, 4.0}, { -0.5, 0.5, 0.0},
    {0.5, 0.5, 3.0}, {1.5, 0.5, 4.0}},
   {{ -1.5, 1.5, -2.0}, { -0.5, 1.5, -2.0},
    {0.5, 1.5, 0.0}, {1.5, 1.5, -1.0}}
};

class vtk441MapperPart1: public vtk441Mapper {
public:
	static vtk441MapperPart1 *New();

	virtual void RenderPiece(vtkRenderer *ren, vtkActor *act) {
		RemoveVTKOpenGLStateSideEffects();
		SetupLight();

    GLUnurbsObj *theNurb;

    int u, v;
    for (u = 0; u < 4; u++) {
       for (v = 0; v < 4; v++) {
          ctlpoints[u][v][0] = 2.0*((GLfloat)u - 1.5);
          ctlpoints[u][v][1] = 2.0*((GLfloat)v - 1.5);

          if ( (u == 1 || u == 2) && (v == 1 || v == 2))
             ctlpoints[u][v][2] = 3.0;
          else
             ctlpoints[u][v][2] = -3.0;
       }
    }

    GLfloat mat_diffuse[] = { 0.7, 0.7, 0.7, 1.0 };
    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat mat_shininess[] = { 100.0 };

    glClearColor (0.0, 0.0, 0.0, 0.0);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_AUTO_NORMAL);
    glEnable(GL_NORMALIZE);

    theNurb = gluNewNurbsRenderer();
    gluNurbsProperty(theNurb, GLU_SAMPLING_TOLERANCE, 25.0);
    gluNurbsProperty(theNurb, GLU_DISPLAY_MODE, GLU_FILL);

    /*gluNurbsCallback(theNurb, GLU_ERROR,
                     (GLvoid (*)()) nurbsError);*/

    GLfloat knots[8] = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
   int i, j;

   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   glPushMatrix();
   glRotatef(330.0, 1.,0.,0.);
   glScalef (0.5, 0.5, 0.5);

   gluBeginSurface(theNurb);
   gluNurbsSurface(theNurb,
                   8, knots, 8, knots,
                   4 * 3, 3, &ctlpoints[0][0][0],
                   4, 4, GL_MAP2_VERTEX_3);
   gluEndSurface(theNurb);

  glPointSize(5.0);
  glDisable(GL_LIGHTING);
  glColor3f(1.0, 1.0, 0.0);
  glBegin(GL_POINTS);
  for (i = 0; i < 4; i++) {
     for (j = 0; j < 4; j++) {
        glVertex3f(ctlpoints[i][j][0],
                   ctlpoints[i][j][1], ctlpoints[i][j][2]);
     }
  }
  glEnd();
  glEnable(GL_LIGHTING);

   glPopMatrix();
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
