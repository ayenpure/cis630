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
#include <vtkModelMetadata.h>

int main(int argc, char *argv[]) {
  vtkPolyDataReader *rdr = vtkPolyDataReader::New();
  rdr->SetFileName("buoyancy.0.vtk");
  vtkPolyData *pd = rdr->GetOutput();
  int numTris = pd->GetNumberOfCells();
  vtkPoints *pts = pd->GetPoints();
  vtkCellArray *cells = pd->GetPolys();
  //vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
  //double *color_ptr = var->GetPointer(0);
  vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("buoyancy");
  float *color_ptr = var->GetPointer(0);
  vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
  //float *normals = n->GetPointer(0);
  vtkModelMetadata *meta = vtkModelMetadata::New();
  vtkIdType npts;
  vtkIdType *ptIds;
  return 0;
}
