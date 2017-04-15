#include <vtkm/io/reader/VTKDataSetReader.h>
#include <vtkm/filter/PointElevation.h>
#include <vtkm/rendering/Actor.h>
#include <vtkm/rendering/Actor.h>
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/MapperRayTracer.h>
#include <vtkm/rendering/Scene.h>
#include <vtkm/rendering/View3D.h>

using std::cout;
using std::endl;

/*vtkm::rendering::View3D *gViewPointer = NULL;

void DisplayCallBack() {
	gViewPointer->Paint();
	glutSwapBuffers();
}*/

int main(int argc, char*argv[]) {
	vtkm::io::reader::VTKDataSetReader reader("hardyglobal.vtk");
	vtkm::cont::DataSet dataset = reader.ReadDataSet();
	vtkm::cont::CoordinateSystem coordinates = dataset.GetCoordinateSystem();
	vtkm::Bounds bounds = coordinates.GetBounds();
	vtkm::filter::PointElevation pefilter;
	pefilter.SetOutputFieldName("elevation");
	pefilter.SetLowPoint(0,bounds.Y.Min,0);
	pefilter.SetHighPoint(0,bounds.Y.Max,0);
	cout << "Starting filter eval" << endl;
	vtkm::filter::ResultField result = pefilter.Execute(dataset, coordinates);
	cout << "Completed filter eval" << endl;
	vtkm::cont::DataSet eval = result.GetDataSet();
	vtkm::rendering::CanvasRayTracer canvas;
	vtkm::rendering::MapperRayTracer mapper;
	/*vtkm::rendering::MapperGL mapper;
	vtkm::rendering::CanvasGL canvas;*/
	vtkm::rendering::Actor actor(dataset.GetCellSet(), dataset.GetCoordinateSystem(), dataset.GetField("wind_velocity_magnitude"));
	vtkm::rendering::Scene scene;
	scene.AddActor(actor);
	vtkm::rendering::View3D view(scene, mapper, canvas);
	vtkm::rendering::Camera camera = view.GetCamera();
	camera.SetPosition(vtkm::make_Vec(0.,0.,0.));
	camera.SetLookAt(vtkm::make_Vec(0.,10.,0.));
	camera.SetViewUp(vtkm::make_Vec(0.,0.,1.));
	camera.SetFieldOfView(90.);
	camera.SetClippingRange(0.01,300.);
	view.SetCamera(camera);
	view.Initialize();
	cout << "Initialize success" << endl;
	view.Paint();
	view.SaveAs("MyVis.ppm");
	return 0;
}
