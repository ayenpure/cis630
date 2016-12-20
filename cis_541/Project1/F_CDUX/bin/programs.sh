./qualitywork 0 &> work0 & ./qualitytriangles 0 &> triangles0 & ./qualitycamera 0 &> camera0

wait
echo "config 0 completed"
rm camera_*

./qualitywork 1 &> work1 & ./qualitytriangles 1 &> triangles1 & ./qualitycamera 1 &> camera1

wait
echo "config 1 completed"
rm camera_*

./qualitywork 2 &> work2 & ./qualitytriangles 2 &> triangles2 & ./qualitycamera 2 &> camera2

wait
echo "config 2 completed"
rm camera_*

./qualitywork 3 &> work3 & ./qualitytriangles 3 &> triangles3 & ./qualitycamera 3 &> camera3

wait
echo "config 3 completed"
rm camera_*

./qualitywork 4 &> work4 & ./qualitytriangles 4 &> triangles4 & ./qualitycamera 4 &> camera4

wait
echo "config 4 completed"
rm camera_*
