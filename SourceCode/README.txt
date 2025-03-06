## Selecting testcase: (MainUnfold.cpp)
 - Choose to run a icosahedron, parabola, or dome by selecting "swtch" accordingly.
 - Modify "mststeps" to determine the number of MST iterations to search for unfolding path.
 - Modify "minim" to toggle between running distortion minimization.

## Selecting unfolding path criteria: (UnfoldingPath.cpp)
 - Select among "metric1/2/3" to determine the criteria to search unfolding path:
   1. Minimum bounding box.
   2. Minimum difference between height and width of bounding box.
   3. Sequential unfolding.

## Select quadratic elements handling criteria: (Body.cpp)
 - Choose geometry at the top to determine which equation to use for midpoint determination for quadratic elements ("parabola"/"refdome").

## Select distortion minimization parameters: (UnfoldTriangles.cpp)
 - Select desired LBFGS parameters.
 - Select how many iterations of unfolding to run.
 - Select when to start printing energy and projected gradients of the LBFGS process.
