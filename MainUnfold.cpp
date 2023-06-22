#include "Unfolding.h"
#include "ParticlePositions.h"
#include "DefectsPositions.h"
#include "Triangulation.h"
#include "UnfoldingPath.h"
#include "UnfoldTriangles.h"
#include "Test.h"

using namespace unfolding;

int main() {
	int numberof_particles = 92;
	string inputparts = "T3exp.vtk";						//Input file for particle positions
	string inputcon = "T3expdata.txt";						//Input file for particle connectivty
	string inputdef = "NULL";
	string inputtri = "NULL";
	int mststeps = 1;
	//Search Parameters for Defects Connectivity
	Real RC = 4.0, DELR = 2.0, DR = 1.0, INITIALR = 1.0;
	ParticlePositions Part(numberof_particles,inputparts,inputcon);
	DefectsPositions Dfct(Part,RC, DELR, DR, INITIALR);
	Triangulation Tri(Part, Dfct, inputtri);
	UnfoldingPath UPath(Tri,mststeps);
	UnfoldTriangles UnTri(UPath, Part);
}