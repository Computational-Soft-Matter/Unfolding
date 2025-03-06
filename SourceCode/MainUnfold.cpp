#include "Unfolding.h"
#include "ParticlePositions.h"
#include "DefectsPositions.h"
#include "Triangulation.h"
#include "UnfoldingPath.h"
#include "UnfoldTriangles.h"
#include "Test.h"
#include <LBFGSB.h>
#include <Eigen>

using namespace unfolding;

int main() {
	//Number of particles in the given 3D structure
	int numberof_particles = 161;	

	//Input file for particle positions
	string inputparts;		

	//Input file for particle connectivty
	string inputcon;	

	//Input file for larger triangle connectivity
	//If this file is given, then the code do not search for defects or build large triangles from defects
	//Otherwise, it should be set to "NULL"
	string inputlargetris;	

	//Input file for defects positions (This should be given if 'inputtriscon' is given)
	//If this file is given, then the code do not search for defects or build large triangles from defects
	//Otherwise, it should be set to "NULL"
	string inputdef;

	//Input file for smaller triangles connectivity
	//Should be set to "NULL" if not given
	string inputsmalltris;

	//Folder name where the output is stored
	string folder;

	//swtch determines which geometry to work with
	// 0 - Icosahedron
	// 1 - Dome
	// 2 - Parabola
	int swtch = 2;

	if (swtch == 0) {										
		/*inputparts = "Tsevenexp.vtk";
		inputcon = "T7expdata.txt";						    
		inputlargetris = "NULL";
		inputdef = "NULL";
		inputsmalltris = "NULL";
		folder = "T7Metric1";*/

		inputparts = "Tthreeexp.vtk";
		inputcon = "Tthreeexpdata.txt";
		inputlargetris = "NULL";
		inputdef = "NULL";
		inputsmalltris = "NULL";
		folder = "T3Metric1";

		/*inputparts = "TsvnParticlesrefined.vtk";
		inputcon = "TsvnParticlesrefineddata.txt";						
		inputlargetris = "Tsvntris.txt";
		inputdef = "Tsvndefects.txt";
		inputsmalltris = "NULL";
		folder = "T7RefinedSpecial";

		inputparts = "Tsvnexp.vtk";						
		inputcon = "T7dataexp.vtk";						
		inputlargetris = "Tsvntris.txt";
		inputdef = "Tsvndefects.txt";
		inputsmalltris = "NULL";
		folder = "T7RefinedNewModel";

		inputparts = "TthreeRefined.vtk";						
		inputcon = "T3refineddata.txt";						
		inputlargetris = "NULL";
		inputdef = "NULL";
		inputsmalltris = "NULL";
		folder = "T3Refined";*/
	}
	else if (swtch == 1) {
		/*inputparts = "RefinedSmallDome.vtk";
		inputcon = "RefinedSmallDomedata.txt";						
		inputlargetris = "IcosTri.txt";
		inputdef = "Icos105Defect.txt";
		inputsmalltris = "NULL";
		folder = "IcosSmallRefinedM1test";
		
		inputparts = "IcosModerateRefined.vtk";						
		inputcon = "IcosModerateRefineddata.txt";						
		inputlargetris = "Icos193Tri2.txt";
		inputdef = "Icos193Defect2.txt";
		inputsmalltris = "NULL";
		folder = "IcosModerateRefinedM1";

		inputparts = "IcosModerateDomeCut.vtk";						
		inputcon = "Icos12Dome193Cut2data.txt";						
		inputlargetris = "Icos193Tri2.txt";
		inputdef = "Icos193Defect2.txt";
		inputsmalltris = "NULL";
		folder = "IcosModerateMetric3";*/

		inputparts = "IcosDomeeksho.vtk";						
		inputcon = "Icos12Dome105data.txt";						
		inputlargetris = "IcosTri.txt";
		inputdef = "Icos105Defect.txt";
		inputsmalltris = "NULL";
		folder = "IcosSmallTestLinear";
	}
	else if (swtch == 2) {
		/*inputparts = "finalthirteen.vtk";
		inputcon = "13finaldata.txt";						
		inputlargetris = "1finaltris.txt";
		inputdef = "1finaldefects.txt";
		inputsmalltris = "NULL";
		folder = "parabolascarM1";*/
		
		inputparts = "finalseven.vtk";
		inputcon = "finalsevendata.txt";						
		inputlargetris = "7finaltris.txt";
		inputdef = "7finaldefects.txt";
		inputsmalltris = "NULL";
		folder = "testrefineparabolamulticutM3";

		/*inputparts = "onefinal.vtk";						
		inputcon = "1finaldata.txt";					
		inputlargetris = "1finaltris.txt";
		inputdef = "1finaldefects.txt";
		inputsmalltris = "NULL";
		folder = "parabolasingletestnorefine";*/

		/*inputparts = "refinedparabola.vtk";
		inputcon = "refinedparaboladata.txt";						
		inputlargetris = "1finaltris.txt";
		inputdef = "1finaldefects.txt";
		inputsmalltris = "NULL";
		folder = "parabolasingletestREFINE";

		inputparts = "finalsevenrefined.vtk";
		inputcon = "finalsevenrefineddata.txt";
		inputlargetris = "7finaltris.txt";
		inputdef = "7finaldefects.txt";
		inputsmalltris = "NULL";
		folder = "multicutparabolarefined";*/
	}
	
	bool Minim = 1;				// 1 if minimization is on, 0 if not
	int mststeps = 1;			// Number of times MST algorithm will be run to find unfolding path
	
	//Search Parameters for Defects Connectivity
	Real RC = 4.0, DELR = 2.0, DR = 1.0, INITIALR = 1.0;
	ParticlePositions Part(numberof_particles,inputparts,inputcon, inputlargetris,inputsmalltris, folder);
	
	// This data is needed only if structure is cut along a scar
	bool scars = 0;			    // 1 if cut along scar, 0 if not
	vector<vector<int>>scartrls;
	scartrls.push_back({ 0,8,72,88,100,28,23,25,21,31,1 });
	scartrls.push_back({ 0,61,10,14,17,141,121,129,113,153,33 });
	scartrls.push_back({ 0,62,69,85,97,142,122,130,114,154,34 });
	scartrls.push_back({ 0,63,70,86,98,143,123,131,115,155,35 });
	scartrls.push_back({ 0,64,71,87,99,144,124,132,116,156,36 });


	DefectsPositions Dfct(Part,RC, DELR, DR, INITIALR, inputdef);
	Triangulation Tri(Part, Dfct, inputlargetris);
	UnfoldingPath UPath(Tri,Part,mststeps);
	UnfoldTriangles UnTri(UPath, Part, Minim, scars, scartrls);
}