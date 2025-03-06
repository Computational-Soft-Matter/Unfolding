#ifndef __Triangulation_h__
#define __Triangulation_h__
#include "Unfolding.h"
#include "ParticlePositions.h"
#include "DefectsPositions.h"
namespace unfolding {
	class Triangulation {
	public:
		ParticlePositions _partpos;
		DefectsPositions _defectpos;
		string _inputtri;						// input file containing the triangles to unfold
		vector<vector<int>>triangles;			// data structure to contain the triangles of defects
		
		Triangulation(ParticlePositions PartPos, DefectsPositions DefPos, string InputTri) : _inputtri(InputTri), _partpos(PartPos), _defectpos(DefPos) {
			buildTriangles();
			getTriangles();
			//shuffleTriangles();
		};

		void buildTriangles();				// Creates triangles among the defects from defect particles' connectivity 
		void getTriangles();				// Prints the triangles
		void shuffleTriangles();			// Randomly shuffles the order of the triangles
	};
}
#endif