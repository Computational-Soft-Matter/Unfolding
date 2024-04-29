#ifndef __Defects_Positions_h__
#define __Defects_Positions_h__
#include "Unfolding.h"
#include "ParticlePositions.h"
namespace unfolding {
	class DefectsPositions {
	public:
		ParticlePositions _partpos;
		int n_def;			   // Number of defects
		Real _rc;              // Cut-off distance for particle search
		Real _delr;	           // Window width
		Real _dr;		       // Radius change
		Real _initialr; 	   // Starting radius
		vector<int>defect_pos;				    // position of particles that has defected valence
		vector<vector<int>>defect_neighbor;	    // neighboring defects of defects
		string _inputdefectfile = "NULL";			//Defects' locations are given in an inputfile
		map<int, int>id_defects;		//Maps the actual defect particle numbers to its position in defect_pos

		DefectsPositions(ParticlePositions PartPos, Real R_c, Real Delr, Real Dr, Real InitialR,string inputdef) :_partpos(PartPos), _rc(R_c), _delr(Delr), _dr(Dr), _initialr(InitialR), _inputdefectfile(inputdef){
			findDefects();
			findConnectivity();
			mapdefectid();
			getNofDefects();
			getDefectsLocations();
		};
		DefectsPositions(ParticlePositions PartPos, string InputF) :_partpos(PartPos), _inputdefectfile(InputF) {
			findDefects();
		};

		void findDefects();
		void findConnectivity();
		void mapdefectid();
		int getNofDefects();
		void getDefectsLocations();
		void getAllDefectsNeighbors();
		void getDefectsNeighborsID(int i);
	};
}
#endif