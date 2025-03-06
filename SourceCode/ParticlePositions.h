#ifndef __Particle_Positions_h__
#define __Particle_Positions_h__
#include "Unfolding.h"
namespace unfolding {
	class ParticlePositions {
	public:
		int n_part;								// Number of particles
		
		//Details of these strings can be found in "MainUnfold.cpp"
		string input_particlefile;
		string input_connectivityfile;
		string input_triangles;
		string input_smalltris;
		string _folder;

		vector<Vector3d> pos;					// position of particles
		vector<vector<int>> pos_neighbor;		// neighbour nodes of target particles
		vector<vector<int>>MinDisGrp;			// data structure to contain shortest path from one particle to other through connecting edges
		vector<vector<int>> Triangles;			// Data structure to store small triangles if given

		ParticlePositions(int n, string inputpart, string inputcon, string inputtriangles,string inputtri, string folder) {
			n_part = n;
			input_particlefile = inputpart;
			input_connectivityfile = inputcon;
			input_triangles = inputtriangles;
			input_smalltris = inputtri;

			//Initializing data structure to run Floyd Warshall Algorithm
			MinDisGrp.resize(n_part, vector<int>(n_part));
			readParticlePositions();
			readParticleConnectivity();
			if (input_triangles != "NULL")readTriangles();
			calcFloyd();
			_folder = folder;
		}

		//Function names here are self explanatory
		void readParticlePositions();
		void readParticleConnectivity();
		void getAllParticles();
		void getParticleID(int i);
		void getParticleNeighbors(int i);
		void readTriangles();

		void calcFloyd();					// Run Floyd Warshall algorithm to calculate distance between each set of particles
		Vector3d getcenter();				// Get center spot - average position of all particles
		void centerparticles();				// Translate all particles to 0,0,0
	};
}

#endif