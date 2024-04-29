#ifndef __Particle_Positions_h__
#define __Particle_Positions_h__
#include "Unfolding.h"
namespace unfolding {
	class ParticlePositions {
	public:
		int n_part;								// Number of particles
		string input_particlefile;
		string input_connectivityfile;
		string input3Dcon;
		string input_triangles;
		vector<Vector3d> pos;					// position of particles
		vector<vector<int>> pos_neighbor;		// neighbour nodes of target particles
		vector<vector<int>>MinDisGrp;			// data structure to contain shortest path from one particle to other through connecting edges
		vector<vector<int>> Triangles;
		ParticlePositions(int n, string inputpart, string inputcon, string input3d, string inputtriangles) {
			n_part = n;
			input_particlefile = inputpart;
			input_connectivityfile = inputcon;
			input_triangles = inputtriangles;
			input3Dcon = input3d;
			MinDisGrp.resize(n_part, vector<int>(n_part));
			readParticlePositions();
			readParticleConnectivity();
			if (input_triangles != "NULL")readTriangles();
			calcFloyd();
		}

		void readParticlePositions();
		void readParticleConnectivity();
		void getAllParticles();
		void getParticleID(int i);
		void getParticleNeighbors(int i);
		void readTriangles();
		void calcFloyd();
	};
}

#endif