#ifndef __Particle_Positions_h__
#define __Particle_Positions_h__
#include "Unfolding.h"
namespace unfolding {
	class ParticlePositions {
	public:
		int n_part;								// Number of particles
		string input_particlefile;
		string input_connectivityfile;
		vector<Vector3d> pos;					// position of particles
		vector<vector<int>> pos_neighbor;		// neighbour nodes of target particles
		vector<vector<int>>MinDisGrp;			// data structure to contain shortest path from one particle to other through connecting edges

		ParticlePositions(int n, string inputpart, string inputcon) {
			n_part = n;
			input_particlefile = inputpart;
			input_connectivityfile = inputcon;
			MinDisGrp.resize(n_part, vector<int>(n_part));
			readParticlePositions();
			readParticleConnectivity();
			calcFloyd();
		}

		void readParticlePositions();
		void readParticleConnectivity();
		void getAllParticles();
		void getParticleID(int i);
		void getParticleNeighbors(int i);
		void calcFloyd();
	};
}

#endif