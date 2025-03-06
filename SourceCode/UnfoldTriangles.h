#ifndef __UNFOLDTRIANGLES_H__
#define __UNFOLDTRIANGLES_H__
#include "UnfoldingPath.h"
#include "ParticlePositions.h"
#include "Body.h"

namespace unfolding {
	class UnfoldTriangles {
	public:
		UnfoldingPath _unfpt;
		ParticlePositions _partpos;

		int counter;											// Counts how many particles have been placed in the 2D template
		int minimcounter;										// Counts how many smaller triangles have been unfolded
		int extid;												// Counts the number of repeated particles are placed in the 2D template
		bool minim;												// If 1, minimization is run
		bool scar;												// If 1, unfolding is performed upon the given scar trail "trals"
		string so = _partpos.input_particlefile;				// String for input particles file
		vector<pair<pair<int, int>, vector<int>>>allpaths;		// Stores previously calulcated shortest paths between two nodes on the 3D structure
		vector<Vector2d> pos2D;									// Position of particles on flat unfolded plane
		vector<Vector2d> posord;								// Position of particles on flat unfolded plane
		map < int, int> ID2D;									// ID of particles in pos2D vector
		map<int, int>ord;										// Maps each node index to its index in unfolded template (nodes can be repeated, which is why this is needed)
		map<int, int>ordrev;									// Reverse mapping of "ord"
		vector<bool>mapped;										// Flags for particles that have been already mapped
		vector<vector<int>>donetrisglobal;						// Stores all the triangles that have been flattened to 2D
		map<vector<int>, bool>maptris;							// Flags for triangles that have been flattened to 2D
		vector<vector<int>>rep_pts;								// Stores all the repeated particles
		vector<vector<int>>untripts;							// Set of points that are mapped at each iteration of large triangle unfolding
		vector<vector<int>>nountripts;							// Set of points that were inside the large triangles but not mapped at each iteration
		vector<vector<int>>trals;								// Trail of scar to cut along	
		
		UnfoldTriangles(UnfoldingPath UNFPT, ParticlePositions PartPos, bool Minim, bool Scar, vector<vector<int>>TRLS) :_unfpt(UNFPT), _partpos(PartPos), minim(Minim), scar(Scar), trals(TRLS) {
			pos2D.resize(_partpos.n_part, { 0,0 });
			mapped.resize(_partpos.n_part, false);
			rep_pts.resize(_partpos.n_part, vector<int>());
			for (int i = 0; i < _partpos.n_part; i++) {
				ID2D[i] = i;
				rep_pts[i].push_back(i);
			}
			extid = _partpos.n_part;
			counter = 0;
			minimcounter = 0;
			Unfold(scar);
			PrintPts();
			PrintTrs();
			PrintOutline();
		}

		TriPtsDS TriPts(vector<int>T);									// Finds the particles that will be unfolded as part of the given large triangle
		TriPtsDS ScarTriPts(vector<int>T, vector<vector<int>>Trls);		// Finds the particles that will be unfolded if cut along the given trail of scar
		void Unfold(bool isScar);										// Runs the unfolding of smaller triangles sequentially
		void PrintPts();												// Prints the vertices of the large triangles
		void PrintTrs();												// Prints the smaller triangles
		void PrintOutline();											// Prints the outline of the unfolding path
		void PrintEnergy();												// Prints the energy and projected gradient
		void find_paths(vector<vector<int> >& pathsn,
			vector<int>& pathn,
			vector<vector<int>>& parentn,
			int nn, int un,int sm,int tm);					// Function to find shortest path between two nodes by backtracking solution from Breadth First Search algorithm
		void CalcDistortion(string s);						// Runs minimization if asked and calculates metrics to compare distortion in the flat template
	};
}
#endif
