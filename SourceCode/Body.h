#pragma once
#pragma once
#ifndef __Body_h__
#define __Body_h__

#include "Unfolding.h"

namespace unfolding {
	class Body {
	public:
		int n;										// Total number of DOF's
		int _npart;									// Total number of particles in the 3D structure
		bool _swtch;								// if 1, then calculates metrics of reference vs deformed; if 0, then sets gradients and energies
		string _MS;									// Filename for reading particle positions 
		Real _Energy;								// Total Energy	
		vector<bool>bndred;							// List of edges at the boundary of the unfolded template
		map<int, int>_ord3D;						// Maps each node index to its index in unfolded template (nodes can be repeated, which is why this is needed)
		map<string, int>_Midpts;					// Maps the mid points of edges where string defines which edge and int defines the index of the mid point
		vector<Vector3d>_MPos3D;					// Positions of particles in the original 3D structure
		vector<double>_PosOrdrd;					// Nodal positions for DOF's
		vector<double>_GrdOrdrd;					// Gradient value for DOF's
		vector<Vector2d>_Positions2D;				// Nodal positions for each unfolded points
		vector<vector<int>>_Tris;					// List of triangles in the model expressed as tuple of vertices
		vector<pair<int, int>>_Edges;				// List of edges in the model expressed as pair of vertices
		vector<int>_EdgeTri;						// List of triangles that each edge belong to
		vector< pair<vector<int>, Real>>DiffTris;	// Container to store area distortion for triangles in reference vs deformed
		vector< pair<vector<int>, Real>>TRCS;		// Container to store trace parameter in Evans energy for triangles in reference vs deformed
		vector< pair<vector<int>, Real>>JS;			// Container to store jacobian parameter in Evans energy for triangles in reference vs deformed

		Body(vector<Vector2d>Positions, vector<vector<int>>Tris, map<int, int>ord3d, map<string,int>midpts, string MainStructure, bool swtch,int npart) :_Positions2D(Positions), _Tris(Tris), _ord3D(ord3d), _Midpts(midpts), _MS(MainStructure), _swtch(swtch),_npart(npart){
			n = _Positions2D.size();				// Number of DOF's 
			_Energy = 0.0;
			_PosOrdrd.resize(2 * n, 0.0);
			_GrdOrdrd.resize(2 * n, 0.0);
			// This section sets nodal positions for DOF's from unfolded particle locations
			for (int i = 0; i < n; i++) {
				_PosOrdrd[2 * i + 0] = _Positions2D[i][0];
				_PosOrdrd[2 * i + 1] = _Positions2D[i][1];
			}
			// This section builds list of unfolded triangles
			map<pair<int, int>, int>done;
			for (int i = 0; i < _Tris.size(); i++) {
				int a1 = _Tris[i][0];
				int b1 = _Tris[i][1];
				int c1 = _Tris[i][2];
				pair<int, int>p1;
				pair<int, int>p2;
				pair<int, int>p3;
				if (a1 > b1)p1 = make_pair(a1, b1);
				else p1 = make_pair(b1, a1);

				if (b1 > c1)p2 = make_pair(b1, c1);
				else p2 = make_pair(c1, b1);

				if (a1 > c1)p3 = make_pair(a1, c1);
				else p3 = make_pair(c1, a1);

				if (done[p1] != 5) {
					done[p1] = 5;
					_Edges.push_back(p1);
					_EdgeTri.push_back(c1);
				}
				if (done[p2] != 5) {
					done[p2] = 5;
					_Edges.push_back(p2);
					_EdgeTri.push_back(a1);
				}
				if (done[p3] != 5) {
					done[p3] = 5;
					_Edges.push_back(p3);
					_EdgeTri.push_back(b1);
				}
			}
			getMPos(_MS);
			if (_swtch == 1) {
				checkDiffTri(1);	// Calls function to calculate metrics
			}
			else {
				setEvansElasticEnrgGrd(1,0,0,0);	// Calls function to set energy and gradients
			}
		};

		// Self explanatory functions
		void getMPos(string _MS);
		double getenrg();
		vector<double> getGrd();
		vector<double> getPos();
		void resetGrd();
		void resetPos(vector<double>pos);
		void resetEnergy();

		void findBndrEdgs();		 // Find which edges are at the boundary of the unfolded template
		Real checkTotalSurfaceTri(); // Calculates the surface area of the given 3D structure mesh

		void setEdgGrd();			// Update gradients and energy of DOF's based on spring energy at each edge 
		void setTriGrd();			// Update gradients and energy of DOF's based on deformation of each triangle 
		void setNeoHookeanGrd2D();	// Update gradients and energy of DOF's based on Neo Hookean Energy defined in 2D 
		void setEvansElasticEnrgGrd(int scheme, bool checkmaterial, bool checkelement, bool checkbody);	// Update gradients and energy of DOF's based on Evans Elastic energy
		
		void checkconsistencyEdgTri(Real h, Real tol);		   // Numerically checks if the gradients set by 'setEdgGrd' and 'setTriGrd' are consistent at the model level
		void checkconsistencyNeoHookean2D(Real h, Real tol);   // Numerically checks if the gradients set by 'setNeoHookeanGrd2D' are consistent at the material and model level
		
		void checkDiffTri(int scheme);			// Function to calculate difference of triangles' various parameters in deformed and reference configuration
		vector<Real> checkDiffEdg(int wt);		// Function to calculate difference of triangle' edge lengths in deformed and reference configuration

		vector< pair<vector<int>, Real>> getdiff();		// Get difference in triangle area for each triangle
		vector< pair<vector<int>, Real>> gettrcs();		// Get value for trace parameter in Evans Elastic Function for each triangle
		vector< pair<vector<int>, Real>> getjs();		// Get value for jacobian parameter in Evans Elastic Function for each triangle
	};
}

#endif