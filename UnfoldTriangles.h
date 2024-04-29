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
		string so = _partpos.input3Dcon;
		vector<pair<vector<int>, vector<int>>>triptsset;
		vector<pair<pair<int, int>, vector<int>>>allpaths;
		vector<Vector2d> pos2D;					// position of particles on flat unfolded plane
		vector<Vector2d> posord;
		map < int, int> ID2D;					// ID of particles in pos2D vector
		map<int, int>ord;
		map<int, int>ordrev;
		vector<bool>mapped;
		vector<vector<int>>donetrisglobal;
		map<vector<int>, bool>maptris;
		vector<vector<int>>rep_pts;
		vector<vector<int>>untripts;
		vector<vector<int>>nountripts;
		int counter;
		int minimcounter;
		int iter;
		int extid;
		bool minim;
		UnfoldTriangles(UnfoldingPath UNFPT, ParticlePositions PartPos, bool Minim) :_unfpt(UNFPT), _partpos(PartPos), minim(Minim) {
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
			Unfold();
			PrintPts();
			//PrintTrs();
			PrintOutline();
			CalcDistortion("ABC");
		}
		TriPtsDS TriPts(vector<int>T);
		void Unfold();
		void PrintPts();
		void PrintTrs();
		void PrintOutline();
		void find_paths(vector<vector<int> >& pathsn,
			vector<int>& pathn,
			vector<vector<int>>& parentn,
			int nn, int un);
		void CalcDistortion(string s);
	};
}
#endif
