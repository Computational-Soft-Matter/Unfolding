#ifndef __UNFOLDTRIANGLES_H__
#define __UNFOLDTRIANGLES_H__
#include "UnfoldingPath.h"
#include "ParticlePositions.h"

namespace unfolding {
	class UnfoldTriangles {
	public:
		UnfoldingPath _unfpt;
		ParticlePositions _partpos;
		vector<pair<vector<int>, vector<int>>>triptsset;
		vector<Vector2d> pos2D;					// position of particles on flat unfolded plane
		map < int, int> ID2D;					// ID of particles in pos2D vector
		vector<bool>mapped;
		vector<vector<int>>rep_pts;
		int extid;
		UnfoldTriangles(UnfoldingPath UNFPT, ParticlePositions PartPos) :_unfpt(UNFPT), _partpos(PartPos) {
			pos2D.resize(_partpos.n_part, { 0,0 });
			mapped.resize(_partpos.n_part, false);
			rep_pts.resize(_partpos.n_part, vector<int>());
			for (int i = 0; i < _partpos.n_part; i++) {
				ID2D[i] = i;
			}
			extid = _partpos.n_part;
			Unfold();
		}
		TriPtsDS TriPts(vector<int>T);
		void Unfold();

	};
}
#endif
