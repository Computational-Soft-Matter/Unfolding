#pragma once
#pragma once
#ifndef __Body_h__
#define __Body_h__

#include "Unfolding.h"

namespace unfolding {
	class Body {
	public:
		Body(vector<Vector2d>Positions, vector<vector<int>>Tris, map<int, int>ord3d, string MainStructure, bool swtch, Real tolTri, Real tolEdg,int npart) :_Positions2D(Positions), _Tris(Tris), _ord3D(ord3d), _MS(MainStructure), _swtch(swtch), _toltri(tolTri), _toledg(tolEdg),_npart(npart){
			n = _Positions2D.size();
			//cout << "DOF: "<< n << endl;
			_Energy = 0.0;
			_PosOrdrd.resize(2 * n, 0.0);
			_GrdOrdrd.resize(2 * n, 0.0);
			for (int i = 0; i < n; i++) {
				_PosOrdrd[2 * i + 0] = _Positions2D[i][0];
				_PosOrdrd[2 * i + 1] = _Positions2D[i][1];
			}
			//print_vector(_PosOrdrd);
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
			//cout << _Edges.size() << " " << _EdgeTri.size() << endl;
			getMPos(_MS);
			if (_swtch == 1) {
				bndred.resize(_Edges.size(), 0);
				findBndrEdgs();
				//checkDiffEdg(_toledg);
				//checkDiffTri(_toltri);
				//checkconsistency(1e-10, 1e-8);
			}
			else {
				//checkconsistencyNeoHookean2D(1e-6, 1e-6);
				setNeoHookeanGrd2D();
				//setEdgGrd();
				//setTriGrd();
			}
		};
		int n;
		int _npart;
		bool _swtch;
		Real _Energy;
		Real _toltri;
		Real _toledg;
		vector<Real>Err;
		string _MS;
		map < int, int>_ord3D;
		vector<bool>bndred;
		vector<Vector3d>_MPos3D;
		vector<double>_PosOrdrd;
		vector<double>_GrdOrdrd;
		vector<Vector2d>_Positions2D;
		vector<vector<int>>_Tris;
		vector<pair<int, int>>_Edges;
		vector<int>_EdgeTri;
		vector<double>posi;
		void getMPos(string _MS);
		vector<double> getGrd();
		vector<double> getPos();
		void resetGrd();
		void setEdgGrd();
		void setNeoHookeanGrd();
		void setNeoHookeanGrd2D();
		pair<vector< pair<vector<int>, Real>>, vector< pair<vector<int>, Real>>> plotNeoHookeanGrd2D();
		void setTriGrd();
		void resetPos(vector<double>pos);
		void resetEnergy();
		double getenrg();
		void checkconsistency(Real h, Real tol);
		void checkconsistencyNeoHookean(Real h, Real tol);
		void checkconsistencyNeoHookean2D(Real h, Real tol);
		vector<Real> checkDiffEdg(int wt);
		vector< pair<vector<int>, Real>> checkDiffTri();
		Real checkTotalSurfaceTri();
		void findBndrEdgs();
	};
}

#endif