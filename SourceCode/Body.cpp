#include "Body.h"

namespace unfolding {
	Real muval = 1.0;			// Parameters for Evans Elastic Energy used across multiple functions
	Real lambdaval = 1.5;		// Parameters for Evans Elastic Energy used across multiple functions
	Real Force = 1.0;			// Defined force for nodal locations used across multiple functions
	Real TrCnst = 1.0;			// A constant value to control the weight of Trace term in Evans Energy Function 
	bool parabola = 1;			// Turned on if working with parabola geometry
	bool refdome = 0;			// Turned on if working with dome geometry
	void Body::getMPos(string _MS) {
		string myText;
		ifstream MyReadFile(_MS);
		while (getline(MyReadFile, myText)) {
			string word;
			istringstream ss(myText);
			Real temp_pos[3] = {};
			int gp = 0;
			while (ss >> word)
			{
				temp_pos[gp] = stof(word);
				gp++;
				if (gp == 3) {
					_MPos3D.push_back({ temp_pos[0], temp_pos[1], temp_pos[2] });
					gp = 0;
				}
			}
		}
	}
	void Body::resetPos(vector<double>pos) {
		for (int i = 0; i < n; i++) {
			_PosOrdrd[2 * i + 0] = pos[2 * i + 0];
			_PosOrdrd[2 * i + 1] = pos[2 * i + 1];
		}
	}
	void Body::resetEnergy() {
		_Energy = 0.0;
	}
	vector<double> Body::getPos() {
		return _PosOrdrd;
	}
	vector<double> Body::getGrd() {
		return _GrdOrdrd;
	}
	void Body::resetGrd() {
		_GrdOrdrd.resize(2 * n, 0.0);
	}
	double Body::getenrg() {
		return _Energy;
	}
	void Body::setEdgGrd() {
		for (int i = 0; i < _Edges.size(); i++) {
			int a = _Edges[i].first;
			int b = _Edges[i].second;
			int c = _EdgeTri[i];
			Real l = dist(_MPos3D[a % _npart], _MPos3D[b % _npart]);
			int am = _ord3D[a];
			int bm = _ord3D[b];
			Real Al = trianglearea(_MPos3D[a % _npart], _MPos3D[b % _npart], _MPos3D[c % _npart]);
			Vector2d A = { _PosOrdrd[2 * am + 0],_PosOrdrd[2 * am + 1] };
			Vector2d B = { _PosOrdrd[2 * bm + 0],_PosOrdrd[2 * bm + 1] };
			Real d = dist(A, B);
			Real enrg = d - l;
			_Energy += 0.5*pow(enrg, 2.0)*Al;
			_GrdOrdrd[2 * am + 0] += Al*(enrg) * (A[0] - B[0]) / d;
			_GrdOrdrd[2 * am + 1] += Al*(enrg) * (A[1] - B[1]) / d;
			_GrdOrdrd[2 * bm + 0] += Al * (enrg) * (B[0] - A[0]) / d;
			_GrdOrdrd[2 * bm + 1] += Al * (enrg) * (B[1] - A[1]) / d;

		}
	}
	void Body::setTriGrd() {
		for (int i = 0; i < _Tris.size(); i++) {
			int a = _Tris[i][0];
			int b = _Tris[i][1];
			int c = _Tris[i][2];
			int am = _ord3D[a];
			int bm = _ord3D[b];
			int cm = _ord3D[c];
			Vector2d A = { _PosOrdrd[2 * am + 0],_PosOrdrd[2 * am + 1] };
			Vector2d B = { _PosOrdrd[2 * bm + 0],_PosOrdrd[2 * bm + 1] };
			Vector2d C = { _PosOrdrd[2 * cm + 0],_PosOrdrd[2 * cm + 1] };
			Real Al = trianglearea(_MPos3D[a % _npart], _MPos3D[b % _npart], _MPos3D[c % _npart]);
			Real Aa = ((B[0] - A[0]) * (C[1] - A[1]) - (B[1] - A[1]) * (C[0] - A[0]));
			Real Ad = 0.5 * (Aa);
			Real enrg = (pow(Ad, 2.0) - pow(Al, 2.0));
			_Energy += pow(enrg, 2.0);
			_GrdOrdrd[2 * am + 0] += 2.0 * enrg * (Ad) * (B[1] - C[1]);
			_GrdOrdrd[2 * am + 1] += 2.0 * enrg * (Ad) * (C[0] - B[0]);
			_GrdOrdrd[2 * bm + 0] += 2.0 * enrg * (Ad) * (C[1] - A[1]);
			_GrdOrdrd[2 * bm + 1] += 2.0 * enrg * (Ad) * (A[0] - C[0]);
			_GrdOrdrd[2 * cm + 0] += 2.0 * enrg * (Ad) * (A[1] - B[1]);
			_GrdOrdrd[2 * cm + 1] += 2.0 * enrg * (Ad) * (B[0] - A[0]);
		}
	}
	vector< pair<vector<int>, Real>> Body::getdiff() {
		return DiffTris;
	}
	vector< pair<vector<int>, Real>> Body::gettrcs() {
		return TRCS;
	}
	vector< pair<vector<int>, Real>> Body::getjs() {
		return JS;
	}
	void Body::findBndrEdgs() {
		for (int i = 0; i < _Edges.size(); i++) {
			int a = _Edges[i].first;
			int b = _Edges[i].second;
			int bn = 0;
			for (int j = 0; j < _Tris.size(); j++) {
				if (vector_search(_Tris[j], a) == 1 && vector_search(_Tris[j], b) == 1) {
					bn++;
				}
			}
			if (bn == 1) {
				bndred[i] = 1;
			}
		}
	}
	Real Body::checkTotalSurfaceTri() {
		Real total = 0.0;
		int countr = 0;
		for (int i = 0; i < _Tris.size(); i++) {
			int a = _Tris[i][0];
			int b = _Tris[i][1];
			int c = _Tris[i][2];
			Real Al = trianglearea(_MPos3D[a % _npart], _MPos3D[b % _npart], _MPos3D[c % _npart]);
			total += Al;
		}
		cout << "Bad Triangles: " << countr << endl;
		return total;
	}
	
	// Helper functions to calculate derivative and energy for NeoHookean2D energy
	MatrixXd Pmatrix2D(MatrixXd F) {
		MatrixXd CC = F.transpose() * F;
		Real I1 = CC.trace();
		Real I1bar = pow(F.determinant(), -1.0) * I1;
		Real I3 = CC.determinant();
		MatrixXd P = 1.0 * 4.0 * (I3 * I3 - 1.0 / (I3 * I3)) * F.transpose().inverse() + pow(F.determinant(), -1.0) * (2.0*F - F.transpose().inverse() * I1);
		//MatrixXd P = 2.0 * F;
		return P;
	}
	Real FEnergy2D(MatrixXd F) {
		MatrixXd CC = F.transpose() * F;
		Real I1 = CC.trace();
		Real I1bar = pow(F.determinant(), -1.0) * I1;
		Real I3 = CC.determinant();
		Real enrg = 1.0*(I3 * I3 + 1.0 / (I3 * I3) - 2.0) + (I1bar - 2.0);
		//Real enrg = I1;
		return enrg;
	}

	// Helper functions to calculate shape functions and their derivatives for linear(3) and quadratic(6) triangluar elements
	vector<Real> NTri6(Vector2d s) {
		vector<Real>_functions;
		_functions.resize(6, 0.0);
		_functions[3] = 4.0 * s(0) * (1.0 - s(0) - s(1));
		_functions[4] = 4.0 * s(0) * s(1);
		_functions[5] = 4.0 * s(1) * (1.0 - s(0) - s(1));
		_functions[0] = (1.0 - s(0) - s(1)) * (1.0 - 2.0 * s(0) - 2.0 * s(1));
		_functions[1] = s(0) * (2.0 * s(0) - 1.0);
		_functions[2] = s(1) * (2.0 * s(1) - 1.0);
		return _functions;
	}
	vector<Vector2d> DNTri6(Vector2d s) {
		vector<Vector2d>_derivatives;
		_derivatives.resize(6, { 0.0,0.0 });
		_derivatives[0] = { 4.0 * (s(0) + s(1)) - 3.0, 4.0 * (s(0) + s(1)) - 3.0 };
		_derivatives[1] = { 4.0 * s(0) - 1.0, 0.0 };
		_derivatives[2] = { 0.0, 4.0 * s(1) - 1.0 };
		_derivatives[3] = { 4.0-8.0*s(0)-4.0*s(1), -4.0 * s(0) };
		_derivatives[4] = { 4.0 * s(1), 4.0 * s(0) };
		_derivatives[5] = { -4.0 * s(1), 4.0 - 8.0 * s(1) - 4.0 * s(0) };
		return _derivatives;
	}
	vector<Real> NTri3(Vector2d s) {
		vector<Real>_functions;
		_functions.resize(3, 0.0);
		_functions[0] = 1.0 - s(0) - s(1);
		_functions[1] = s(0);
		_functions[2] = s(1);
		return _functions;
	}
	vector<Vector2d> DNTri3(Vector2d s) {
		vector<Vector2d>_derivatives;
		_derivatives.resize(3, { 0.0,0.0 });
		_derivatives[0] = { -1.0, -1.0 };
		_derivatives[1] = { 1.0, 0.0 };
		_derivatives[2] = { 0.0, 1.0 };
		return _derivatives;
	}

	// Helper functions to calculate basis vectors and metric tensor for deformed and reference triangles used in Evans Elastic Energy 
	vector<Vector3d> basiss(vector<Vector3d>nodes, vector<Vector2d>DN) {
		vector<Vector3d> basis;
		basis.resize(2, { 0.0,0.0,0.0 });
		for (int nb = 0; nb < nodes.size(); nb++) {
			basis[0] += DN[nb](0) * nodes[nb];
			basis[1] += DN[nb](1) * nodes[nb];
		}
		return basis;
	}
	MatrixXd metric(vector<Vector3d>_a) {
		MatrixXd _metricTensor(2, 2);
		_metricTensor(0, 0) = _a[0].dot(_a[0]);
		_metricTensor(0, 1) = _a[0].dot(_a[1]);
		_metricTensor(1, 1) = _a[1].dot(_a[1]);
		_metricTensor(1, 0) = _metricTensor(0, 1);
		return _metricTensor;
	}
	MatrixXd metricinv(MatrixXd MT) {
		MatrixXd _metricTensorinv(2, 2);
		Real adj = MT(0, 0) * MT(1, 1) - MT(0, 1) * MT(1, 0);
		_metricTensorinv(0, 0) = MT(1,1)/adj;
		_metricTensorinv(0, 1) = -MT(0, 1) / adj;
		_metricTensorinv(1, 1) = MT(0,0) / adj;
		_metricTensorinv(1, 0) = _metricTensorinv(0, 1);
		return _metricTensorinv;
	}
	Real sqrtdetmetric(MatrixXd M) {
		return sqrt(M(0, 0) * M(1, 1) - M(0, 1) * M(1, 0));
	}
	vector<Vector3d> duals(vector<Vector3d>_a,MatrixXd Minv) {
		vector<Vector3d>dual;
		dual.resize(2, { 0.0,0.0,0.0 });
		dual[0] = Minv(0, 0) * _a[0] + Minv(0, 1) * _a[1];
		dual[1] = Minv(1, 0) * _a[0] + Minv(1, 1) * _a[1];
		return dual;
	}

	// Functions to set Evans Elastic Energy using given deformation gradient
	pair<MatrixXd, Real> EnergyEvans(MatrixXd F, string f, Real _mu, Real _kS,bool print) {
		MatrixXd C(3, 3), P(3, 3),M(3,3);
		for (int I = 0; I < 3; I++) {
			for (int L = 0; L < 3; L++) {
				C(I, L) = 0.0;
				P(I, L) = 0.0;
				M(I, L) = 0.0;
			}
		}
		C = F.transpose() * F;
		Real _trC = C(0, 0) + C(1, 1) + C(2, 2);
		Real trCSquare = 0.0;
		for (int I = 0; I < 3; I++) {
			for (int k = 0; k < 3; k++) {
				trCSquare += C(I, k) * C(k, I);
			}
		}
		Real _J = sqrt((_trC * _trC - trCSquare) / 2.0);
		Real W = _kS * pow(_J - 1.0, 2.0) / 2.0 + TrCnst*0.5 * _mu * (_trC / _J - 2.0);

		if (f == "Energy") {
			return make_pair(P, W);
		}
		else if (f == "Both") {
			P = ((_trC * F - F * C) * _kS * (_J - 1.0) / _J) + TrCnst *((_mu / _J) * F - (_trC * F - F * C)*(0.5 * _mu * _trC / (_J * _J))/_J);
			return make_pair(P, W);
		}
	}
	pair<Real, Real> EnergyEvansEnergies(MatrixXd F, string f, Real _mu, Real _kS) {
		MatrixXd C(3, 3), P(3, 3);
		for (int I = 0; I < 3; I++) {
			for (int J = 0; J < 3; J++) {
				C(I, J) = 0.0;
				P(I, J) = 0.0;
			}
		}
		C = F.transpose() * F;
		Real _trC = C(0, 0) + C(1, 1) + C(2, 2);
		Real trCSquare = 0.0;
		for (int I = 0; I < 3; I++) {
			for (int k = 0; k < 3; k++) {
				trCSquare += C(I, k) * C(k, I);
			}
		}
		Real _J = sqrt((_trC * _trC - trCSquare) / 2.0);
		Real W1 = _kS * pow(_J - 1.0, 2.0) / 2.0;
		Real W2 = TrCnst*0.5 * _mu * (_trC / _J - 2.0);
		return make_pair(W1, W2);
	}

	// Function to set Hooke's energy using given deofrmation gradient
	pair<MatrixXd, Real> EnergyHooke(MatrixXd F, string f, Real mu, Real lambda, bool print) {
		MatrixXd Id(3,3),U(3, 3), EPS(3, 3),S(3,3);
		for (int I = 0; I < 3; I++) {
			for (int L = 0; L < 3; L++) {
				U(I, L) = 0.0;
				S(I, L) = 0.0;
				EPS(I, L) = 0.0;
				if (I == L && I!=2) {
					Id(I, L) = 1.0;
				}
				else {
					Id(I, L) = 0.0;
				}
			}
		}
		U = F - Id;
		//cout << U << endl;
		EPS = 0.5 * (U + U.transpose());
		//cout << EPS << endl;
		Real ETRC = EPS(0, 0) + EPS(1, 1) + EPS(2, 2);
		for (int I = 0; I < 3; I++) {
			for (int L = 0; L < 3; L++) {
				S(I, L) = 2.0 * mu * EPS(I, L);
				if (I == L) {
					S(I, L) += lambda * ETRC;
				}
			}
		}
		Real W = 0.0;
		for (int I = 0; I < 3; I++) {
			for (int L = 0; L < 3; L++) {
				W += S(I, L) * EPS(I, L);
			}
		}
		W = W * 0.5;
		//cout << W << endl;
		if (f == "Energy") {
			return make_pair(S, W);
		}
		else if (f == "Both") {
			return make_pair(S, W);
		}
	}
	pair<MatrixXd, Real> EnergyHookeEPS(MatrixXd EPS, string f, Real mu, Real lambda, bool print) {
		MatrixXd S(3, 3);
		for (int I = 0; I < 3; I++) {
			for (int L = 0; L < 3; L++) {
				S(I, L) = 0.0;
			}
		}
		Real ETRC = EPS(0, 0) + EPS(1, 1) + EPS(2, 2);
		for (int I = 0; I < 3; I++) {
			for (int L = 0; L < 3; L++) {
				S(I, L) = 2.0 * mu * EPS(I, L);
				if (I == L) {
					S(I, L) += lambda * ETRC;
				}
			}
		}
		Real W = 0.0;
		for (int I = 0; I < 3; I++) {
			for (int L = 0; L < 3; L++) {
				W += S(I, L) * EPS(I, L);
			}
		}
		W = W * 0.5;
		//cout << W << endl;
		if (f == "Energy") {
			return make_pair(S, W);
		}
		else if (f == "Both") {
			return make_pair(S, W);
		}
	}

	// Function to check material level consistency for Evans Elastic Energy
	void checkEvansMaterialConsistency(MatrixXd F, Real h, Real tol,Real _mu,Real _kS,string Enrg) {
		for (int t = 0; t < 3; t++) {
			for (int j = 0; j < 3; j++) {
				MatrixXd Fp = F;
				MatrixXd Fm = F;
				pair<MatrixXd, Real> P = EnergyEvans(F, "Both", _mu, _kS,1);
				if(Enrg=="Hooke") P = EnergyHooke(F, "Both", muval, lambdaval, 1);
				//cout << P.second << endl;
				Fp(t, j) += h;
				Fm(t, j) -= h;
				pair<MatrixXd,Real> Ep = EnergyEvans(Fp, "Both",_mu,_kS,1);
				if (Enrg == "Hooke") Ep = EnergyHooke(Fp, "Both", muval, lambdaval, 1);
				pair<MatrixXd, Real> Em = EnergyEvans(Fm, "Both", _mu, _kS,1);
				if (Enrg == "Hooke") Em = EnergyHooke(Fm, "Both", muval, lambdaval, 1);
				Real errM = (Ep.second - Em.second) / (2.0 * h) - P.first(t, j);
				if (abs(errM) > tol)cout << "problem Material " << abs(errM) << endl;
			}
		}
	}
	
	void Body::setNeoHookeanGrd2D() {
		for (int i = 0; i < _Tris.size(); i++) {
			int a = _Tris[i][0];
			int b = _Tris[i][1];
			int c = _Tris[i][2];
			Vector3d D, E, FF;
			vector<pair<int, int>>ptpairs;
			pair<int, int>pr1, pr2, pr3;
			if (a < b)pr1 = make_pair(a, b);
			else make_pair(b, a);
			if (a < c)pr2 = make_pair(a, c);
			else make_pair(c, a);
			if (b < c)pr3 = make_pair(b, c);
			else make_pair(c, b);
			ptpairs.push_back(pr1);
			ptpairs.push_back(pr2);
			ptpairs.push_back(pr3);
			vector<int>pms;
			for (int ptk = 0; ptk < 3; ptk++) {
				int pt1 = ptpairs[ptk].first;
				int pt2 = ptpairs[ptk].second;
				string ptnow = to_string(pt1) + "a" + to_string(pt2);
				pms.push_back(_Midpts[ptnow]);
			}
			int am = _ord3D[a];
			int bm = _ord3D[b];
			int cm = _ord3D[c];
			Real Al = abs(trianglearea(_MPos3D[a % _npart], _MPos3D[b % _npart], _MPos3D[c % _npart]));
			Real b1, c2, c1;
			Vector3d A = { _PosOrdrd[2 * am + 0],_PosOrdrd[2 * am + 1],0.0 };
			Vector3d B = { _PosOrdrd[2 * bm + 0],_PosOrdrd[2 * bm + 1],0.0 };
			Vector3d C = { _PosOrdrd[2 * cm + 0],_PosOrdrd[2 * cm + 1],0.0 };
			D = { _PosOrdrd[2 * pms[0] + 0],_PosOrdrd[2 * pms[0] + 1],0.0 };
			E = { _PosOrdrd[2 * pms[1] + 0],_PosOrdrd[2 * pms[1] + 1],0.0 };
			FF = { _PosOrdrd[2 * pms[2] + 0],_PosOrdrd[2 * pms[2] + 1],0.0 };
			cout << A << "PT"<< endl;
			cout << B << "PT" << endl;
			cout << C << "PT" << endl;
			cout << D << "PT" << endl;
			cout << E << "PT" << endl;
			cout << FF << "PT" << endl;
			Vector3d AB = B - A;
			Vector3d AC = C - A;
			Vector3d A3d = { _MPos3D[a % _npart][0],_MPos3D[a % _npart][1],_MPos3D[a % _npart][2] };
			Vector3d B3d = { _MPos3D[b % _npart][0],_MPos3D[b % _npart][1],_MPos3D[b % _npart][2] };
			Vector3d C3d = { _MPos3D[c % _npart][0],_MPos3D[c % _npart][1],_MPos3D[c % _npart][2] };
			Vector3d AB3d = B3d - A3d;
			Vector3d AC3d = C3d - A3d;
			Real fact = 1.0;
			if ((AB.cross(AC)[2]) / abs(AB.cross(AC)[2]) != (AB3d.cross(AC3d)[0]) / abs(AB3d.cross(AC3d)[0])) {
				//fact = -1.0;
			}
			MatrixXd F(2, 2);
			Real N1x, N1y, N2x, N2y, N3x, N3y;
			bool change = 0;
			//Real fact = 1.0;
			if (AB.cross(AC)[2] < 0) {
				fact = -1.0;
			}
			vector<Vector3d>ABC = { A,B,C };
			Real dism = dist(_MPos3D[a % _npart], _MPos3D[b % _npart]);
			Real dis1 = dist(_MPos3D[a % _npart], _MPos3D[c % _npart]);
			Real dis2 = dist(_MPos3D[b % _npart], _MPos3D[c % _npart]);
			Real ah = (dis1 * dis1 - dis2 * dis2 + dism * dism) / (2 * dism);
			Real ha = sqrt(dis1 * dis1 - ah * ah);
			Vector3d rootpos = { 0.0,0.0,0.0 };
			Vector3d leverpos = {dism,0.0,0.0};
			Vector3d PP0 = { rootpos[0] + ah * (leverpos[0] - rootpos[0]) / dism,rootpos[1] + ah * (leverpos[1] - rootpos[1]) / dism,0 };
			Vector3d PP1 = { PP0[0] + ha * (leverpos[1] - rootpos[1]) / dism,PP0[1] - ha * (leverpos[0] - rootpos[0]) / dism ,0 };
			Vector3d PP2 = { PP0[0] - ha * (leverpos[1] - rootpos[1]) / dism,PP0[1] + ha * (leverpos[0] - rootpos[0]) / dism ,0 };
			if (fact == 1.0) {
				if (PP1[1] > 0.0) {
					c1 = PP1[0];
					c2 = PP1[1];
				}
				else {
					c1 = PP2[0];
					c2 = PP2[1];
				}
			}
			else {
				if (PP1[1] < 0.0) {
					c1 = PP1[0];
					c2 = PP1[1];
				}
				else {
					c1 = PP2[0];
					c2 = PP2[1];
				}
			}
			b1 = dism;
			//b1 = 0.34;
			//c1 = 0.34 / 2.0;
			//c2 = fact * b1 * sqrt(3.0) / 2.0;
			//c2 = fact*2.0 * Al / b1;
			//c1 = sqrt(pow(dist(_MPos3D[a % _npart], _MPos3D[c % _npart]), 2.0) - c2 * c2);
			N1x = -(1.0 / b1);
			N1y = c1 / (b1 * c2) - (1 / c2);
			N2x = -N1x;
			N2y = -(c1 / (b1 * c2));
			N3x = 0.0;
			N3y = 1.0 / c2;
			F(0, 0) = N1x * A[0] + N2x * B[0];
			F(0, 1) = N1y * A[0] + N2y * B[0] + N3y * C[0];
			F(1, 0) = N1x * A[1] + N2x * B[1];
			F(1, 1) = N1y * A[1] + N2y * B[1] + N3y * C[1];
			_Energy += FEnergy2D(F) * Al;
			MatrixXd P = Pmatrix2D(F);
			vector<vector<Real>>Derivatives;
			Derivatives.push_back({ 0.0,0.0 });
			Derivatives.push_back({ 0.0,0.0 });
			Derivatives.push_back({ 0.0,0.0 });
			Derivatives[0][0] = Al * (P(0, 0) * N1x + P(0, 1) * N1y);
			Derivatives[0][1] = Al * (P(1, 0) * N1x + P(1, 1) * N1y);
			Derivatives[1][0] = Al * (P(0, 0) * N2x + P(0, 1) * N2y);
			Derivatives[1][1] = Al * (P(1, 0) * N2x + P(1, 1) * N2y);
			Derivatives[2][0] = Al * (P(0, 0) * N3x + P(0, 1) * N3y);
			Derivatives[2][1] = Al * (P(1, 0) * N3x + P(1, 1) * N3y);
			if (change == 1) {
				_GrdOrdrd[2 * am + 0] += Derivatives[0][0];
				_GrdOrdrd[2 * am + 1] += Derivatives[0][1];
				_GrdOrdrd[2 * bm + 0] += Derivatives[2][0];
				_GrdOrdrd[2 * bm + 1] += Derivatives[2][1];
				_GrdOrdrd[2 * cm + 0] += Derivatives[1][0];
				_GrdOrdrd[2 * cm + 1] += Derivatives[1][1];
			}
			else {
				_GrdOrdrd[2 * am + 0] += Derivatives[0][0];
				_GrdOrdrd[2 * am + 1] += Derivatives[0][1];
				_GrdOrdrd[2 * bm + 0] += Derivatives[1][0];
				_GrdOrdrd[2 * bm + 1] += Derivatives[1][1];
				_GrdOrdrd[2 * cm + 0] += Derivatives[2][0];
				_GrdOrdrd[2 * cm + 1] += Derivatives[2][1];
			}
		}
	}
	void Body::setEvansElasticEnrgGrd(int scheme,bool checkmaterial,bool checkelement, bool checkbody) {
		// Checks material level consistency
		if (checkmaterial) {
			MatrixXd F(3, 3);
			F(0, 0) = 1.25;
			F(0, 1) = 0.2;
			F(0, 2) = 0.0;
			F(1, 0) = 0.2;
			F(1, 1) = 0.925;
			F(1, 2) = 0.0;
			F(2, 0) = 0.0;
			F(2, 1) = 0.0;
			F(2, 2) = 0.0;
			checkEvansMaterialConsistency(F, 1e-7, 1e-6, 1.0, 1.0,"Hoo");
		}
		
		// Necessary members to define force at nodes
		Real pts = sqrt(_npart);
		if (scheme)pts += (pts-1);
		Real L = 4.0 / (pts-1.0);
		set<int>doneforcing;

		// Iterates through trianglular elements to set energy and gradients
		for (int trr = 0; trr < _Tris.size(); trr++) {
			Real _mu = 1.0;
			Real _kS = 1.0;
			int a = _Tris[trr][0];
			int b = _Tris[trr][1];
			int c = _Tris[trr][2];
			Vector3d D, E, FF;
			vector<pair<int, int>>ptpairs;
			pair<int, int>pr1, pr2, pr3;
			if (a < b)pr1 = make_pair(a, b);
			else pr1 = make_pair(b, a);
			if (b < c)pr2 = make_pair(b, c);
			else pr2 = make_pair(c, b);
			if (a < c)pr3 = make_pair(a, c);
			else pr3 = make_pair(c, a);
			ptpairs.push_back(pr1);
			ptpairs.push_back(pr2);
			ptpairs.push_back(pr3);
			vector<int>pms;
			for (int ptk = 0; ptk < 3; ptk++) {
				int pt1 = ptpairs[ptk].first;
				int pt2 = ptpairs[ptk].second;
				string ptnow = to_string(pt1) + "a" + to_string(pt2);
				pms.push_back(_Midpts[ptnow]);
			}
			int am = _ord3D[a];
			int bm = _ord3D[b];
			int cm = _ord3D[c];
			vector<int>ords;
			ords = { am,bm,cm };
			vector<int>pntpos;
			vector<Vector2d>_quadpoints;
			vector<Real>_quadweights;
			vector<Vector3d>nodescur;
			vector<Vector3d>nodesref;
			vector<Real>N;
			vector<Vector2d>DN;
			Vector3d A = { _PosOrdrd[2 * am + 0],_PosOrdrd[2 * am + 1],0.0 };
			Vector3d B = { _PosOrdrd[2 * bm + 0],_PosOrdrd[2 * bm + 1],0.0 };
			Vector3d C = { _PosOrdrd[2 * cm + 0],_PosOrdrd[2 * cm + 1],0.0 };
			/*_PosOrdrd[2 * pms[0] + 0] = (A[0] + B[0]) / 2.0;
			_PosOrdrd[2 * pms[0] + 1] = (A[1] + B[1]) / 2.0;
			_PosOrdrd[2 * pms[1] + 0] = (B[0] + C[0]) / 2.0;
			_PosOrdrd[2 * pms[1] + 1] = (B[1] + C[1]) / 2.0;
			_PosOrdrd[2 * pms[2] + 0] = (A[0] + C[0]) / 2.0;
			_PosOrdrd[2 * pms[2] + 1] = (A[1] + C[1]) / 2.0;*/
			D = { _PosOrdrd[2 * pms[0] + 0],_PosOrdrd[2 * pms[0] + 1],0.0 };
			E = { _PosOrdrd[2 * pms[1] + 0],_PosOrdrd[2 * pms[1] + 1],0.0 };
			FF = { _PosOrdrd[2 * pms[2] + 0],_PosOrdrd[2 * pms[2] + 1],0.0 };
			vector<Vector3d>NPs; 
			Vector3d NP1, NP2, NP3;
			NP1 = (_MPos3D[a % _npart] + _MPos3D[b % _npart]) / 2.0;
			NP2 = (_MPos3D[b % _npart] + _MPos3D[c % _npart]) / 2.0;
			NP3 = (_MPos3D[c % _npart] + _MPos3D[a % _npart]) / 2.0;
			NPs = {NP1,NP2,NP3};
			
			// This section defines force at node level (element level)
			for (int IP = 0; IP < 3; IP++) {
				Real factor = L;
				if (scheme == 0) {
					if (_MPos3D[_Tris[trr][IP] % _npart][0] == 4.0 && doneforcing.find(ords[IP]) == doneforcing.end()) {
						doneforcing.insert(ords[IP]);
						if (_MPos3D[_Tris[trr][IP] % _npart][1] == 0.0 || _MPos3D[_Tris[trr][IP] % _npart][1] == 4.0)factor = L / 2.0;
						//cout << ords[IP] <<" "<< Force*factor << endl;
						Real Finalforce = Force * factor;
						Vector2d A = { _MPos3D[_Tris[trr][IP] % _npart][0],_MPos3D[_Tris[trr][IP] % _npart][1] };
						Vector2d B = { _PosOrdrd[2 * ords[IP]],_PosOrdrd[2 * ords[IP] + 1] };
						//_Energy += -Finalforce * (_PosOrdrd[2 * ords[IP]] - _MPos3D[_Tris[trr][IP] % _npart][0]);
						//_GrdOrdrd[2 * ords[IP]] += -Finalforce;
					}
				}
				factor = L;
				if (scheme) {
					if (_MPos3D[_Tris[trr][IP] % _npart][0] == 4.0 && doneforcing.find(ords[IP]) == doneforcing.end()) {
						doneforcing.insert(ords[IP]);
						if (_MPos3D[_Tris[trr][IP] % _npart][1] == 0.0 || _MPos3D[_Tris[trr][IP] % _npart][1] == 4.0)factor = L / 6.0;
						else factor = L / 3.0;
						//cout << ords[IP] <<" "<< Force*factor << endl;
						Real Finalforce = Force * factor;
						Vector2d A = { _MPos3D[_Tris[trr][IP] % _npart][0],_MPos3D[_Tris[trr][IP] % _npart][1] };
						Vector2d B = { _PosOrdrd[2 * ords[IP]],_PosOrdrd[2 * ords[IP] + 1] };
						//_Energy += -Finalforce * (_PosOrdrd[2 * ords[IP]] - _MPos3D[_Tris[trr][IP] % _npart][0]);
						//_GrdOrdrd[2 * ords[IP]] += -Finalforce;
					}
					if (NPs[IP][0] == 4.0 && doneforcing.find(pms[IP]) == doneforcing.end()) {
						doneforcing.insert(pms[IP]);
						//cout << NPs[IP][1]<<" scheme " << endl;
						
						factor = 2.0*L / 3.0;
						//cout <<pms[IP]<< " "<< Force * factor << endl;
						Real Finalforce = Force * factor;
						Vector2d A = { NPs[IP][0], NPs[IP][1] };
						Vector2d B = { _PosOrdrd[2 * pms[IP]],_PosOrdrd[2 * pms[IP] + 1] };
						//_Energy += -Finalforce * (_PosOrdrd[2 * pms[IP]] - NPs[IP][0]);
						//_GrdOrdrd[2 * pms[IP]] += -Finalforce;
					}
				}
			}
			
			// Calculates midpoints for parabola, dome cases
			
			if (parabola) {
				NP1[0] = (NP1[1] * NP1[1] + NP1[2] * NP1[2]) / 8.0;
				NP2[0] = (NP2[1] * NP2[1] + NP2[2] * NP2[2]) / 8.0;
				NP3[0] = (NP3[1] * NP3[1] + NP3[2] * NP3[2]) / 8.0;
			}
			vector<Vector3d>finalmids;
			for (int mids = 0; mids < 3; mids++) {
				Vector3d NPnow = NPs[mids];
				Real Radius = sqrt(NPnow[0] * NPnow[0] + NPnow[1] * NPnow[1] + NPnow[2] * NPnow[2]);
				Real theta = acos(NPnow[2] / Radius);
				Real phi = atan2(NPnow[1], NPnow[0]);
				finalmids.push_back({2.0*sin(theta)*cos(phi),2.0 * sin(theta) * sin(phi) ,2.0 * cos(theta) });
			}
			Real Aref = abs(trianglearea(_MPos3D[a % _npart], _MPos3D[b % _npart], _MPos3D[c % _npart]));
			Real dism = dist(_MPos3D[a % _npart], _MPos3D[b % _npart]);
			Real dis1 = dist(_MPos3D[a % _npart], _MPos3D[c % _npart]);
			Real dis2 = dist(_MPos3D[b % _npart], _MPos3D[c % _npart]);
			Real ah = (dis1 * dis1 - dis2 * dis2 + dism * dism) / (2 * dism);
			Real ha = sqrt(dis1 * dis1 - ah * ah);
			Vector3d rootpos = { 0.0,0.0,0.0 };
			Vector3d leverpos = { dism,0.0,0.0 };
			Vector3d PP0 = { rootpos[0] + ah * (leverpos[0] - rootpos[0]) / dism,rootpos[1] + ah * (leverpos[1] - rootpos[1]) / dism,0 };
			Vector3d PP1 = { PP0[0] + ha * (leverpos[1] - rootpos[1]) / dism,PP0[1] - ha * (leverpos[0] - rootpos[0]) / dism ,0 };
			Vector3d PP2 = { PP0[0] - ha * (leverpos[1] - rootpos[1]) / dism,PP0[1] + ha * (leverpos[0] - rootpos[0]) / dism ,0 };
			int pts;
			if (scheme == 0) {
				pntpos = { 2 * am + 0 ,2 * am + 1,2 * bm + 0,2 * bm + 1,2 * cm + 0,2 * cm + 1 };
				nodescur = { A,B,C };
				pts = 3;
				nodesref = { _MPos3D[a % _npart],_MPos3D[b % _npart],_MPos3D[c % _npart] };
				_quadpoints.resize(1, { 1.0 / 3.0, 1.0 / 3.0 });
				_quadweights.resize(1, 0.5);
			}
			else if (scheme == 1) {
				pntpos = { 2 * am + 0 ,2 * am + 1,2 * bm + 0,2 * bm + 1,2 * cm + 0,2 * cm + 1 ,2*pms[0]+0,2 * pms[0] + 1,2 * pms[1] + 0,2 * pms[1] + 1,2 * pms[2] + 0,2 * pms[2] + 1 };
				nodescur = { A,B,C,D,E,FF };
				pts = 6;
				if(refdome)nodesref = { _MPos3D[a % _npart],_MPos3D[b % _npart],_MPos3D[c % _npart],finalmids[0],finalmids[1],finalmids[2] };	// For domes
				else nodesref = { _MPos3D[a % _npart],_MPos3D[b % _npart],_MPos3D[c % _npart],NP1,NP2,NP3 };			// For parabola, icosahedron
				_quadpoints.resize(3, { 0.0,0.0 }); 
				_quadpoints[0] = { 1.0 / 6.0, 2.0 / 3.0 };
				_quadpoints[1] = { 1.0 / 6.0, 1.0 / 6.0 };
				_quadpoints[2] = { 2.0 / 3.0, 1.0 / 6.0 };
				_quadweights.resize(3, 1.0/6.0);
			}
			vector<Vector3d> basis, dual, refbasis, refdual;
			basis.resize(2, { 0.0,0.0,0.0 });
			dual.resize(2, { 0.0,0.0,0.0 });
			refbasis.resize(2, { 0.0,0.0,0.0 });
			refdual.resize(2, { 0.0,0.0,0.0 });
			MatrixXd F(3, 3), P(3, 3), Placeholder(3, 3), Ucomm(3,3), Ufin(3,3);
			// Iterating through quadrature points to update energy and gradients
			for (int qds = 0; qds < _quadpoints.size(); qds++) {
				vector<Vector3d>_n;
				_n.resize(3, { 0.0,0.0,0.0 });
				if (scheme == 0) {
					N = NTri3(_quadpoints[qds]);
					DN = DNTri3(_quadpoints[qds]);
				}
				if (scheme == 1) {
					N = NTri6(_quadpoints[qds]);
					DN = DNTri6(_quadpoints[qds]);
				}
				MatrixXd Iso(3,3),EPS(3, 3);

				// Reference Configuration Basis and Duals
				refbasis = basiss(nodesref, DN);
				MatrixXd metrictensorref(2, 2);
				metrictensorref = metric(refbasis);
				Real G = sqrtdetmetric(metrictensorref);
				MatrixXd metrictensorinvref(2, 2);
				metrictensorinvref = metricinv(metrictensorref);
				refdual = duals(refbasis, metrictensorinvref);
				

				// Deformed Configuration Basis and Duals
				basis = basiss(nodescur, DN);
				MatrixXd metrictensor(2, 2); 
				metrictensor = metric(basis);
				Real g = sqrtdetmetric(metrictensor);
				MatrixXd metrictensorinv(2, 2);
				metrictensorinv = metricinv(metrictensor);
				dual = duals(basis, metrictensorinv);
				Real weight = _quadweights[qds] * G;
				for (int I = 0; I < 3; I++) {
					for (int J = 0; J < 3; J++) {
						F(I, J) = 0.0;
						P(I, J) = 0.0;
						EPS(I, J) = 0.0;
						Ucomm(I, J) = 0.0;
						Iso(I, J) = 0.0;
						Ufin(I, J) = 0, 0;
						Placeholder(I, J) = 0.0;
					}
				}
				for (int alpha = 0; alpha < 2; alpha++) {
					for (int I = 0; I < 3; I++) {
						for (int J = 0; J < 3; J++) {
							F(I, J) += basis[alpha][I] * refdual[alpha][J];
						}
					}
				}

				pair<MatrixXd, Real> Res = EnergyEvans(F, "Both", _mu, _kS,0);
				//pair<MatrixXd, Real> Res = EnergyHooke(F, "Both", muval, lambdaval, 0);
				//pair<MatrixXd, Real> Res = EnergyHookeEPS(EPS, "Both", muval, lambdaval, 0);
				
				P = Res.first;
				_Energy += weight*Res.second;

				for (int alpha = 0; alpha < 2; alpha++)
					_n[alpha] += P * refdual[alpha];
				vector<Real>Derivatives;
				Derivatives.resize(2 * nodescur.size(), 0.0);
				
				for (int nb = 0; nb < nodescur.size(); nb++) {
					for (int alpha = 0; alpha < 2; alpha++) {
						Vector3d ftmp;
						ftmp = { 0.0,0.0,0.0 };
						
						ftmp = _n[alpha] * DN[nb](alpha) * weight;
						Derivatives[2*nb] += ftmp[0];
						Derivatives[2*nb+1] += ftmp[1];
					}
				}
				for (int nb = 0; nb < nodescur.size(); nb++) {
					_GrdOrdrd[pntpos[2 * nb]] += Derivatives[2 * nb];
					_GrdOrdrd[pntpos[2 * nb + 1]] += Derivatives[2 * nb + 1];						
				}
				
				Real h = 1e-8;
				Real tol = 1e-5;
				// This section checks element level consistency
				if (checkelement) {
					for (int nb = 0; nb < nodescur.size(); nb++) {
						for (int ps = 0; ps < 2; ps++) {
							vector<Vector3d> basisP, dualP, basisM, dualM;
							basisP.resize(2, { 0.0,0.0,0.0 });
							dualP.resize(2, { 0.0,0.0,0.0 });
							basisM.resize(2, { 0.0,0.0,0.0 });
							dualM.resize(2, { 0.0,0.0,0.0 });
							vector<Vector3d>newnodes = nodescur;
							newnodes[nb][ps] += h;
							basisP = basiss(newnodes, DN);
							MatrixXd metrictensorP(2, 2);
							metrictensorP = metric(basisP);
							Real gP = sqrtdetmetric(metrictensorP);
							MatrixXd metrictensorinvP(2, 2);
							metrictensorinvP = metricinv(metrictensorP);
							dualP = duals(basisP, metrictensorinvP);
							newnodes[nb][ps] -= 2.0*h;
							basisM = basiss(newnodes, DN);
							MatrixXd metrictensorM(2, 2);
							metrictensorM = metric(basisM);
							Real gM = sqrtdetmetric(metrictensorM);
							MatrixXd metrictensorinvM(2, 2);
							metrictensorinvM = metricinv(metrictensorM);
							dualM = duals(basisM, metrictensorinvM);
							MatrixXd FP(3, 3), FM(3, 3), PP(3,3),PM(3,3);
							FP = Placeholder;
							FM = Placeholder;
							PP = Placeholder;
							PM = Placeholder;
							for (int alpha = 0; alpha < 2; alpha++) {
								for (int I = 0; I < 3; I++) {
									for (int J = 0; J < 3; J++) {
										FP(I, J) += basisP[alpha][I] * refdual[alpha][J];
										FM(I, J) += basisM[alpha][I] * refdual[alpha][J];
									}
								}
							}
							//pair<MatrixXd, Real> ResP = EnergyHooke(FP, "Both", muval, lambdaval, 0);
							//pair<MatrixXd, Real> ResM = EnergyHooke(FM, "Both", muval, lambdaval, 0);
							pair<MatrixXd, Real> ResP = EnergyEvans(FP, "Both", _mu, _kS,0);
							pair<MatrixXd, Real> ResM = EnergyEvans(FM, "Both", _mu, _kS,0);
							Real EP = weight * ResP.second;
							Real EM = weight * ResM.second;
							Real errE = (EP - EM) / (2.0 * h) - Derivatives[2 * nb+ps];
							if (abs(errE) > tol) {
								cout << "problem Element " << setprecision(10) << abs(errE) << endl;
								cout << a << " " << b << " " << c << endl;
								cout << A << endl;
								cout << B << endl;
								cout << C << endl;
								cout << D << endl;
								cout << E << endl;
								cout << FF << endl;
								cout << _MPos3D[a % _npart] << endl;
								cout << _MPos3D[b % _npart] << endl;
								cout << _MPos3D[c % _npart] << endl;
								cout << NP1 << endl;
								cout << NP2 << endl;
								cout << NP3 << endl;
							}
						}
					}
				}
			}		
		}
		// This section checks body level consistency
		if (checkbody) {
			for (int N = 0; N < _npart; N++) {
				Real hb = 1e-7;
				Real tolb = 1e-5;
				for (int K = 0; K < 2; K++) {
					Real EnergyP = 0.0;
					Real EnergyM = 0.0;
					_PosOrdrd[2 * N + K] += hb;
					for (int trr = 0; trr < _Tris.size(); trr++) {
						int a = _Tris[trr][0];
						int b = _Tris[trr][1];
						int c = _Tris[trr][2];
						vector<int>Ts;
						Ts = { a % _npart,b % _npart,c % _npart };
						Vector3d D, E, FF;
						vector<pair<int, int>>ptpairs;
						pair<int, int>pr1, pr2, pr3;
						if (a < b)pr1 = make_pair(a, b);
						else pr1 = make_pair(b, a);
						if (b < c)pr2 = make_pair(b, c);
						else pr2 = make_pair(c, b);
						if (a < c)pr3 = make_pair(a, c);
						else pr3 = make_pair(c, a);
						ptpairs.push_back(pr1);
						ptpairs.push_back(pr2);
						ptpairs.push_back(pr3);
						vector<int>pms;
						for (int ptk = 0; ptk < 3; ptk++) {
							int pt1 = ptpairs[ptk].first;
							int pt2 = ptpairs[ptk].second;
							string ptnow = to_string(pt1) + "a" + to_string(pt2);
							pms.push_back(_Midpts[ptnow]);
						}
						int am = _ord3D[a];
						int bm = _ord3D[b];
						int cm = _ord3D[c];
						vector<int>pntpos;

						vector<Vector2d>_quadpoints;
						vector<Real>_quadweights;
						vector<Vector3d>nodescur;
						vector<Vector3d>nodesref;
						vector<Real>N;
						vector<Vector2d>DN;
						Vector3d A = { _PosOrdrd[2 * am + 0],_PosOrdrd[2 * am + 1],0.0 };
						Vector3d B = { _PosOrdrd[2 * bm + 0],_PosOrdrd[2 * bm + 1],0.0 };
						Vector3d C = { _PosOrdrd[2 * cm + 0],_PosOrdrd[2 * cm + 1],0.0 };
						D = { _PosOrdrd[2 * pms[0] + 0],_PosOrdrd[2 * pms[0] + 1],0.0 };
						E = { _PosOrdrd[2 * pms[1] + 0],_PosOrdrd[2 * pms[1] + 1],0.0 };
						FF = { _PosOrdrd[2 * pms[2] + 0],_PosOrdrd[2 * pms[2] + 1],0.0 };

						vector<Vector3d>NPs;
						Vector3d NP1, NP2, NP3;
						NP1 = (_MPos3D[a % _npart] + _MPos3D[b % _npart]) / 2.0;
						NP2 = (_MPos3D[b % _npart] + _MPos3D[c % _npart]) / 2.0;
						NP3 = (_MPos3D[c % _npart] + _MPos3D[a % _npart]) / 2.0;
						NPs = { NP1,NP2,NP3 };
						bool parabola = 0;
						if (parabola) {
							NP1[0] = (NP1[1] * NP1[1] + NP1[2] * NP1[2]) / 8.0;
							NP2[0] = (NP2[1] * NP2[1] + NP2[2] * NP2[2]) / 8.0;
							NP3[0] = (NP3[1] * NP3[1] + NP3[2] * NP3[2]) / 8.0;
						}
						vector<Vector3d>finalmids;
						for (int mids = 0; mids < 3; mids++) {
							Vector3d NPnow = NPs[mids];
							//cout << NPnow << endl;
							Real Radius = sqrt(NPnow[0] * NPnow[0] + NPnow[1] * NPnow[1] + NPnow[2] * NPnow[2]);
							Real theta = acos(NPnow[2] / Radius);
							Real phi = atan2(NPnow[1], NPnow[0]);
							//cout << Radius << " " << theta << " " << " " << phi << endl;
							finalmids.push_back({ 2.0 * sin(theta) * cos(phi),2.0 * sin(theta) * sin(phi) ,2.0 * cos(theta) });
							//cout << 2.0 * sin(theta) * cos(phi) << " " << 2.0 * sin(theta) * sin(phi) << " " << 2.0 * cos(theta) << endl;
						}
						Real Aref = 0.0;
						if (scheme == 0) {
							pntpos = { 2 * am + 0 ,2 * am + 1,2 * bm + 0,2 * bm + 1,2 * cm + 0,2 * cm + 1 };
							nodescur = { A,B,C };
							nodesref = { _MPos3D[a % _npart],_MPos3D[b % _npart],_MPos3D[c % _npart] };
							_quadpoints.resize(1, { 1.0 / 3.0, 1.0 / 3.0 });
							_quadweights.resize(1, 0.5);
						}
						else if (scheme == 1) {
							pntpos = { 2 * am + 0 ,2 * am + 1,2 * bm + 0,2 * bm + 1,2 * cm + 0,2 * cm + 1 ,2 * pms[0] + 0,2 * pms[0] + 1,2 * pms[1] + 0,2 * pms[1] + 1,2 * pms[2] + 0,2 * pms[2] + 1 };
							nodescur = { A,B,C,D,E,FF };
							if(refdome)nodesref = { _MPos3D[a % _npart],_MPos3D[b % _npart],_MPos3D[c % _npart],finalmids[0],finalmids[1],finalmids[2] };
							else nodesref = { _MPos3D[a % _npart],_MPos3D[b % _npart],_MPos3D[c % _npart],NP1,NP2,NP3 };
							_quadpoints.resize(3, { 0.0,0.0 });
							_quadpoints[0] = { 1.0 / 6.0, 2.0 / 3.0 };
							_quadpoints[1] = { 1.0 / 6.0, 1.0 / 6.0 };
							_quadpoints[2] = { 2.0 / 3.0, 1.0 / 6.0 };
							//_quadpoints[0] = { 0.5, 0.5 };
							//_quadpoints[1] = { 0.0, 0.5 };
							//_quadpoints[2] = { 0.5, 0.0 };
							_quadweights.resize(3, 1.0 / 6.0);
						}
						MatrixXd F(3, 3);
						vector<Vector3d> basis, dual, refbasis, refdual;
						basis.resize(2, { 0.0,0.0,0.0 });
						dual.resize(2, { 0.0,0.0,0.0 });
						refbasis.resize(2, { 0.0,0.0,0.0 });
						refdual.resize(2, { 0.0,0.0,0.0 });
						Aref = 0.0;
						Real Adef = 0.0;
						Real _Jtotal = 0.0;
						Real TRCtotal = 0.0;
						for (int qds = 0; qds < _quadpoints.size(); qds++) {
							vector<Vector3d>_n;
							_n.resize(3, { 0.0,0.0,0.0 });
							if (scheme == 0) {
								N = NTri3(_quadpoints[qds]);
								DN = DNTri3(_quadpoints[qds]);
							}
							if (scheme == 1) {
								N = NTri6(_quadpoints[qds]);
								DN = DNTri6(_quadpoints[qds]);
							}
							// Reference Configuration Basis and Duals
							refbasis = basiss(nodesref, DN);
							MatrixXd metrictensorref(2, 2);
							metrictensorref = metric(refbasis);
							Real G = sqrtdetmetric(metrictensorref);
							MatrixXd metrictensorinvref(2, 2);
							metrictensorinvref = metricinv(metrictensorref);
							refdual = duals(refbasis, metrictensorinvref);


							// Deformed Configuration Basis and Duals
							basis = basiss(nodescur, DN);
							MatrixXd metrictensor(2, 2);
							metrictensor = metric(basis);
							Real g = sqrtdetmetric(metrictensor);
							MatrixXd metrictensorinv(2, 2);
							metrictensorinv = metricinv(metrictensor);
							dual = duals(basis, metrictensorinv);
							for (int I = 0; I < 3; I++) {
								for (int J = 0; J < 3; J++) {
									F(I, J) = 0.0;
								}
							}
							for (int alpha = 0; alpha < 2; alpha++) {
								for (int I = 0; I < 3; I++) {
									for (int J = 0; J < 3; J++) {
										F(I, J) += basis[alpha][I] * refdual[alpha][J];
									}
								}
							}
							Real weight = _quadweights[qds] * G;
							//pair<MatrixXd, Real> ResE = EnergyEvans(F, "Both", 1.0, 1.0, 0);
							pair<Real, Real> ResE = EnergyEvansEnergies(F, "Both", 1.0, 1.0);
							//Energy1 += weight * (ResE.first);
							//Energy2 += weight * (ResE.second);
							//pair<MatrixXd, Real> ResE = EnergyHooke(F, "Both", muval, lambdaval, 0);
							//EnergyTotal = Energy1 + Energy2;
							EnergyP += weight * ResE.second;
						}
					}
					// This section applies body force on predefined nodes of the model
					for (int i = 0; i < _npart; i++) {
						if (_MPos3D[i][0] == 4.0) {
							Real factorb = pow(0.5, abs(_MPos3D[i][1] - 2.0) / L);
							Real Finalforceb = Force * factorb;
							Vector2d A = { _MPos3D[i][0],_MPos3D[i][1] };
							Vector2d B = { _PosOrdrd[2 * _ord3D[i]],_PosOrdrd[2 * _ord3D[i] + 1] };
							//EnergyP += -Finalforceb * L * (_PosOrdrd[2 * _ord3D[i]] - _MPos3D[i][0]);
						}
					}
					_PosOrdrd[2 * N + K] -= 2.0*hb;
					for (int trr = 0; trr < _Tris.size(); trr++) {
						int a = _Tris[trr][0];
						int b = _Tris[trr][1];
						int c = _Tris[trr][2];
						vector<int>Ts;
						Ts = { a % _npart,b % _npart,c % _npart };
						Vector3d D, E, FF;
						vector<pair<int, int>>ptpairs;
						pair<int, int>pr1, pr2, pr3;
						if (a < b)pr1 = make_pair(a, b);
						else pr1 = make_pair(b, a);
						if (b < c)pr2 = make_pair(b, c);
						else pr2 = make_pair(c, b);
						if (a < c)pr3 = make_pair(a, c);
						else pr3 = make_pair(c, a);
						ptpairs.push_back(pr1);
						ptpairs.push_back(pr2);
						ptpairs.push_back(pr3);
						vector<int>pms;
						for (int ptk = 0; ptk < 3; ptk++) {
							int pt1 = ptpairs[ptk].first;
							int pt2 = ptpairs[ptk].second;
							string ptnow = to_string(pt1) + "a" + to_string(pt2);
							pms.push_back(_Midpts[ptnow]);
						}
						int am = _ord3D[a];
						int bm = _ord3D[b];
						int cm = _ord3D[c];
						vector<int>pntpos;

						vector<Vector2d>_quadpoints;
						vector<Real>_quadweights;
						vector<Vector3d>nodescur;
						vector<Vector3d>nodesref;
						vector<Real>N;
						vector<Vector2d>DN;
						Vector3d A = { _PosOrdrd[2 * am + 0],_PosOrdrd[2 * am + 1],0.0 };
						Vector3d B = { _PosOrdrd[2 * bm + 0],_PosOrdrd[2 * bm + 1],0.0 };
						Vector3d C = { _PosOrdrd[2 * cm + 0],_PosOrdrd[2 * cm + 1],0.0 };
						D = { _PosOrdrd[2 * pms[0] + 0],_PosOrdrd[2 * pms[0] + 1],0.0 };
						E = { _PosOrdrd[2 * pms[1] + 0],_PosOrdrd[2 * pms[1] + 1],0.0 };
						FF = { _PosOrdrd[2 * pms[2] + 0],_PosOrdrd[2 * pms[2] + 1],0.0 };

						vector<Vector3d>NPs;
						Vector3d NP1, NP2, NP3;
						NP1 = (_MPos3D[a % _npart] + _MPos3D[b % _npart]) / 2.0;
						NP2 = (_MPos3D[b % _npart] + _MPos3D[c % _npart]) / 2.0;
						NP3 = (_MPos3D[c % _npart] + _MPos3D[a % _npart]) / 2.0;
						NPs = { NP1,NP2,NP3 };
						// This section calculates projection of midpoints to parabola surface
						bool parabola = 0;
						if (parabola) {
							NP1[0] = (NP1[1] * NP1[1] + NP1[2] * NP1[2]) / 8.0;
							NP2[0] = (NP2[1] * NP2[1] + NP2[2] * NP2[2]) / 8.0;
							NP3[0] = (NP3[1] * NP3[1] + NP3[2] * NP3[2]) / 8.0;
						}
						// This section calculates projection of midpoints to dome surface
						vector<Vector3d>finalmids;
						for (int mids = 0; mids < 3; mids++) {
							Vector3d NPnow = NPs[mids];
							//cout << NPnow << endl;
							Real Radius = sqrt(NPnow[0] * NPnow[0] + NPnow[1] * NPnow[1] + NPnow[2] * NPnow[2]);
							Real theta = acos(NPnow[2] / Radius);
							Real phi = atan2(NPnow[1], NPnow[0]);
							//cout << Radius << " " << theta << " " << " " << phi << endl;
							finalmids.push_back({ 2.0 * sin(theta) * cos(phi),2.0 * sin(theta) * sin(phi) ,2.0 * cos(theta) });
							//cout << 2.0 * sin(theta) * cos(phi) << " " << 2.0 * sin(theta) * sin(phi) << " " << 2.0 * cos(theta) << endl;
						}
						Real Aref = 0.0;
						if (scheme == 0) {
							pntpos = { 2 * am + 0 ,2 * am + 1,2 * bm + 0,2 * bm + 1,2 * cm + 0,2 * cm + 1 };
							nodescur = { A,B,C };
							nodesref = { _MPos3D[a % _npart],_MPos3D[b % _npart],_MPos3D[c % _npart] };
							_quadpoints.resize(1, { 1.0 / 3.0, 1.0 / 3.0 });
							_quadweights.resize(1, 0.5);
						}
						else if (scheme == 1) {
							pntpos = { 2 * am + 0 ,2 * am + 1,2 * bm + 0,2 * bm + 1,2 * cm + 0,2 * cm + 1 ,2 * pms[0] + 0,2 * pms[0] + 1,2 * pms[1] + 0,2 * pms[1] + 1,2 * pms[2] + 0,2 * pms[2] + 1 };
							nodescur = { A,B,C,D,E,FF };
							if(refdome)nodesref = { _MPos3D[a % _npart],_MPos3D[b % _npart],_MPos3D[c % _npart],finalmids[0],finalmids[1],finalmids[2] };
							else nodesref = { _MPos3D[a % _npart],_MPos3D[b % _npart],_MPos3D[c % _npart],NP1,NP2,NP3 };
							_quadpoints.resize(3, { 0.0,0.0 });
							_quadpoints[0] = { 1.0 / 6.0, 2.0 / 3.0 };
							_quadpoints[1] = { 1.0 / 6.0, 1.0 / 6.0 };
							_quadpoints[2] = { 2.0 / 3.0, 1.0 / 6.0 };
							//_quadpoints[0] = { 0.5, 0.5 };
							//_quadpoints[1] = { 0.0, 0.5 };
							//_quadpoints[2] = { 0.5, 0.0 };
							_quadweights.resize(3, 1.0 / 6.0);
						}
						MatrixXd F(3, 3);
						vector<Vector3d> basis, dual, refbasis, refdual;
						basis.resize(2, { 0.0,0.0,0.0 });
						dual.resize(2, { 0.0,0.0,0.0 });
						refbasis.resize(2, { 0.0,0.0,0.0 });
						refdual.resize(2, { 0.0,0.0,0.0 });
						Aref = 0.0;
						Real Adef = 0.0;
						Real _Jtotal = 0.0;
						Real TRCtotal = 0.0;
						for (int qds = 0; qds < _quadpoints.size(); qds++) {
							vector<Vector3d>_n;
							_n.resize(3, { 0.0,0.0,0.0 });
							if (scheme == 0) {
								N = NTri3(_quadpoints[qds]);
								DN = DNTri3(_quadpoints[qds]);
							}
							if (scheme == 1) {
								N = NTri6(_quadpoints[qds]);
								DN = DNTri6(_quadpoints[qds]);
							}
							// Reference Configuration Basis and Duals
							refbasis = basiss(nodesref, DN);
							MatrixXd metrictensorref(2, 2);
							metrictensorref = metric(refbasis);
							Real G = sqrtdetmetric(metrictensorref);
							MatrixXd metrictensorinvref(2, 2);
							metrictensorinvref = metricinv(metrictensorref);
							refdual = duals(refbasis, metrictensorinvref);


							// Deformed Configuration Basis and Duals
							basis = basiss(nodescur, DN);
							MatrixXd metrictensor(2, 2);
							metrictensor = metric(basis);
							Real g = sqrtdetmetric(metrictensor);
							MatrixXd metrictensorinv(2, 2);
							metrictensorinv = metricinv(metrictensor);
							dual = duals(basis, metrictensorinv);
							for (int I = 0; I < 3; I++) {
								for (int J = 0; J < 3; J++) {
									F(I, J) = 0.0;
								}
							}
							for (int alpha = 0; alpha < 2; alpha++) {
								for (int I = 0; I < 3; I++) {
									for (int J = 0; J < 3; J++) {
										F(I, J) += basis[alpha][I] * refdual[alpha][J];
									}
								}
							}
							Real weight = _quadweights[qds] * G;
							//pair<MatrixXd, Real> ResE = EnergyEvans(F, "Both", 1.0, 1.0, 0);
							pair<Real, Real> ResE = EnergyEvansEnergies(F, "Both", 1.0, 1.0);
							//Energy1 += weight * (ResE.first);
							//Energy2 += weight * (ResE.second);
							//pair<MatrixXd, Real> ResE = EnergyHooke(F, "Both", muval, lambdaval, 0);
							//EnergyTotal = Energy1 + Energy2;
							EnergyM += weight * ResE.second;
						}
					}
					// This section applies body force on predefined nodes of the model
					for (int i = 0; i < _npart; i++) {
						if (_MPos3D[i][0] == 4.0) {
							Real factorm = pow(0.5, abs(_MPos3D[i][1] - 2.0) / L);
							Real Finalforcem = Force * factorm;
							Vector2d A = { _MPos3D[i][0],_MPos3D[i][1] };
							Vector2d B = { _PosOrdrd[2 * _ord3D[i]],_PosOrdrd[2 * _ord3D[i] + 1] };
							//EnergyM += -Finalforcem * L * (_PosOrdrd[2 * _ord3D[i]] - _MPos3D[i][0]);
						}
					}
					Real errBody = (EnergyP - EnergyM) / (2.0 * hb) - _GrdOrdrd[2 * N + K];
					if (abs(errBody) > tolb) {
						cout << "Body Level Consistency Failed" << endl;;
					}
					_PosOrdrd[2 * N + K] += hb;
				}
			}
		}
		return;
	}

	void Body::checkconsistencyEdgTri(Real h, Real tol) {
		vector<Real>Err;
		for (int i = 0; i < _Edges.size(); i++) {
			int a = _Edges[i].first;
			int b = _Edges[i].second;
			int c = _EdgeTri[i];
			Real l = dist(_MPos3D[a % _npart], _MPos3D[b % _npart]);
			Real Al = trianglearea(_MPos3D[a % _npart], _MPos3D[b % _npart], _MPos3D[c % _npart]);
			int am = _ord3D[a];
			int bm = _ord3D[b];
			Vector2d A = { _PosOrdrd[2 * am + 0],_PosOrdrd[2 * am + 1] };
			Vector2d B = { _PosOrdrd[2 * bm + 0],_PosOrdrd[2 * bm + 1] };
			Real d = dist(A, B);
			Real enrg = d - l;
			_Energy += 0.5 * Al * pow(enrg, 2.0);
			_GrdOrdrd[2 * am + 0] += Al * (enrg) * (A[0] - B[0]) / d;
			_GrdOrdrd[2 * am + 1] += Al * (enrg) * (A[1] - B[1]) / d;
			_GrdOrdrd[2 * bm + 0] += Al * (enrg) * (B[0] - A[0]) / d;
			_GrdOrdrd[2 * bm + 1] += Al * (enrg) * (B[1] - A[1]) / d;
			Vector2d Xplus = { _PosOrdrd[2 * am + 0] + h,_PosOrdrd[2 * am + 1] };
			Vector2d Xminus = { _PosOrdrd[2 * am + 0] - h,_PosOrdrd[2 * am + 1] };
			Real dplus = dist(Xplus, B);
			Real dminus = dist(Xminus, B);
			Real Eplus = 0.5 * pow(dplus - l, 2.0);
			Real Eminus = 0.5 * pow(dminus - l, 2.0);
			Real error = ((Eplus - Eminus) / (2.0 * h)) - (enrg) * (A[0] - B[0]) / d;
			Err.push_back(error);
			if ((abs(error)) > tol)cout << "problem" << endl;
		}
		Real Toter = 0;
		Real maxFault = -10;
		for (int i = 0; i < Err.size(); i++) {
			Toter = Err[i] * Err[i];
		}
		Toter = sqrt(Toter);
	}
	void Body::checkconsistencyNeoHookean2D(Real h, Real tol) {
		for (int i = 0; i < _Tris.size(); i++) {
			int a = _Tris[i][0];
			int b = _Tris[i][1];
			int c = _Tris[i][2];
			int am = _ord3D[a];
			int bm = _ord3D[b];
			int cm = _ord3D[c];
			Real Al = abs(trianglearea(_MPos3D[a % _npart], _MPos3D[b % _npart], _MPos3D[c % _npart]));
			Real b1, c2, c1;
			Vector3d A = { _PosOrdrd[2 * am + 0],_PosOrdrd[2 * am + 1],0.0 };
			Vector3d B = { _PosOrdrd[2 * bm + 0],_PosOrdrd[2 * bm + 1],0.0 };
			Vector3d C = { _PosOrdrd[2 * cm + 0],_PosOrdrd[2 * cm + 1],0.0 };
			Vector3d AB = B - A;
			Vector3d AC = C - A;
			MatrixXd F(2, 2);
			Real N1x, N1y, N2x, N2y, N3x, N3y;
			bool change = 0;
			Real fact = 0.0;
			if (AB.cross(AC)[2] < 0) {
				fact = -1.0;
			}
			vector<Vector3d>ABC = { A,B,C };
			Real dism = dist(_MPos3D[a % _npart], _MPos3D[b % _npart]);
			Real dis1 = dist(_MPos3D[a % _npart], _MPos3D[c % _npart]);
			Real dis2 = dist(_MPos3D[b % _npart], _MPos3D[c % _npart]);
			Real ah = (dis1 * dis1 - dis2 * dis2 + dism * dism) / (2 * dism);
			Real ha = sqrt(dis1 * dis1 - ah * ah);
			Vector3d rootpos = { 0.0,0.0,0.0 };
			Vector3d leverpos = { dism,0.0,0.0 };
			Vector3d PP0 = { rootpos[0] + ah * (leverpos[0] - rootpos[0]) / dism,rootpos[1] + ah * (leverpos[1] - rootpos[1]) / dism,0 };
			Vector3d PP1 = { PP0[0] + ha * (leverpos[1] - rootpos[1]) / dism,PP0[1] - ha * (leverpos[0] - rootpos[0]) / dism ,0 };
			Vector3d PP2 = { PP0[0] - ha * (leverpos[1] - rootpos[1]) / dism,PP0[1] + ha * (leverpos[0] - rootpos[0]) / dism ,0 };
			if (fact == 1.0) {
				if (PP1[1] > 0.0) {
					c1 = PP1[0];
					c2 = PP1[1];
				}
				else {
					c1 = PP2[0];
					c2 = PP2[1];
				}
			}
			else {
				if (PP1[1] < 0.0) {
					c1 = PP1[0];
					c2 = PP1[1];
				}
				else {
					c1 = PP2[0];
					c2 = PP2[1];
				}
			}
			b1 = dism;
			N1x = -(1.0 / b1);
			N1y = c1 / (b1 * c2) - (1 / c2);
			N2x = -N1x;
			N2y = -(c1 / (b1 * c2));
			N3x = 0.0;
			N3y = 1.0 / c2;
			F(0, 0) = N1x * A[0] + N2x * B[0];
			F(0, 1) = N1y * A[0] + N2y * B[0] + N3y * C[0];
			F(1, 0) = N1x * A[1] + N2x * B[1];
			F(1, 1) = N1y * A[1] + N2y * B[1] + N3y * C[1];
			_Energy += FEnergy2D(F) * Al;
			MatrixXd P = Pmatrix2D(F);
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					MatrixXd Fp = F;
					MatrixXd Fm = F;
					Fp(i, j) += h;
					Fm(i, j) -= h;
					Real Ep = FEnergy2D(Fp);
					Real Em = FEnergy2D(Fm);
					Real errM = (Ep - Em) / (2.0 * h) - P(i, j);
					if (abs(errM) > tol)cout << "problem Material" << setprecision(10) << abs(errM) << endl;
				}
			}
			vector<vector<Real>>Derivatives;
			Derivatives.push_back({ 0.0,0.0 });
			Derivatives.push_back({ 0.0,0.0 });
			Derivatives.push_back({ 0.0,0.0 });
			Derivatives[0][0] = Al * (P(0, 0) * N1x + P(0, 1) * N1y);
			Derivatives[0][1] = Al * (P(1, 0) * N1x + P(1, 1) * N1y);
			Derivatives[1][0] = Al * (P(0, 0) * N2x + P(0, 1) * N2y);
			Derivatives[1][1] = Al * (P(1, 0) * N2x + P(1, 1) * N2y);
			Derivatives[2][0] = Al * (P(0, 0) * N3x + P(0, 1) * N3y);
			Derivatives[2][1] = Al * (P(1, 0) * N3x + P(1, 1) * N3y);
			for (int i = 0; i < 3; i++) {
				Vector3d Now = ABC[i];
				for (int j = 0; j < 2; j++) {
					ABC[i][j] += h;
					MatrixXd Fplus(2, 2);
					Fplus(0, 0) = N1x * ABC[0][0] + N2x * ABC[1][0];
					Fplus(0, 1) = N1y * ABC[0][0] + N2y * ABC[1][0] + N3y * ABC[2][0];
					Fplus(1, 0) = N1x * ABC[0][1] + N2x * ABC[1][1];
					Fplus(1, 1) = N1y * ABC[0][1] + N2y * ABC[1][1] + N3y * ABC[2][1];
					ABC[i][j] -= 2.0 * h;
					MatrixXd Fminus(2, 2);
					Fminus(0, 0) = N1x * ABC[0][0] + N2x * ABC[1][0];
					Fminus(0, 1) = N1y * ABC[0][0] + N2y * ABC[1][0] + N3y * ABC[2][0];
					Fminus(1, 0) = N1x * ABC[0][1] + N2x * ABC[1][1];
					Fminus(1, 1) = N1y * ABC[0][1] + N2y * ABC[1][1] + N3y * ABC[2][1];
					Real Eplus = Al * (FEnergy2D(Fplus));
					Real Eminus = Al * (FEnergy2D(Fminus));
					ABC[i] = Now;
					Real errE = (Eplus - Eminus) / (2.0 * h) - Derivatives[i][j];
					if (abs(errE) > tol)cout << "problem Element" << setprecision(10) << abs(errE) << endl;
				}
			}
			if (change == 1) {
				_GrdOrdrd[2 * am + 0] += Derivatives[0][0];
				_GrdOrdrd[2 * am + 1] += Derivatives[0][1];
				_GrdOrdrd[2 * bm + 0] += Derivatives[2][0];
				_GrdOrdrd[2 * bm + 1] += Derivatives[2][1];
				_GrdOrdrd[2 * cm + 0] += Derivatives[1][0];
				_GrdOrdrd[2 * cm + 1] += Derivatives[1][1];
			}
			else {
				_GrdOrdrd[2 * am + 0] += Derivatives[0][0];
				_GrdOrdrd[2 * am + 1] += Derivatives[0][1];
				_GrdOrdrd[2 * bm + 0] += Derivatives[1][0];
				_GrdOrdrd[2 * bm + 1] += Derivatives[1][1];
				_GrdOrdrd[2 * cm + 0] += Derivatives[2][0];
				_GrdOrdrd[2 * cm + 1] += Derivatives[2][1];
			}
		}
	}

	vector<Real> Body::checkDiffEdg(int wt) {
		vector<Real>difs;
		int countr = 0;
		int countrb = 0;
		for (int i = 0; i < _Edges.size(); i++) {
			int a = _Edges[i].first;
			int b = _Edges[i].second;
			Real l = dist(_MPos3D[a % _npart], _MPos3D[b % _npart]);
			int am = _ord3D[a];
			int bm = _ord3D[b];

			Vector2d A = { _PosOrdrd[2 * am + 0],_PosOrdrd[2 * am + 1] };
			Vector2d B = { _PosOrdrd[2 * bm + 0],_PosOrdrd[2 * bm + 1] };

			Real d = dist(A, B);

			Real enrg = (pow(d, 2.0) - pow(l, 2.0));
			Real pd = 100 * abs(l - d) / l;
			if (wt == 0) {
				difs.push_back(pd);
			}
			else if (wt == 1) {
				if (bndred[i] == 1)difs.push_back(pd);
			}
			else {
				if (bndred[i] == 0)difs.push_back(pd);
			}
		}
		cout << "Bad Edges: Inner :  " << countr << endl;
		cout << "Bad Edges: Boundary :  " << countrb << endl;
		return difs;
	}
	void Body::checkDiffTri(int scheme) {
		vector< pair<vector<int>, Real>> Diffs;
		vector< pair<vector<int>, Real>> Js;
		vector< pair<vector<int>, Real>> Trcs;
		int countr = 0;
		Real TotalActual=0.0;
		Real TotalDeformed=0.0;
		Real TotalDistortion = 0.0;
		Real Energy1 = 0.0;
		Real Energy2 = 0.0;
		Real EnergyTotal = 0.0;
		Real FinalJ = 0.0;
		Real FinalTrc = 0.0;
		Real pts = sqrt(_npart);
		if (scheme)pts += (pts - 1);
		Real L = 4.0 / (pts - 1.0);
		set<int>doneforcing;
		for (int trr = 0; trr < _Tris.size(); trr++) {
			int a = _Tris[trr][0];
			int b = _Tris[trr][1];
			int c = _Tris[trr][2];
			vector<int>Ts;
			Ts = { a % _npart,b % _npart,c % _npart };
			Vector3d D, E, FF;
			vector<pair<int, int>>ptpairs;
			pair<int, int>pr1, pr2, pr3;
			if (a < b)pr1 = make_pair(a, b);
			else pr1 = make_pair(b, a);
			if (b < c)pr2 = make_pair(b, c);
			else pr2 = make_pair(c, b);
			if (a < c)pr3 = make_pair(a, c);
			else pr3 = make_pair(c, a);
			ptpairs.push_back(pr1);
			ptpairs.push_back(pr2);
			ptpairs.push_back(pr3);
			vector<int>pms;
			for (int ptk = 0; ptk < 3; ptk++) {
				int pt1 = ptpairs[ptk].first;
				int pt2 = ptpairs[ptk].second;
				string ptnow = to_string(pt1) + "a" + to_string(pt2);
				pms.push_back(_Midpts[ptnow]);
			}
			int am = _ord3D[a];
			int bm = _ord3D[b];
			int cm = _ord3D[c];
			vector<int>pntpos;
			vector<int>ords;
			ords = { am,bm,cm };
			vector<Vector2d>_quadpoints;
			vector<Real>_quadweights;
			vector<Vector3d>nodescur;
			vector<Vector3d>nodesref;
			vector<Real>N;
			vector<Vector2d>DN;
			Vector3d A = { _PosOrdrd[2 * am + 0],_PosOrdrd[2 * am + 1],0.0 };
			Vector3d B = { _PosOrdrd[2 * bm + 0],_PosOrdrd[2 * bm + 1],0.0 };
			Vector3d C = { _PosOrdrd[2 * cm + 0],_PosOrdrd[2 * cm + 1],0.0 };
			D = { _PosOrdrd[2 * pms[0] + 0],_PosOrdrd[2 * pms[0] + 1],0.0 };
			E = { _PosOrdrd[2 * pms[1] + 0],_PosOrdrd[2 * pms[1] + 1],0.0 };
			FF = { _PosOrdrd[2 * pms[2] + 0],_PosOrdrd[2 * pms[2] + 1],0.0 };
			
			vector<Vector3d>NPs;
			Vector3d NP1, NP2, NP3;
			NP1 = (_MPos3D[a % _npart] + _MPos3D[b % _npart]) / 2.0;
			NP2 = (_MPos3D[b % _npart] + _MPos3D[c % _npart]) / 2.0;
			NP3 = (_MPos3D[c % _npart] + _MPos3D[a % _npart]) / 2.0;
			NPs = { NP1,NP2,NP3 };

			for (int IP = 0; IP < 3; IP++) {
				Real factor = L;
				if (_MPos3D[_Tris[trr][IP] % _npart][0] == 4.0 && doneforcing.find(ords[IP]) == doneforcing.end()) {
					doneforcing.insert(ords[IP]);
					if (_MPos3D[_Tris[trr][IP] % _npart][1] == 0.0 || _MPos3D[_Tris[trr][IP] % _npart][1] == 4.0)factor = L / 2;
					//cout << ords[IP] <<" "<< factor << endl;
					Real Finalforce = Force * factor;
					Vector2d A = { _MPos3D[_Tris[trr][IP]][0],_MPos3D[_Tris[trr][IP]][1] };
					Vector2d B = { _PosOrdrd[2 * ords[IP]],_PosOrdrd[2 * ords[IP] + 1] };
					//EnergyTotal += -Finalforce * (_PosOrdrd[2 * ords[IP]] - _MPos3D[_Tris[trr][IP]][0]);
					//_GrdOrdrd[2 * ords[IP]] += -Finalforce;
				}
				factor = L;
				if (scheme) {
					if (NPs[IP][0] == 4.0 && doneforcing.find(pms[IP]) == doneforcing.end()) {
						doneforcing.insert(pms[IP]);
						//cout << NPs[IP] << endl;
						//cout << pms[IP] << endl;
						if (NPs[IP][1] == 0.0 || NPs[IP][1] == 4.0)factor = L / 2;
						//cout << factor << endl;
						Real Finalforce = Force * factor;
						Vector2d A = { NPs[IP][0], NPs[IP][1] };
						Vector2d B = { _PosOrdrd[2 * pms[IP]],_PosOrdrd[2 * pms[IP] + 1] };
						//EnergyTotal += -Finalforce * (_PosOrdrd[2 * pms[IP]] - NPs[IP][0]);
						//_GrdOrdrd[2 * pms[IP]] += -Finalforce;
					}
				}
			}

			if (parabola) {
				NP1[0] = (NP1[1] * NP1[1] + NP1[2] * NP1[2]) / 8.0;
				NP2[0] = (NP2[1] * NP2[1] + NP2[2] * NP2[2]) / 8.0;
				NP3[0] = (NP3[1] * NP3[1] + NP3[2] * NP3[2]) / 8.0;
			}
			vector<Vector3d>finalmids;
			for (int mids = 0; mids < 3; mids++) {
				Vector3d NPnow = NPs[mids];
				Real Radius = sqrt(NPnow[0] * NPnow[0] + NPnow[1] * NPnow[1] + NPnow[2] * NPnow[2]);
				Real theta = acos(NPnow[2] / Radius);
				Real phi = atan2(NPnow[1], NPnow[0]);
				finalmids.push_back({ 2.0 * sin(theta) * cos(phi),2.0 * sin(theta) * sin(phi) ,2.0 * cos(theta) });
			}
			
			if (scheme == 0) {
				pntpos = { 2 * am + 0 ,2 * am + 1,2 * bm + 0,2 * bm + 1,2 * cm + 0,2 * cm + 1 };
				nodescur = { A,B,C };
				nodesref = { _MPos3D[a % _npart],_MPos3D[b % _npart],_MPos3D[c % _npart] };
				_quadpoints.resize(1, { 1.0 / 3.0, 1.0 / 3.0 });
				_quadweights.resize(1, 0.5);
			}
			else if (scheme == 1) {
				pntpos = { 2 * am + 0 ,2 * am + 1,2 * bm + 0,2 * bm + 1,2 * cm + 0,2 * cm + 1 ,2 * pms[0] + 0,2 * pms[0] + 1,2 * pms[1] + 0,2 * pms[1] + 1,2 * pms[2] + 0,2 * pms[2] + 1 };
				nodescur = { A,B,C,D,E,FF };
				if(refdome)nodesref = { _MPos3D[a % _npart],_MPos3D[b % _npart],_MPos3D[c % _npart],finalmids[0],finalmids[1],finalmids[2] };
				else nodesref = { _MPos3D[a % _npart],_MPos3D[b % _npart],_MPos3D[c % _npart],NP1,NP2,NP3 };
				_quadpoints.resize(3, { 0.0,0.0 });
				_quadpoints[0] = { 1.0 / 6.0, 2.0 / 3.0 };
				_quadpoints[1] = { 1.0 / 6.0, 1.0 / 6.0 };
				_quadpoints[2] = { 2.0 / 3.0, 1.0 / 6.0 };
				_quadweights.resize(3, 1.0 / 6.0);
			}
			MatrixXd F(3, 3);
			vector<Vector3d> basis, dual, refbasis, refdual;
			basis.resize(2, { 0.0,0.0,0.0 });
			dual.resize(2, { 0.0,0.0,0.0 });
			refbasis.resize(2, { 0.0,0.0,0.0 });
			refdual.resize(2, { 0.0,0.0,0.0 });


			Real Aref = 0.0;
			Real Adef = 0.0;
			Real _Jtotal = 0.0;
			Real TRCtotal = 0.0;
			for (int qds = 0; qds < _quadpoints.size(); qds++) {
				vector<Vector3d>_n;
				_n.resize(3, { 0.0,0.0,0.0 });
				if (scheme == 0) {
					N = NTri3(_quadpoints[qds]);
					DN = DNTri3(_quadpoints[qds]);
				}
				if (scheme == 1) {
					N = NTri6(_quadpoints[qds]);
					DN = DNTri6(_quadpoints[qds]);
				}
				// Reference Configuration Basis and Duals
				refbasis = basiss(nodesref, DN);
				MatrixXd metrictensorref(2, 2);
				metrictensorref = metric(refbasis);
				Real G = sqrtdetmetric(metrictensorref);
				MatrixXd metrictensorinvref(2, 2);
				metrictensorinvref = metricinv(metrictensorref);
				refdual = duals(refbasis, metrictensorinvref);


				// Deformed Configuration Basis and Duals
				basis = basiss(nodescur, DN);
				MatrixXd metrictensor(2, 2);
				metrictensor = metric(basis);
				Real g = sqrtdetmetric(metrictensor);
				MatrixXd metrictensorinv(2, 2);
				metrictensorinv = metricinv(metrictensor);
				dual = duals(basis, metrictensorinv);

				// Weights
				Aref += double(_quadweights[qds] * G);
				Adef += double(_quadweights[qds] * g);
				MatrixXd Cc(3, 3);
				for (int I = 0; I < 3; I++) {
					for (int J = 0; J < 3; J++) {
						F(I, J) = 0.0;
						Cc(I, J) = 0.0;
					}
				}
				for (int alpha = 0; alpha < 2; alpha++) {
					for (int I = 0; I < 3; I++) {
						for (int J = 0; J < 3; J++) {
							F(I, J) += basis[alpha][I] * refdual[alpha][J];
						}
					}
				}
				Real weight = _quadweights[qds] * G;
				pair<Real, Real> ResE = EnergyEvansEnergies(F, "Both", 1.0, 1.0);
				Energy1 += weight * (ResE.first);
				Energy2 += weight * (ResE.second);
				EnergyTotal = Energy1 + Energy2;
				//EnergyTotal += weight*ResE.second;
				
				Cc = F.transpose() * F;
				Real _trC = Cc(0, 0) + Cc(1, 1) + Cc(2, 2);
				Real trCSquare = 0.0;
				for (int I = 0; I < 3; I++) {
					for (int k = 0; k < 3; k++) {
						trCSquare += Cc(I, k) * Cc(k, I);
					}
				}
				Real _Jfinal = sqrt((_trC * _trC - trCSquare) / 2.0);
				FinalJ += ((_Jfinal-1.0)* (_Jfinal - 1.0))*weight;
				FinalTrc += ((_trC / _Jfinal - 2.0) * (_trC / _Jfinal - 2.0)) * weight;
				_Jtotal += (_Jfinal - 1.0) * weight;
				TRCtotal += (_trC / _Jfinal - 2.0) * weight;
			}
			
			TotalDistortion += abs(Adef-Aref);
			TotalActual += Aref;
			TotalDeformed += Adef;
			Real Val = ((Adef/Aref)-1.0);
			Diffs.push_back(make_pair(Ts, Val));
			Trcs.push_back(make_pair(Ts, TRCtotal/Aref));
			Js.push_back(make_pair(Ts, _Jtotal/Aref));
		}
		cout << TotalDistortion << " " << TotalActual << " " << TotalDeformed << endl;
		cout << sqrt(FinalJ/TotalActual) << " " << sqrt(FinalTrc / TotalActual) << " " << EnergyTotal << endl;
		DiffTris = Diffs;
		TRCS = Trcs;
		JS = Js;
	}

	
}