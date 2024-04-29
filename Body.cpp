#include "Body.h"

namespace unfolding {
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

	void Body::setEdgGrd() {
		for (int i = 0; i < _Edges.size(); i++) {
			int a = _Edges[i].first;
			int b = _Edges[i].second;
			int c = _EdgeTri[i];
			Real l = dist(_MPos3D[a % _npart], _MPos3D[b % _npart]);
			int am = _ord3D[a];
			int bm = _ord3D[b];
			Real Al = trianglearea(_MPos3D[a % _npart], _MPos3D[b % _npart], _MPos3D[c % _npart]);
			//cout << a <<" "<< b <<" "<< c << endl;
			Vector2d A = { _PosOrdrd[2 * am + 0],_PosOrdrd[2 * am + 1] };
			Vector2d B = { _PosOrdrd[2 * bm + 0],_PosOrdrd[2 * bm + 1] };

			Real d = dist(A, B);

			Real enrg = d - l;
			//cout << d << " "<<l << endl;
			_Energy += 0.5*pow(enrg, 2.0)*Al;

			_GrdOrdrd[2 * am + 0] += Al*(enrg) * (A[0] - B[0]) / d;
			_GrdOrdrd[2 * am + 1] += Al*(enrg) * (A[1] - B[1]) / d;
			_GrdOrdrd[2 * bm + 0] += Al * (enrg) * (B[0] - A[0]) / d;
			_GrdOrdrd[2 * bm + 1] += Al * (enrg) * (B[1] - A[1]) / d;

		}
		//print_vector(_GrdOrdrd);
	}

	void Body::resetGrd() {
		_GrdOrdrd.resize(2 * n, 0.0);
	}

	double Body::getenrg() {
		return _Energy;
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
	MatrixXd Pmatrix(MatrixXd F) {
		MatrixXd CC = F.transpose() * F;
		Real I1 = CC.trace();
		Real I1bar = pow(F.determinant(), -2.0 / 3.0) * I1;
		Real I3 = CC.determinant();
		MatrixXd P = 4.0 * (I3 * I3 - 1.0 / (I3 * I3)) * F.transpose().inverse() + 2.0 * pow(F.determinant(), -2.0 / 3.0) * (F - F.transpose().inverse() * I1 / 3.0);
		//MatrixXd P = 2.0 * F;
		return P;
	}
	Real FEnergy(MatrixXd F) {
		MatrixXd CC = F.transpose() * F;
		Real I1 = CC.trace();
		Real I1bar = pow(F.determinant(), -2.0 / 3.0) * I1;
		Real I3 = CC.determinant();
		Real enrg = (I3 * I3 + 1.0 / (I3 * I3) - 2.0) + (I1bar - 3.0);
		//Real enrg = I1;
		return enrg;
	}
	MatrixXd Pmatrix2D(MatrixXd F) {
		MatrixXd CC = F.transpose() * F;
		Real I1 = CC.trace();
		Real I1bar = pow(F.determinant(), -1.0) * I1;
		Real I3 = CC.determinant();
		MatrixXd P = 4.0 * (I3 * I3 - 1.0 / (I3 * I3)) * F.transpose().inverse() + pow(F.determinant(), -1.0) * (2.0*F - F.transpose().inverse() * I1);
		//MatrixXd P = 2.0 * F;
		return P;
	}
	Real FEnergy2D(MatrixXd F) {
		MatrixXd CC = F.transpose() * F;
		Real I1 = CC.trace();
		Real I1bar = pow(F.determinant(), -1.0) * I1;
		Real I3 = CC.determinant();
		Real enrg = (I3 * I3 + 1.0 / (I3 * I3) - 2.0) + (I1bar - 2.0);
		//Real enrg = I1;
		return enrg;
	}
	void Body::setNeoHookeanGrd() {
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
			MatrixXd F(3, 3);
			Real N1x, N1y, N2x, N2y, N3x, N3y;
			if (AB.cross(AC)[2] < 0) {
				Vector3d Temp = B;
				B = C;
				C = Temp;
			}
			vector<Vector3d>ABC = { A,B,C };
			b1 = dist(_MPos3D[a % _npart], _MPos3D[b % _npart]);
			c2 = 2.0 * Al / b1;
			c1 = sqrt(pow(dist(_MPos3D[a % _npart], _MPos3D[c % _npart]), 2.0) - c2 * c2);
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
			F(0, 2) = 0.0;
			F(1, 2) = 0.0;
			F(2, 0) = 0.0;
			F(2, 1) = 0.0;
			F(2, 2) = 1.0;
			_Energy += FEnergy(F) * Al;
			MatrixXd P = Pmatrix(F);
			_GrdOrdrd[2 * am + 0] += Al * (P(0, 0) * N1x + P(0, 1) * N1y);
			_GrdOrdrd[2 * am + 1] += Al * (P(1, 0) * N1x + P(1, 1) * N1y);
			_GrdOrdrd[2 * bm + 0] += Al * (P(0, 0) * N2x + P(0, 1) * N2y);
			_GrdOrdrd[2 * bm + 1] += Al * (P(1, 0) * N2x + P(1, 1) * N2y);
			_GrdOrdrd[2 * cm + 0] += Al * (P(0, 0) * N3x + P(0, 1) * N3y);
			_GrdOrdrd[2 * cm + 1] += Al * (P(1, 0) * N3x + P(1, 1) * N3y);
		}
	}
	void Body::setNeoHookeanGrd2D() {
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
			Real fact = 1.0;
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
			b1 = 0.34;
			c1 = 0.34 / 2.0;
			c2 = fact * b1 * sqrt(3.0) / 2.0;
			//c2 = fact*2.0 * Al / b1;
			//c1 = sqrt(pow(dist(_MPos3D[a % _npart], _MPos3D[c % _npart]), 2.0) - c2 * c2);
			//cout << b1 << " " << c1 << " " << c2 << endl;
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
			//cout << F << endl;
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
	pair<vector< pair<vector<int>, Real>>, vector< pair<vector<int>, Real>>> Body::plotNeoHookeanGrd2D() {
		vector< pair<vector<int>, Real>> I3s;
		vector< pair<vector<int>, Real>> I1s;
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
			Real fact = 1.0;
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
			b1 = 0.34;
			c1 = 0.34 / 2.0;
			c2 = fact * b1 * sqrt(3.0) / 2.0;
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
			vector<int>Ts;
			Ts = { a % _npart,b % _npart,c % _npart };
			MatrixXd CC = F.transpose() * F;
			Real I1 = CC.trace();
			Real I1bar = pow(F.determinant(), -1.0) * I1;
			Real I3 = CC.determinant();
			I3s.push_back(make_pair(Ts, I3));
			I1s.push_back(make_pair(Ts, I1bar));
		}
		return make_pair(I3s, I1s);
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

			//_Energy += pow(enrg, 2.0);
		}
		//cout << "Bad Edges: Inner :  " << countr << endl;
		//cout << "Bad Edges: Boundary :  " << countrb << endl;
		return difs;
	}

	vector< pair<vector<int>, Real>> Body::checkDiffTri() {
		vector< pair<vector<int>, Real>> Diffs;
		int countr = 0;
		for (int i = 0; i < _Tris.size(); i++) {
			int a = _Tris[i][0];
			int b = _Tris[i][1];
			int c = _Tris[i][2];
			Real Al = trianglearea(_MPos3D[a % _npart], _MPos3D[b % _npart], _MPos3D[c % _npart]);

			int am = _ord3D[a];
			int bm = _ord3D[b];
			int cm = _ord3D[c];

			Vector3d A = { _PosOrdrd[2 * am + 0],_PosOrdrd[2 * am + 1] ,0.0 };
			Vector3d B = { _PosOrdrd[2 * bm + 0],_PosOrdrd[2 * bm + 1] ,0.0 };
			Vector3d C = { _PosOrdrd[2 * cm + 0],_PosOrdrd[2 * cm + 1] ,0.0 };

			Real At = trianglearea(A, B, C);
			vector<int>Ts;
			Ts = { a % _npart,b % _npart,c % _npart };
			Al = 0.050056;
			Real Val = ((At - Al) / Al);
			Diffs.push_back(make_pair(Ts, Val));
		}
		return Diffs;
	}
	Real Body::checkTotalSurfaceTri() {
		Real total = 0.0;
		int countr = 0;
		for (int i = 0; i < _Tris.size(); i++) {
			int a = _Tris[i][0];
			int b = _Tris[i][1];
			int c = _Tris[i][2];
			Real Al = trianglearea(_MPos3D[a % _npart], _MPos3D[b % _npart], _MPos3D[c % _npart]);
			total += Al; /* / Al);*/
		}
		//cout << "Bad Triangles: " << countr << endl;
		return total;
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


	void Body::checkconsistency(Real h, Real tol) {
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
			//cout << l << " " << d << endl;
			//Real enrg = (pow(d, 2.0) - pow(l, 2.0));
			//_Energy += pow(enrg, 2.0);

			Real enrg = d - l;
			//cout << d << " "<<l << endl;
			_Energy += 0.5 * Al * pow(enrg, 2.0);

			//_GrdOrdrd[2 * am + 0] += (enrg) * pow(d,-0.5) * (A[0] - B[0]);
			//_GrdOrdrd[2 * am + 1] += (enrg) * pow(d,-0.5)* (A[1] - B[1]);
			//_GrdOrdrd[2 * bm + 0] += (enrg) * pow(d, -0.5) * (B[0] - A[0]);
			//_GrdOrdrd[2 * bm + 1] += (enrg) * pow(d, -0.5) * (B[1] - A[1]);

			_GrdOrdrd[2 * am + 0] += Al * (enrg)*(A[0] - B[0])/d;
			_GrdOrdrd[2 * am + 1] += Al * (enrg)*(A[1] - B[1])/d;
			_GrdOrdrd[2 * bm + 0] += Al * (enrg)*(B[0] - A[0])/d;
			_GrdOrdrd[2 * bm + 1] += Al * (enrg)*(B[1] - A[1])/d;

			Vector2d Xplus = { _PosOrdrd[2 * am + 0] + h,_PosOrdrd[2 * am + 1] };
			Vector2d Xminus = { _PosOrdrd[2 * am + 0] - h,_PosOrdrd[2 * am + 1] };

			Real dplus = dist(Xplus, B);
			Real dminus = dist(Xminus, B);
			Real Eplus = 0.5*pow(dplus-l, 2.0);
			Real Eminus = 0.5*pow(dminus - l, 2.0);

			Real error = ((Eplus - Eminus) / (2.0 * h)) - (enrg) * (A[0] - B[0]) / d;
			Err.push_back(error);
			//cout << h << " " << Xplus << " x " << Xminus << " x " << A << " a " << dplus << " " << dminus << " " << Eplus << " " << Eminus << " " << enrg<<" " << 4.0 * (enrg) * (A[0] - B[0]) << endl;
			//cout << abs(error) << " e " << i << " i " << _Edges.size() << endl;
			if ((abs(error)) > tol)cout << "problem" << endl;
		}
		Real Toter = 0;
		Real maxFault = -10;
		for (int i = 0; i < Err.size(); i++) {
			Toter = Err[i] * Err[i];
			//Toter = max(abs(Err[i]),maxFault);
		}
		Toter = sqrt(Toter);
		//cout << Toter << endl;
	}
	void Body::checkconsistencyNeoHookean2D(Real h, Real tol) {
		for (int i = 0; i < _Tris.size(); i++) {
			int a = _Tris[i][0];
			int b = _Tris[i][1];
			int c = _Tris[i][2];
			//cout << a << " " << b << " " << c << endl;
			int am = _ord3D[a];
			int bm = _ord3D[b];
			int cm = _ord3D[c];
			Real Al = abs(trianglearea(_MPos3D[a % _npart], _MPos3D[b % _npart], _MPos3D[c % _npart]));
			//cout << am << " " << bm << " " << cm << endl;
			//cout << 2 * am << " " << 2 * am + 1 << " " << 2 * bm << " " << 2 * bm + 1 << " " << endl;
			Real b1, c2, c1;
			Vector3d A = { _PosOrdrd[2 * am + 0],_PosOrdrd[2 * am + 1],0.0 };
			Vector3d B = { _PosOrdrd[2 * bm + 0],_PosOrdrd[2 * bm + 1],0.0 };
			Vector3d C = { _PosOrdrd[2 * cm + 0],_PosOrdrd[2 * cm + 1],0.0 };
			Vector3d AB = B - A;
			Vector3d AC = C - A;
			MatrixXd F(2,2);
			Real N1x, N1y, N2x, N2y, N3x, N3y;
			bool change = 0;
			Real fact = 0.0;
			if (AB.cross(AC)[2] < 0) {
				/*Vector3d Temp = B;
				B = C;
				C = Temp;
				change = 1;*/
				fact = -1.0;
			}
			vector<Vector3d>ABC = { A,B,C };
			Real dism = ceil(dist(_MPos3D[a % _npart], _MPos3D[b % _npart]) * 1000.0) / 1000.0;
			Real dis1 = ceil(dist(_MPos3D[a % _npart], _MPos3D[c % _npart]) * 1000.0) / 1000.0;
			Real dis2 = ceil(dist(_MPos3D[b % _npart], _MPos3D[c % _npart]) * 1000.0) / 1000.0;
			Real ah = (dis1 * dis1 - dis2 * dis2 + dism * dism) / (2 * dism);
			Real ha = sqrt(dis1 * dis1 - ah * ah);
			Vector3d rootpos = { 0.0,0.0,0.0 };
			Vector3d leverpos = { dism,0.0,0.0 };
			Vector3d PP0 = { rootpos[0] + ah * (leverpos[0] - rootpos[0]) / dism,rootpos[1] + ah * (leverpos[1] - rootpos[1]) / dism,0 };
			Vector3d PP1 = { PP0[0] + ha * (leverpos[1] - rootpos[1]) / dism,PP0[1] - ha * (leverpos[0] - rootpos[0]) / dism ,0 };
			Vector3d PP2 = { PP0[0] - ha * (leverpos[1] - rootpos[1]) / dism,PP0[1] + ha * (leverpos[0] - rootpos[0]) / dism ,0 };
			//cout << dism << " " << dis1 << " " << dis2 << endl;
			//cout << PP1 << " second " << PP2 << endl;
			//cout << fact << endl;
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
			//cout << b1 << " " << c1 << " " << c2 << endl;
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
			//cout << F << endl;
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
	void Body::checkconsistencyNeoHookean(Real h, Real tol) {
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
			MatrixXd F(3,3);
			Real N1x, N1y, N2x, N2y, N3x, N3y;
			if (AB.cross(AC)[2] <0) {
				Vector3d Temp = B;
				B = C;
				C = Temp;
			}
			vector<Vector3d>ABC = { A,B,C };
			b1 = dist(_MPos3D[a % _npart], _MPos3D[b % _npart]);
			c2 = 2.0 * Al / b1;
			c1 = sqrt(pow(dist(_MPos3D[a % _npart], _MPos3D[c % _npart]), 2.0) - c2 * c2);
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
			F(0, 2) = 0.0;
			F(1, 2) = 0.0;
			F(2, 0) = 0.0;
			F(2, 1) = 0.0;
			F(2, 2) = 1.0;
			//F(0, 0) = 2.0;
			//F(0, 1) = 0.0;
			//F(1, 0) = 0.0;
			//F(1, 1) = 1.0;
			//cout << enrg << endl;
			_Energy += FEnergy2D(F) * Al;
			MatrixXd P=Pmatrix2D(F);
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					MatrixXd Fp = F;
					MatrixXd Fm = F;
					Fp(i, j) += h;
					Fm(i, j) -= h;
					Real Ep = FEnergy2D(Fp);
					Real Em = FEnergy2D(Fm);
					Real errM = (Ep - Em) / (2.0 * h) - P(i, j);
					if (abs(errM) > tol)cout << "problem Material" << setprecision(10)<< abs(errM) << endl;
				}
			}
			vector<vector<Real>>Derivatives;
			Derivatives.push_back({ 0.0,0.0 });
			Derivatives.push_back({ 0.0,0.0 });
			Derivatives.push_back({ 0.0,0.0 });
			Derivatives[0][0]=  Al * (P(0, 0) * N1x + P(0, 1) * N1y);
			Derivatives[0][1] = Al * (P(1, 0) * N1x + P(1, 1) * N1y);
			Derivatives[1][0] = Al * (P(0, 0) * N2x + P(0, 1) * N2y);
			Derivatives[1][1] = Al * (P(1, 0) * N2x + P(1, 1) * N2y);
			Derivatives[2][0] = Al * (P(0, 0) * N3x + P(0, 1) * N3y);
			Derivatives[2][1] = Al * (P(1, 0) * N3x + P(1, 1) * N3y);
			for (int i = 0; i < 3; i++) {
				Vector3d Now = ABC[i];
				for (int j = 0; j < 2; j++) {
					ABC[i][j] += h;
					MatrixXd Fplus(3,3);
					Fplus(0, 0) = N1x * ABC[0][0] + N2x * ABC[1][0];
					Fplus(0, 1) = N1y * ABC[0][0] + N2y * ABC[1][0] + N3y * ABC[2][0];
					Fplus(1, 0) = N1x * ABC[0][1] + N2x * ABC[1][1];
					Fplus(1, 1) = N1y * ABC[0][1] + N2y * ABC[1][1] + N3y * ABC[2][1];
					Fplus(0, 2) = 0.0;
					Fplus(1, 2) = 0.0;
					Fplus(2, 0) = 0.0;
					Fplus(2, 1) = 0.0;
					Fplus(2, 2) = 1.0;
					ABC[i][j] -= 2.0*h;
					MatrixXd Fminus(3, 3);
					Fminus(0, 0) = N1x * ABC[0][0] + N2x * ABC[1][0];
					Fminus(0, 1) = N1y * ABC[0][0] + N2y * ABC[1][0] + N3y * ABC[2][0];
					Fminus(1, 0) = N1x * ABC[0][1] + N2x * ABC[1][1];
					Fminus(1, 1) = N1y * ABC[0][1] + N2y * ABC[1][1] + N3y * ABC[2][1];
					Fminus(0, 2) = 0.0;
					Fminus(1, 2) = 0.0;
					Fminus(2, 0) = 0.0;
					Fminus(2, 1) = 0.0;
					Fminus(2, 2) = 1.0;
					Real Eplus = Al * (FEnergy2D(Fplus));
					Real Eminus = Al * (FEnergy2D(Fminus));
					ABC[i] = Now;
					Real errE = (Eplus - Eminus) / (2.0 * h) - Derivatives[i][j];
					if (abs(errE) > tol)cout << "problem Element" << setprecision(10) << abs(errE) << endl;
				}
			}
			_GrdOrdrd[2 * am + 0] += Derivatives[0][0];
			_GrdOrdrd[2 * am + 1] += Derivatives[0][1];
			_GrdOrdrd[2 * bm + 0] += Derivatives[1][0];
			_GrdOrdrd[2 * bm + 1] += Derivatives[1][1];
			_GrdOrdrd[2 * cm + 0] += Derivatives[2][0];
			_GrdOrdrd[2 * cm + 1] += Derivatives[2][1];


		}
	}
}