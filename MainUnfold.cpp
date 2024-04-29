#include "Unfolding.h"
#include "ParticlePositions.h"
#include "DefectsPositions.h"
#include "Triangulation.h"
#include "UnfoldingPath.h"
#include "UnfoldTriangles.h"
#include "Test.h"
#include <LBFGSB.h>
#include <Eigen>

using namespace unfolding;
vector<vector<int>>Tris;
vector<Real>EnergyF;
map<int, int>ordr;
int iter = 1;
int nprt;
string s;
VectorXd computeGrd(VectorXd x) {
	vector<Vector2d>Pos;
	int n = (x.size());
	Pos.resize(int(n / 2), { 0.0,0.0 });
	for (int i = 0; i<int(n / 2); i++) {
		Pos[i][0] = x[2 * i + 0];
		Pos[i][1] = x[2 * i + 1];
	}
	Body bdcomp(Pos, Tris, ordr, s, 0, 0, 0, nprt);
	vector<double>gr = bdcomp.getGrd();
	VectorXd GR;
	GR.resize(n);
	for (int i = 0; i < n; i++) {
		GR[i] = gr[i];
	}
	return GR;
}
double computeEnrg(VectorXd x) {
	vector<Vector2d>Pos;
	int n = (x.size());
	Pos.resize(int(n / 2), { 0.0,0.0 });
	for (int i = 0; i<int(n / 2); i++) {
		Pos[i][0] = x[2 * i + 0];
		Pos[i][1] = x[2 * i + 1];
	}
	Body bdcomp(Pos, Tris, ordr, s, 0, 0, 0, nprt);
	//cout << Tris.size() << " " << Pos.size() << endl;
	double enrg = bdcomp.getenrg();
	return enrg;
}
double foo(const VectorXd& x, VectorXd& grad)
{
	cout << "iter" << iter << endl;
	iter++;
	double f;
	f = computeEnrg(x);
	cout << "enrg = " << f << endl;
	EnergyF.push_back(f);
	grad = computeGrd(x);
	//cout << "grad = " << grad << endl;
	return f;
}
int main() {
	int numberof_particles = 376;
	/*string inputparts = "parabolacirlcetestexp.vtk";						//Input file for particle positions
	string inputcon = "parabolacirlcetestexpdata.txt";						//Input file for particle connectivty
	string input3D = "parabolacirlcetestexp.vtk";
	string inputtriscon = "NULL";
	string inputdef = "parabolacircledefects.txt";
	string inputtri = "parabolacircletris.txt";*/

	string inputparts = "subdiviParaexp.vtk";						//Input file for particle positions
	string inputcon = "subdiviParatrisexp.vtk";						//Input file for particle connectivty
	string input3D = "subdiviParaexp.vtk";
	string inputtriscon = "subdiviTris.txt";
	string inputdef = "parabolacircledefects.txt";
	string inputtri = "parabolacircletris.txt";

	bool Minim = 1;
	int mststeps = 10;
	//Search Parameters for Defects Connectivity
	Real RC = 4.0, DELR = 2.0, DR = 1.0, INITIALR = 1.0;
	cout << inputtriscon << endl;
	ParticlePositions Part(numberof_particles,inputparts,inputcon,input3D, inputtriscon);
	vector<Vector2d>Positions;
	map<int, int>ord3d;
	int np = Part.n_part;
	for (int i = 0; i < np; i++) {
		Positions.push_back({ Part.pos[i][1],Part.pos[i][2] });
		ord3d[i] = i;
	}
	map<vector<int>, bool>donetris;
	
	/*for (int i = 0; i < np; i++) {
		cout << i<<" "<<Part.pos_neighbor[i].size() << endl;
		for (int j=0; j < Part.pos_neighbor[i].size(); j++) {
			if(Part.pos_neighbor[i].size() == 4 && j == 3)continue;
			vector<int>T;
			T.push_back(i);
			T.push_back(Part.pos_neighbor[i][j]);
			T.push_back(Part.pos_neighbor[i][(j+1)% Part.pos_neighbor[i].size()]);
			sort(T.begin(), T.end());
			if (donetris[T] == 0) {
				Tris.push_back(T);
				donetris[T] = 1;
				//print_vector(T);
			}
		}
	}*/
	Tris = Part.Triangles;
	Body bd(Positions, Tris, ord3d, inputparts, 0, 0, 0, np);
	//bd.checkconsistencyNeoHookean(1e-7, 1e-6);
	VectorXd xo;
	xo.resize(np * 2);
	vector<double>Pos2;
	Pos2.clear();
	Pos2 = bd.getPos();
	cout << "HEY" << endl;
	cout << bd.getenrg() << "Energy" << endl;
	for (int l = 0; l < 2 * np; l++) {
		xo[l] = Pos2[l];
	}
	ordr = ord3d;
	s = inputparts;
	nprt = np;
	LBFGSBParam<double> param;
	param.epsilon = 1e-15;
	param.delta = 1e-15;
	LBFGSBSolver<double> solver(param);
	double fx;
	iter = 1;
	VectorXd lb = VectorXd::Constant(2 * np,  -10000.0);
	VectorXd ub = VectorXd::Constant(2 * np,  10000.0);
	lb[0] = 0.0;
	ub[0] = 0.0;
	lb[1] = 0.0;
	ub[1] = 0.0;
	lb[3] = 0.0;
	ub[3] = 0.0;
	int niter = solver.minimize(foo, xo, fx, lb, ub);
	EnergyF.push_back(fx);
	for (int l = 0; l < np; l++) {
		Positions[l][0] = xo[2 * l + 0];
		Positions[l][1] = xo[2 * l + 1];
	}
	
	Body bd2(Positions, Tris, ord3d, inputparts, 0, 0, 0, np);
	vector< pair<vector<int>, Real>> difs;
	difs = bd2.checkDiffTri();
	pair<vector< pair<vector<int>, Real>>, vector< pair<vector<int>, Real>>> resis;
	resis = bd2.plotNeoHookeanGrd2D();
	//bd2.checkconsistency(1e-6, 1e-6);
	ofstream f;
	f.open("TestNeoHookeanParabola" + Part.input_particlefile);
	f << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS ";
	f << Positions.size();
	f << " float\n";
	for (int i = 0; i < Positions.size(); i++) {
		f << Positions[i][0] << " " << Positions[i][1] << " " << 0 << endl;
	}
	f << endl;
	f << "CELLS ";
	f << Tris.size() << " " << 4 * Tris.size() << endl;
	for (int i = 0; i < Tris.size(); i++) {
		f << 3 << " " << Tris[i][0] << " " << Tris[i][1] << " " << Tris[i][2] << endl;
	}
	f << endl;
	f << "CELL_TYPES ";
	f << Tris.size() << endl;
	for (int i = 0; i < Tris.size(); i++) {
		f << 5 << endl;
	}
	f << "\nCELL_DATA " << Tris.size() << "\n";
	f << "SCALARS neighbour float\n";
	f << "LOOKUP_TABLE default\n";
	Real totaldistortion = 0.0;
	for (int i = 0; i < Tris.size(); i++) {
		vector<int>maint = { 0,0,0 };
		for (int j = 0; j < 3; j++) {
			maint[j] = Tris[i][j];
		}
		sort(maint.begin(), maint.end());
		print_vector(maint);
		bool found = false;
		for (int j = 0; j < difs.size(); j++) {
			vector<int>ts = difs[j].first;
			sort(ts.begin(), ts.end());
			if (maint == ts) {
				f << difs[j].second << endl;
				cout << difs[j].second << endl;
				found = true;
				break;
			}
		}
	}
	f.close();

	ofstream f1;
	f1.open("TestNeoHookeanParabolaI1bar" + Part.input_particlefile);
	f1 << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS ";
	f1 << Positions.size();
	f1 << " float\n";
	for (int i = 0; i < Positions.size(); i++) {
		f1 << Positions[i][0] << " " << Positions[i][1] << " " << 0 << endl;
	}
	f1 << endl;
	f1 << "CELLS ";
	f1 << Tris.size() << " " << 4 * Tris.size() << endl;
	for (int i = 0; i < Tris.size(); i++) {
		f1 << 3 << " " << Tris[i][0] << " " << Tris[i][1] << " " << Tris[i][2] << endl;
	}
	f1 << endl;
	f1 << "CELL_TYPES ";
	f1 << Tris.size() << endl;
	for (int i = 0; i < Tris.size(); i++) {
		f1 << 5 << endl;
	}
	f1 << "\nCELL_DATA " << Tris.size() << "\n";
	f1 << "SCALARS neighbour float\n";
	f1 << "LOOKUP_TABLE default\n";
	for (int i = 0; i < Tris.size(); i++) {
		vector<int>maint = { 0,0,0 };
		for (int j = 0; j < 3; j++) {
			maint[j] = Tris[i][j];
		}
		sort(maint.begin(), maint.end());
		print_vector(maint);
		bool found = false;
		for (int j = 0; j < resis.second.size(); j++) {
			vector<int>ts = resis.second[j].first;
			sort(ts.begin(), ts.end());
			if (maint == ts) {
				f1 << resis.second[j].second << endl;
				found = true;
				break;
			}
		}
	}
	f1.close();

	ofstream f3;
	f3.open("TestNeoHookeanParabolaI3" + Part.input_particlefile);
	f3 << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS ";
	f3 << Positions.size();
	f3 << " float\n";
	for (int i = 0; i < Positions.size(); i++) {
		f3 << Positions[i][0] << " " << Positions[i][1] << " " << 0 << endl;
	}
	f3 << endl;
	f3 << "CELLS ";
	f3 << Tris.size() << " " << 4 * Tris.size() << endl;
	for (int i = 0; i < Tris.size(); i++) {
		f3 << 3 << " " << Tris[i][0] << " " << Tris[i][1] << " " << Tris[i][2] << endl;
	}
	f3 << endl;
	f3 << "CELL_TYPES ";
	f3 << Tris.size() << endl;
	for (int i = 0; i < Tris.size(); i++) {
		f3 << 5 << endl;
	}
	f3 << "\nCELL_DATA " << Tris.size() << "\n";
	f3 << "SCALARS neighbour float\n";
	f3 << "LOOKUP_TABLE default\n";
	for (int i = 0; i < Tris.size(); i++) {
		vector<int>maint = { 0,0,0 };
		for (int j = 0; j < 3; j++) {
			maint[j] = Tris[i][j];
		}
		sort(maint.begin(), maint.end());
		print_vector(maint);
		bool found = false;
		for (int j = 0; j < resis.first.size(); j++) {
			vector<int>ts = resis.first[j].first;
			sort(ts.begin(), ts.end());
			if (maint == ts) {
				f3 << resis.first[j].second << endl;
				found = true;
				break;
			}
		}
	}
	f3.close();
	for (int i = 0; i < resis.first.size(); i++) {
		cout << difs[i].second / resis.first[i].second << endl;
	}
	return 0;
	/*
	vector<Vector2d>Positions;
	Positions.push_back({ 0.0,0.0 });
	Positions.push_back({ 0.0,1.0 });
	Positions.push_back({ 1.0,1.0 });
	Positions.push_back({ 1.0,0.0 });
	Tris.push_back({ 0,1,2 });
	Tris.push_back({ 0,2,3 });
	
	map<int, int>ord3d;
	ord3d[0] = 0;
	ord3d[1] = 1;
	ord3d[2] = 2;
	ord3d[3] = 3;
	int np = 4;
	string inp = "testtris.vtk";
	Body bd(Positions, Tris, ord3d, inp, 0, 0, 0, np);
	bd.checkconsistencyNeoHookean2D(1e-6, 1e-6);
	VectorXd xo;
	xo.resize(np * 2);
	vector<double>Pos2;
	Pos2.clear();
	Pos2 = bd.getPos();
	cout << bd.getenrg() << "Energy" << endl;
	for (int l = 0; l < 2 * np; l++) {
		xo[l] = Pos2[l];
	}
	ordr = ord3d;
	s = inp;
	nprt = np;
	LBFGSBParam<double> param;
	param.epsilon = 1e-12;
	param.max_iterations = 100;
	LBFGSBSolver<double> solver(param);
	double fx;
	iter = 1;
	VectorXd lb = VectorXd::Constant(2*np, -1000.0);
	VectorXd ub = VectorXd::Constant(2 * np, 1000.0);
	lb[0] = 0.0;
	ub[0] = 0.0;
	lb[1] = 0.0;
	ub[1] = 0.0;
	int niter = solver.minimize(foo, xo, fx,lb,ub);
	EnergyF.push_back(fx);
	for (int l = 0; l < np; l++) {
		Positions[l][0] = xo[2 * l + 0];
		Positions[l][1] = xo[2 * l + 1];
	}

	Body bd2(Positions, Tris, ord3d, inp, 0, 0, 0, np);
	vector< pair<vector<int>, Real>> difs;
	difs = bd2.checkDiffTri();

	ofstream f;
	f.open("TestNeoHookeanRectangle" + Part.input_particlefile);
	f << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS ";
	f << Positions.size();
	f << " float\n";
	for (int i = 0; i < Positions.size(); i++) {
		f << Positions[i][0] << " " << Positions[i][1] << " " << 0 << endl;
	}
	f << endl;
	f << "CELLS ";
	f << Tris.size() << " " << 4 * Tris.size() << endl;
	for (int i = 0; i < Tris.size(); i++) {
		f << 3 << " " << Tris[i][0] << " " << Tris[i][1] << " " << Tris[i][2] << endl;
	}
	f << endl;
	f << "CELL_TYPES ";
	f << Tris.size() << endl;
	for (int i = 0; i < Tris.size(); i++) {
		f << 5 << endl;
	}
	f << "\nCELL_DATA " << Tris.size() << "\n";
	f << "SCALARS neighbour float\n";
	f << "LOOKUP_TABLE default\n";
	Real totaldistortion = 0.0;
	for (int i = 0; i < Tris.size(); i++) {
		vector<int>maint = { 0,0,0 };
		for (int j = 0; j < 3; j++) {
			maint[j] = Tris[i][j];
		}
		sort(maint.begin(), maint.end());
		//print_vector(maint);
		bool found = false;
		for (int j = 0; j < difs.size(); j++) {
			vector<int>ts = difs[j].first;
			sort(ts.begin(), ts.end());
			if (maint == ts) {
				f << difs[j].second << endl;
				found = true;
				break;
			}
		}
	}
	f.close();
	return 0;*/
	DefectsPositions Dfct(Part,RC, DELR, DR, INITIALR, inputdef);
	Triangulation Tri(Part, Dfct, inputtri);
	UnfoldingPath UPath(Tri,Part,mststeps);
	UnfoldTriangles UnTri(UPath, Part, Minim);


}