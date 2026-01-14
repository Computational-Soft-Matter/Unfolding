#include "Unfolding.h"
#include "ParticlePositions.h"
#include "DefectsPositions.h"
#include "Triangulation.h"
#include "UnfoldingPath.h"
#include "UnfoldTriangles.h"
#include <LBFGSB.h>
#include <Eigen>

using namespace unfolding;

//.................Helper Functions for Distortion Minimization..............//
vector<vector<int>>Tris;
vector<Real>EnergyF;
map<int, int>ordr;
map<string, int>stordr;
int iter = 1;
int nprt;
string s;
pair<double, VectorXd> computes(VectorXd x) {
	vector<Vector2d>Pos;
	int n = (x.size());
	Pos.resize(int(n / 2), { 0.0,0.0 });
	for (int i = 0; i<int(n / 2); i++) {
		Pos[i][0] = x[2 * i + 0];
		Pos[i][1] = x[2 * i + 1];
	}
	Body bdcomp(Pos, Tris, ordr, stordr, s, 0, nprt);
	double enrg = bdcomp.getenrg();
	vector<double>gr = bdcomp.getGrd();
	VectorXd GR;
	GR.resize(n);
	for (int i = 0; i < n; i++) {
		GR[i] = gr[i];
	}
	return make_pair(enrg, GR);
}
double foo(const VectorXd& x, VectorXd& grad)
{
	iter++;
	pair<double, VectorXd>computeds = computes(x);
	double f;
	f = computeds.first;
	cout << iter << "*************************************** " << f << endl;
	EnergyF.push_back(f);
	grad = computeds.second;
	return f;
}

int main() {
	//Number of particles in the given 3D structure
	int numberof_particles = 212;

	//Input file for particle positions
	string inputparts;

	//Input file for particle connectivty
	string inputcon;

	//Input file for larger triangle connectivity
	//If this file is given, then the code do not search for defects or build large triangles from defects
	//Otherwise, it should be set to "NULL"
	string inputlargetris;

	//Input file for defects positions (This should be given if 'inputtriscon' is given)
	//If this file is given, then the code do not search for defects or build large triangles from defects
	//Otherwise, it should be set to "NULL"
	string inputdef;

	//Input file for smaller triangles connectivity
	//Should be set to "NULL" if not given
	string inputsmalltris;

	//Folder name where the output is stored
	string folder;

	bool Minim = 1;				// 1 if minimization is on, 0 if not
	int mststeps = 7500;			// Number of times MST algorithm will be run to find unfolding path
	

	//.............Input Specifications...................//
	inputparts = "SPHEREABFFLAT.vtk";
	inputcon = "ABF3Dtris.txt";
	inputlargetris = "NULL";
	inputdef = "NULL";
	inputsmalltris = "NULL";
	folder = "Starbust"; 

	/*inputparts = "hexparflat.vtk";
	inputcon = "ABF3Dtris.txt";
	inputlargetris = "NULL";
	inputdef = "NULL";
	inputsmalltris = "NULL";
	folder = "Starbust";*/

	//Search Parameters for Defects Connectivity
	Real RC = 4.0, DELR = 2.0, DR = 1.0, INITIALR = 1.0;
	ParticlePositions Part(numberof_particles, inputparts, inputcon, inputlargetris, inputsmalltris, folder);
	
	vector<Vector2d>Positions;
	map<int, int>ord3d;
	int np = Part.n_part;
	
	for (int i = 0; i < np; i++) {
		Positions.push_back({ Part.pos[i][0],Part.pos[i][1] });
		ord3d[i] = i;
	}

	
	map<string, int>mid3d;
	map<vector<int>, bool>donetris;
	string _text;
	ifstream _readFile("SPHEREABFTris.txt");
	int indexnow = 0;
	while (getline(_readFile, _text)) {
		Tris.push_back(vector<int>());
		string word;
		istringstream s(_text);
		while (s >> word)
		{
			Tris[indexnow].push_back(stod(word));
		}
		indexnow++;
	}
	_readFile.close();
	cout << Tris.size() << " Tris " << endl;
	set<string>midptstrs;
	//....................Building the midpoints for quadratic formulation................//
	for (int DTG = 0; DTG < Tris.size(); DTG++) {
		int p1 = Tris[DTG][0];
		int p2 = Tris[DTG][1];
		int p3 = Tris[DTG][2];
		vector<pair<int, int>>ptpairs;
		pair<int, int>pr1, pr2, pr3;
		if (p1 < p2)pr1 = make_pair(p1, p2);
		else pr1 = make_pair(p2, p1);
		if (p2 < p3)pr2 = make_pair(p2, p3);
		else pr2 = make_pair(p3, p2);
		if (p1 < p3)pr3 = make_pair(p1, p3);
		else pr3 = make_pair(p3, p1);

		ptpairs.push_back(pr1);
		ptpairs.push_back(pr2);
		ptpairs.push_back(pr3);
		for (int ptk = 0; ptk < 3; ptk++) {
			int pt1 = ptpairs[ptk].first;
			int pt2 = ptpairs[ptk].second;
			string ptnow = to_string(pt1) + "a" + to_string(pt2);
			if (midptstrs.find(ptnow) == midptstrs.end()) {
				midptstrs.insert(ptnow);
				int pm1 = ord3d[pt1];
				int pm2 = ord3d[pt2];
				Vector2d midpm = (Positions[pm1] + Positions[pm2]) / 2.0;
				Positions.push_back(midpm);
				mid3d[ptnow] = Positions.size() - 1;
			}
		}
	}
	string hex3d = "SPHEREABF3D.vtk";
	Body bd(Positions, Tris, ord3d,mid3d, hex3d, 1, np);
	VectorXd xo;
	xo.resize(Positions.size() * 2);
	vector<double>Pos2;
	Pos2.clear();
	Pos2 = bd.getPos();
	cout << bd.getenrg() << "Energy" << endl;
	for (int l = 0; l < 2 * Positions.size(); l++) {
		xo[l] = Pos2[l];
	}
	ordr = ord3d;
	s = hex3d;
	stordr = mid3d;
	nprt = np;
	
	//.................LBFGS details...............//
	LBFGSBParam<double> param;
	param.epsilon = 1e-8;
	param.delta = 1e-8;
	param.epsilon_rel = 1e-7;
	LBFGSBSolver<double> solver(param);
	double fx;
	iter = 0;
	vector<int>fix = { 0 };
	VectorXd lb = VectorXd::Constant(2 * Positions.size(), -1000.0);
	VectorXd ub = VectorXd::Constant(2 * Positions.size(),  1000.0);
	for (int i = 0; i < fix.size(); i++) {
		lb[2 * fix[i] + 0] = Positions[fix[i]][0];
		ub[2 * fix[i] + 0] = Positions[fix[i]][0];
		lb[2 * fix[i] + 1] = Positions[fix[i]][1];
		ub[2 * fix[i] + 1] = Positions[fix[i]][1];
	}
	int niter = solver.minimize(foo, xo, fx, lb, ub);
	EnergyF.push_back(fx);
	for (int l = 0; l < Positions.size(); l++) {
		Positions[l][0] = xo[2 * l + 0];
		Positions[l][1] = xo[2 * l + 1];
	}
	Body bd2(Positions, Tris, ord3d,mid3d, hex3d, 1, np);
	vector< pair<vector<int>, Real>> difs;
	difs = bd2.getdiff();

	//.........................File writing.......................//
	ofstream f;
	f.open(folder+"/Minimized" + Part.input_particlefile);
	f << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS ";
	f << Positions.size();
	f << " float\n";
	for (int i = 0; i < Positions.size(); i++) {
		f << Positions[i][0] << " " << Positions[i][1] << " " << 0 << endl;
	}
	f << endl;
	f << "CELLS ";
	f << Tris.size() << " " << 7 * Tris.size() << endl;
	for (int DTG = 0; DTG < Tris.size(); DTG++) {
		vector<int>pms;
		int p1 = Tris[DTG][0];
		int p2 = Tris[DTG][1];
		int p3 = Tris[DTG][2];
		vector<pair<int, int>>ptpairs;
		pair<int, int>pr1, pr2, pr3;
		if (p1 < p2)pr1 = make_pair(p1, p2);
		else pr1 = make_pair(p2, p1);
		if (p2 < p3)pr2 = make_pair(p2, p3);
		else pr2 = make_pair(p3, p2);
		if (p1 < p3)pr3 = make_pair(p1, p3);
		else pr3 = make_pair(p3, p1);

		ptpairs.push_back(pr1);
		ptpairs.push_back(pr2);
		ptpairs.push_back(pr3);
		for (int ptk = 0; ptk < 3; ptk++) {
			int pt1 = ptpairs[ptk].first;
			int pt2 = ptpairs[ptk].second;
			string ptnow = to_string(pt1) + "a" + to_string(pt2);
			pms.push_back(mid3d[ptnow]);
		}
		f << 6 << " " << Tris[DTG][0] << " " << Tris[DTG][1] << " " << Tris[DTG][2] << " " << pms[0] << " " << pms[1] << " " << pms[2] << endl;
	}
	f << endl;
	f << "CELL_TYPES ";
	f << Tris.size() << endl;
	for (int i = 0; i < Tris.size(); i++) {
		//f << 5 << endl;
		f << 22 << endl;
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
				//cout << difs[j].second << endl;
				found = true;
				break;
			}
		}
	}
	f.close();
	return 0;
}