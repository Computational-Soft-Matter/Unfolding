#include "UnfoldTriangles.h"
#include <LBFGSB.h>
namespace unfolding {
	// Parameters for LBFGS minimization
	Real param_epsilon = 1e-6;
	Real param_epsilon_rel = 1e-5;
	Real param_delta = 1e-6;
	//Real param_m = 10;
	//Real param_max_linesearch = 40;
	//Real param_ftol = 1e-5;
	//Real param_wolfe = 0.5;
	//Real param_max_submin = 20;
	//Real param_max_iterations = 500;

	int stopminim = 10000;          // Run minimization if number of smaller triangles is less than stopminim (initially set at value greater than number of total smaller triangles)
	int iter = 1;					// Number of iterations to keep track of solver
	bool shouldprint = 0;			// If 1, then print energy and projected gradients
	int startprint = 279;			// Iteration number from when to start printing energy and projected gradients
									// To print on just the last step: T7 = 419, T3 = 179, Parabola=279, Dome=334

	// Parameters to send to 'Body.h/cpp' to build model
	int nprt;						
	string s;						
	vector<vector<int>>Tris;
	vector<Real>EnergyF;
	vector<Real>PGrad;
	vector<pair<int, Vector3d>>midplots;
	map<int, int>ordr;
	map<string, int>MIDP;	

	// Minimization functions to compute energy, gradients
	pair<double,VectorXd> computes(VectorXd x) {
		vector<Vector2d>Pos;
		int n = (x.size());
		Pos.resize(int(n / 2), { 0.0,0.0 });
		for (int i = 0; i<int(n / 2); i++) {
			Pos[i][0] = x[2 * i + 0];
			Pos[i][1] = x[2 * i + 1];
		}
		Body bdcomp(Pos, Tris, ordr,MIDP, s, 0, nprt);
		double enrg = bdcomp.getenrg();
		vector<double>gr = bdcomp.getGrd();
		VectorXd GR;
		GR.resize(n);
		for (int i = 0; i < n; i++) {
			GR[i] = gr[i];
		}
		return make_pair(enrg,GR);
	}
	double foo(const VectorXd& x, VectorXd& grad)
	{
		iter++;
		pair<double, VectorXd>computeds=computes(x);
		double f,maxdf=0.0;
		f = computeds.first;
		if (shouldprint) {
			cout << iter << "*************************************** " << f << endl;
			VectorXd LocalGrad = computeds.second;
			for (int i = 0; i < LocalGrad.size(); i++) {
				if (fabs(LocalGrad[i]) > maxdf)maxdf = fabs(LocalGrad[i]);
			}
			cout << "Infinity Norm of GR: " << maxdf << endl;
		}
		grad = computeds.second;
		return f;
	}

	void UnfoldTriangles::PrintEnergy() {
		ofstream f;
		f.open(_partpos._folder + "/Energy" + _partpos.input_particlefile);
		for (int i = 0; i < EnergyF.size(); i++) {
			f << EnergyF[i] << endl;
		}
		f.close();
		ofstream fg;
		fg.open(_partpos._folder + "/PGrad" + _partpos.input_particlefile);
		for (int i = 0; i < PGrad.size(); i++) {
			fg << PGrad[i] << endl;
		}
		fg.close();
	}
	void UnfoldTriangles::PrintPts(){
		ofstream f;
		f.open("Unfolded"+_partpos.input_particlefile);
		f << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS ";
		f << pos2D.size();
		f << " float\n";
		for (int i = 0; i < pos2D.size(); i++) {
			f << pos2D[i][0] << " " << pos2D[i][1] << " " << 0 << endl;
		}
		f << endl;
		f.close();
	}
	void UnfoldTriangles::PrintTrs() {
		int totaltris = 0;
		vector<int>notprint;
		int problemco = 0;
		for (int i = 0; i < donetrisglobal.size(); i++) {
			if (vector_search(notprint, i) == 1) {
				continue;
			}
			vector<int>maint = { 0,0,0 };
			for (int j = 0; j < 3; j++) {
				maint[j] = donetrisglobal[i][j] % _partpos.n_part;
			}
			sort(maint.begin(), maint.end());
			bool problem=0;
			for (int j = 0; j < donetrisglobal.size(); j++) {
				if (i != j) {
					vector<int>maing = { 0,0,0 };
					for (int r = 0; r < 3; r++) {
						maing[r] = donetrisglobal[j][r] % _partpos.n_part;
					}
					if (maint[0] == maing[0] && maint[1] == maing[1] && maint[2] == maing[2]) {
						problem = 1;
						problemco++;
						notprint.push_back(j);
					}
				}
			}
			totaltris++;
		}
		ofstream f;
		f.open("UnfoldedTris" + _partpos.input_particlefile);
		f << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS ";
		f << pos2D.size();
		f << " float\n";
		for (int i = 0; i < pos2D.size(); i++) {
			f << pos2D[i][0] << " " << pos2D[i][1] << " " << 0 << endl;
		}
		f << endl;
		f << "CELLS ";
		f << totaltris << " " << 4 * totaltris << endl;
		for (int i = 0; i < donetrisglobal.size(); i++) {
			if (vector_search(notprint, i) != 1) {
				f << 3 << " " << ID2D[donetrisglobal[i][0]] << " " << ID2D[donetrisglobal[i][1]] << " " << ID2D[donetrisglobal[i][2]] << endl;
			}
		}
		f << endl;
		f << "CELL_TYPES ";
		f << totaltris << endl;
		for (int i = 0; i < totaltris; i++) {
			f << 5 << endl;
		}
		f.close();
	}
	void UnfoldTriangles::PrintOutline() {
		map<int, int>donepts;
		vector<int>allpts;
		vector<pair<int, int>>edges;
		map<int, int>IDnow;
		int id = 0;
		for (int i = 0; i <_unfpt.plottedfinalpts.size(); i++) {
			vector<int>now = _unfpt.plottedfinalpts[i];
			for (int j = 0; j < 3; j++) {
				if (donepts[now[j]] != -5) {
					donepts[now[j]] = -5;
					IDnow[now[j]] = id;
					id++;
					allpts.push_back(now[j]);
				}
			}
			edges.push_back(make_pair(now[0], now[1]));
			edges.push_back(make_pair(now[0], now[2]));
			edges.push_back(make_pair(now[2], now[1]));
		}
		ofstream f;
		f.open("UnfoldedLines" + _partpos.input_particlefile);
		f << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS ";
		f << allpts.size();
		f << " float\n";
		for (int i = 0; i < allpts.size(); i++) {
			f << pos2D[ID2D[allpts[i]]][0] << " " << pos2D[ID2D[allpts[i]]][1] << " " << 0 << endl;
		}
		f << endl;
		f << "LINES ";
		f << edges.size() << " " << 3 * edges.size() << endl;
		for (int i = 0; i < edges.size(); i++) {
			f << 2 << " " <<IDnow[edges[i].first] << " " << IDnow[edges[i].second] <<  endl;
		}
		f.close();
	}

	void UnfoldTriangles::find_paths(vector<vector<int> >& pathsn,
		vector<int>& pathn,
		vector<vector<int>>& parentn,
		int nn, int un, int sm, int tm)
	{
		if (un == -1) {
			pathsn.push_back(pathn);
			return;
		}
		for (int par : parentn[un]) {
			pathn.push_back(un);
			find_paths(pathsn, pathn, parentn,
				nn, par, sm, tm);
			pathn.pop_back();
		}
	}

	TriPtsDS UnfoldTriangles::TriPts(vector<int>T) {
		vector<int>final_vert;
		vector<vector<int>>trails;
		vector<Vector3d>plots;
		vector<int>initial_vert;
		queue<int>path;
		bool snglcut = 0;
		int n1 = T[0], n2 = T[1], n3 = T[2];
		vector<pair<int, int>>tr_edge;
		pair<int, int>path1 = make_pair(n1, n2);
		pair<int, int>path2 = make_pair(n1, n3);
		pair<int, int>path3 = make_pair(n2, n3);
		tr_edge.push_back(path1);
		tr_edge.push_back(path2);
		tr_edge.push_back(path3);
		vector<vector<int>>pathsnow;
		vector<bool>visitedpath;
		// Running BFS/ Finding previously discovered shortest path between two nodes in the graph
		for (int j = 0; j < 3; j++) {
			visitedpath.clear();
			int source, target;
			visitedpath.resize(_partpos.n_part, 0);
			bool refinement = 0;
			if (refinement) {
				if (tr_edge[j].first < tr_edge[j].second) {
					source = tr_edge[j].first;
					target = tr_edge[j].second;
				}
				else {
					source = tr_edge[j].first;
					target = tr_edge[j].second;
				}
			}
			else {
				source = tr_edge[j].first;
				target = tr_edge[j].second;
			}
			int smain = source;
			int tmain = target;
			int othernode;
			if (j == 0)othernode = n3;
			else if (j == 1)othernode = n2;
			else othernode = n1;
			vector<int>possiblepath;
			vector<vector<int>>parents;
			for (int i = 0; i < _partpos.n_part; i++) {
				parents.push_back(vector<int>());
			}
			parents[source] = { -1 };
			vector<int>graphdis(_partpos.n_part, 99999);
			graphdis[source] = 0;
			for (int k = 0; k < allpaths.size(); k++) {
				if ((allpaths[k].first.first == source && allpaths[k].first.second == target) || (allpaths[k].first.first == source && allpaths[k].first.second == target)) {
					possiblepath = allpaths[k].second;
					break;
				}
			}
			vector<int>trail;
			if (possiblepath.size() > 0) {				// If shortest path previously found then this part runs
				for (int k = 0; k < possiblepath.size(); k++) {
					if (visitedpath[possiblepath[k]] == false) {
						final_vert.push_back(possiblepath[k]);
						initial_vert.push_back(possiblepath[k]);
						path.push(possiblepath[k]);
						visitedpath[possiblepath[k]] = true;
					}
				}
				trail = possiblepath;
				pathsnow.push_back(possiblepath);
			}
			else {									// If shortest path previously not found then this part runs BFS algorithm
				map<int, int>pathfinder;
				pathfinder[source] = -1;
				queue<int>bfs;
				vector<bool>visited;
				visited.resize(_partpos.n_part, false);
				visited[source] = true;
				bfs.push(source);
				bool stop = false;
				while (!bfs.empty() && !stop) {
					source = bfs.front();
					bfs.pop();
					int sourcesize = _partpos.pos_neighbor[source].size();
					for (int k = 0; k < sourcesize; k++) {
						int adj = _partpos.pos_neighbor[source][k];
						if (adj >= _partpos.n_part)continue;
						if (graphdis[adj] > graphdis[source] + 1) {
							graphdis[adj] = graphdis[source] + 1;
							bfs.push(adj);
							parents[adj].clear();
							parents[adj].push_back(source);
						}
						else if (graphdis[adj] == graphdis[source] + 1) {
							parents[adj].push_back(source);
						}
					}
				}
				vector<int>targetparent = { target };
				vector<vector<int>>mulpaths;
				vector<int>singlepath;
				find_paths(mulpaths, singlepath, parents, _partpos.n_part, target, smain, tmain);
				int biggestdist = -1;
				int pathdorkar = 0;
				Real mindis = 9999.0;
				for (int k = 0; k < mulpaths.size(); k++) {
					Vector3d B = _partpos.pos[mulpaths[k][0]];
					Vector3d C = _partpos.pos[mulpaths[k][mulpaths[k].size()-1]];
					Real totdis = 0.0;
					for (int kk = 0; kk < mulpaths[k].size(); kk++) {
						Vector3d A = _partpos.pos[mulpaths[k][kk]];
						Vector3d BC = B - C;
						Vector3d AC = A - C;
						Vector3d crossProduct = BC.cross(AC);
						totdis += (crossProduct.norm() / BC.norm());
					}
					if (totdis < mindis) {
						mindis = totdis;
						pathdorkar = k;
					}
				}
				int pathc = 0;
				target = mulpaths[pathdorkar][pathc];
				pair<int, int>pathmain;
				pathmain = make_pair(tr_edge[j].first, tr_edge[j].second);
				while (pathc != mulpaths[pathdorkar].size()) {
					target = mulpaths[pathdorkar][pathc];
					if (visitedpath[target] == false) {
						path.push(target);
						if (vector_search(final_vert, target) == 0)
							final_vert.push_back(target);
						initial_vert.push_back(target);
						visitedpath[target] = true;
						trail.push_back(target);
					}
					pathc++;
					target = pathfinder[target];
				}
				allpaths.push_back(make_pair(tr_edge[j], trail));
			}
			trails.push_back(trail);
		}
		// This section goes through the shortest path trails to find a set of particles adjacent to the particles in the trail and inside the large triangle
		for (int rr = 0; rr < trails.size(); rr++) {
			int tester;
			for (int r = 0; r < 3; r++) {
				if (T[r] != tr_edge[rr].first && T[r] != tr_edge[rr].second) {
					tester = T[r];
				}
			}
			for (int r = 0; r < trails[rr].size(); r++) {
				int trailnow = trails[rr][r];
				if (tester == 0) {
					if (vector_search(_unfpt._tri._defectpos.defect_pos, trailnow) == 1)continue;
				}
				else {
					if (vector_search(T, trailnow) == 1)continue;
				}
				int vnsize = _partpos.pos_neighbor[trailnow].size();
				vector<vector<int>>donesets;
				donesets.resize(2, vector<int>());
				int setnow = 0;
				int posnow = 0;
				int totalfns = 0;
				for (int rrr = 0; rrr < vnsize; rrr++) {
					int trailneig = _partpos.pos_neighbor[trailnow][rrr];
					for (int rrk = 0; rrk < trails.size(); rrk++) {
						if (vector_search(trails[rrk], trailneig) == 1)
							totalfns++;
					}
					if (vector_search(trails[rr], trailneig) == 1) {
						posnow = 1 - posnow;
					}
					else {
						donesets[posnow].push_back(trailneig);
					}
				}
				if (totalfns >= 3 && tester != 0)continue;
				vector<Real>distdone;
				for (int l = 0; l < 2; l++) {
					distdone.push_back(0.0);
					if (donesets[l].size() == 0) {
						distdone[l] = 100;
						continue;
					}
					for (int ll = 0; ll < donesets[l].size(); ll++) {
						if (donesets[l][ll] > _partpos.n_part)continue;
						distdone[l] += (Real)_partpos.MinDisGrp[tester][donesets[l][ll]];
					}
					distdone[l] = distdone[l] / (Real)(donesets[l].size());
				}


				/********Only for antenna*/
				if (tester == 0 && snglcut) {
					for (int ll = 0; ll < donesets[0].size(); ll++) {
						if (vector_search(final_vert, donesets[0][ll]) == 0 && donesets[0][ll] < _partpos.n_part) {
							path.push(donesets[0][ll]);
							final_vert.push_back(donesets[0][ll]);
						}
					}
					for (int ll = 0; ll < donesets[1].size(); ll++) {
						if (vector_search(final_vert, donesets[1][ll]) == 0 && donesets[1][ll] < _partpos.n_part) {
							path.push(donesets[1][ll]);
							final_vert.push_back(donesets[1][ll]);
						}
					}
				}
				else {
					if (distdone[0] < distdone[1]) {
						for (int ll = 0; ll < donesets[0].size(); ll++) {
							if (vector_search(final_vert, donesets[0][ll]) == 0 && donesets[0][ll] < _partpos.n_part) {
								path.push(donesets[0][ll]);
								final_vert.push_back(donesets[0][ll]);
							}
						}
					}
					else {
						for (int ll = 0; ll < donesets[1].size(); ll++) {
							if (vector_search(final_vert, donesets[1][ll]) == 0 && donesets[1][ll] < _partpos.n_part) {
								path.push(donesets[1][ll]);
								final_vert.push_back(donesets[1][ll]);
							}
						}
					}
				}
			}
		}
		Real At = trianglearea(_partpos.pos[n1], _partpos.pos[n2], _partpos.pos[n3]);
		Vector3d v12(_partpos.pos[n2][0] - _partpos.pos[n1][0], _partpos.pos[n2][1] - _partpos.pos[n1][1], _partpos.pos[n2][2] - _partpos.pos[n1][2]);
		Vector3d v13(_partpos.pos[n3][0] - _partpos.pos[n1][0], _partpos.pos[n3][1] - _partpos.pos[n1][1], _partpos.pos[n3][2] - _partpos.pos[n1][2]);
		Vector3d nml = v12.cross(v13);
		Real norm = nml.norm();
		Vector3d n = nml / norm;
		Real d = _partpos.pos[n1][0];
		Real e = _partpos.pos[n1][1];
		Real f = _partpos.pos[n1][2];
		Real pa = (_partpos.pos[n2][1] - _partpos.pos[n1][1]) * (_partpos.pos[n3][2] - _partpos.pos[n1][2]) - (_partpos.pos[n3][1] - _partpos.pos[n1][1]) * (_partpos.pos[n2][2] - _partpos.pos[n1][2]);
		Real pb = (_partpos.pos[n2][2] - _partpos.pos[n1][2]) * (_partpos.pos[n3][0] - _partpos.pos[n1][0]) - (_partpos.pos[n3][0] - _partpos.pos[n1][0]) * (_partpos.pos[n2][1] - _partpos.pos[n1][1]);
		Real pc = (_partpos.pos[n2][0] - _partpos.pos[n1][0]) * (_partpos.pos[n3][1] - _partpos.pos[n1][1]) - (_partpos.pos[n3][0] - _partpos.pos[n1][0]) * (_partpos.pos[n2][1] - _partpos.pos[n1][1]);
		Real pd = -(pa * _partpos.pos[n1][0] + pb * _partpos.pos[n1][1] + pc * _partpos.pos[n1][2]);
		for (int j = 0; j < final_vert.size(); j++) {
			int vn_ngh = final_vert[j];
			if (vn_ngh > _partpos.n_part)continue;
			Real x = _partpos.pos[vn_ngh][0];
			Real y = _partpos.pos[vn_ngh][1];
			Real z = _partpos.pos[vn_ngh][2];
			Real t = (n(0) * d - n(0) * x + n(1) * e - n(1) * y + n(2) * f - n(2) * z) / (n(0) * n(0) + n(1) * n(1) + n(2) * n(2));
			Vector3d proj = { x + t * n(0),y + t * n(1),z + t * n(2) };
			plots.push_back(_partpos.pos[vn_ngh]);
		}
		vector<int>testvec;
		vector<Vector3d>plotsbackup;
		visitedpath.clear();
		visitedpath.resize(_partpos.n_part, false);
		visitedpath[n1] = 1;
		visitedpath[n2] = 1;
		visitedpath[n3] = 1;
		//This section finds all the particles that are inside the large triangle
		while (!path.empty()) {
			int verticenow = path.front();
			path.pop();
			if (verticenow > _partpos.n_part)continue;
			int vnsize = _partpos.pos_neighbor[verticenow].size();
			for (int k = 0; k < vnsize; k++) {
				int vn_ngh = _partpos.pos_neighbor[verticenow][k];
				if (vn_ngh >= _partpos.n_part)continue;
				if (visitedpath[vn_ngh] == false) {
					if (verticenow != n1 && verticenow != n2 && verticenow != n3)visitedpath[vn_ngh] = true;
					Real x = _partpos.pos[vn_ngh][0];
					Real y = _partpos.pos[vn_ngh][1];
					Real z = _partpos.pos[vn_ngh][2];
					Real t = (n(0) * d - n(0) * x + n(1) * e - n(1) * y + n(2) * f - n(2) * z) / (n(0) * n(0) + n(1) * n(1) + n(2) * n(2));
					Vector3d proj = { x + t * n(0),y + t * n(1),z + t * n(2) };
					Real At_s = trianglearea(proj, _partpos.pos[n1], _partpos.pos[n2]) + trianglearea(proj, _partpos.pos[n3], _partpos.pos[n2]) + trianglearea(proj, _partpos.pos[n1], _partpos.pos[n3]);
					float a = ceil(At_s * 100.0) / 100.0;
					float b = ceil(At * 100.0) / 100.0;
					if (a == b) {
						visitedpath[vn_ngh] = true;
						path.push(vn_ngh);
						if (vector_search(final_vert, vn_ngh) == 0)final_vert.push_back(vn_ngh);
						plots.push_back(_partpos.pos[vn_ngh]);
					}
				}
			}
		}
		TriPtsDS tempDS;
		tempDS.Vrt = final_vert;
		tempDS.Plts = plots;
		return tempDS;

	}
	TriPtsDS UnfoldTriangles::ScarTriPts(vector<int>T,vector<vector<int>>Trls) {
		vector<int>final_vert;
		vector<vector<int>>trails=Trls;
		vector<int>donttouch;
		for (int i = 0; i < Trls.size(); i++) {
			for (int j = 0; j < Trls[i].size(); j++) {
				donttouch.push_back(Trls[i][j]);
			}
		}
		vector<Vector3d>plots;
		queue<int>path;
		int tester;
		int n1 = T[0], n2 = T[1], n3 = T[2];
		vector<pair<int, int>>tr_edge;
		pair<int, int>path1 = make_pair(n1, n2);
		pair<int, int>path2 = make_pair(n1, n3);
		pair<int, int>path3 = make_pair(n2, n3);
		tr_edge.push_back(path1);
		tr_edge.push_back(path2);
		tr_edge.push_back(path3);
		int thistrail;
		for (int rr = 0; rr < trails.size(); rr++) {
			if (vector_search(trails[rr], n1) == 1 && vector_search(trails[rr], n2) == 1)thistrail = rr;
		}
		for (int r = 0; r < 3; r++) {
			if (T[r] != n1 && T[r] != n2) {
				tester = T[r];
			}
		}
		vector<bool>visitedpath;
		visitedpath.clear();
		visitedpath.resize(_partpos.n_part, false);
		for (int r = 0; r < trails[thistrail].size(); r++) {
			int trailnow = trails[thistrail][r];
			final_vert.push_back(trailnow);
			visitedpath[trailnow] = 1;
			int vnsize = _partpos.pos_neighbor[trailnow].size();
			vector<vector<int>>donesets;
			donesets.resize(2, vector<int>());
			int setnow = 0;
			int posnow = 0;
			int totalfns = 0;
			bool isproblem = 0;
			for (int rrr = 0; rrr < vnsize; rrr++) {
				int trailneig = _partpos.pos_neighbor[trailnow][rrr];
				for (int rrk = 0; rrk < trails.size(); rrk++) {
					if (vector_search(trails[rrk], trailneig) == 1)
						totalfns++;
				}
				if (vector_search(trails[thistrail], trailneig) == 1) {
					if (isproblem == 0) {
						isproblem = 1;
						posnow = 1 - posnow;
					}
				}
				else {
					donesets[posnow].push_back(trailneig);
					isproblem = 0;
				}
			}
			if (totalfns >= 3 || trailnow==28 || trailnow == 143 || trailnow == 142 || trailnow == 141 || trailnow == 144)continue;
			vector<Real>distdone;
			for (int l = 0; l < 2; l++) {
				distdone.push_back(0.0);
				if (donesets[l].size() == 0) {
					distdone[l] = 100;
					continue;
				}
				for (int ll = 0; ll < donesets[l].size(); ll++) {
					if (donesets[l][ll] > _partpos.n_part)continue;
					distdone[l] += (Real)_partpos.MinDisGrp[tester][donesets[l][ll]];
				}
				distdone[l] = distdone[l] / (Real)(donesets[l].size());
			}
			if (distdone[0] < distdone[1]) {
				for (int ll = 0; ll < donesets[0].size(); ll++) {
					if (vector_search(final_vert, donesets[0][ll]) == 0 && donesets[0][ll] < _partpos.n_part && vector_search(donttouch, donesets[0][ll])==0) {
						path.push(donesets[0][ll]);
						final_vert.push_back(donesets[0][ll]);
						visitedpath[donesets[0][ll]] = 1;
					}
				}
			}
			else {
				for (int ll = 0; ll < donesets[1].size(); ll++) {
					if (vector_search(final_vert, donesets[1][ll]) == 0 && donesets[1][ll] < _partpos.n_part && vector_search(donttouch, donesets[1][ll]) == 0) {
						path.push(donesets[1][ll]);
						final_vert.push_back(donesets[1][ll]);
						visitedpath[donesets[1][ll]] = 1;
					}
				}
			}	
		}
		visitedpath[n1] = 1;
		visitedpath[n2] = 1;
		while (!path.empty()) {
			int verticenow = path.front();
			path.pop();
			if (verticenow > _partpos.n_part)continue;
			int vnsize = _partpos.pos_neighbor[verticenow].size();
			for (int k = 0; k < vnsize; k++) {
				int vn_ngh = _partpos.pos_neighbor[verticenow][k];
				if (vn_ngh >= _partpos.n_part)continue;
				if (visitedpath[vn_ngh] == false) {
					visitedpath[vn_ngh] = true;
					if(vector_search(donttouch,vn_ngh)==0)path.push(vn_ngh);
					if (vector_search(final_vert, vn_ngh) == 0)final_vert.push_back(vn_ngh);
				}
			}
		}
		for (int j = 0; j < final_vert.size(); j++) {
			int vn_ngh = final_vert[j];
			if (vn_ngh > _partpos.n_part)continue;
			plots.push_back(_partpos.pos[vn_ngh]);
		}
		TriPtsDS tempDS;
		tempDS.Vrt = final_vert;
		tempDS.Plts = plots;
		return tempDS;
	}
	void UnfoldTriangles::Unfold(bool isScar) {
		int totaltris = 0;
		int tricounter = 0;
		vector<vector<int>>TrianglesAll;
		map<vector<int>, bool>isTriangle;
		if (_partpos.input_triangles != "NULL") {
			TrianglesAll = _partpos.Triangles;
			for (int i = 0; i < TrianglesAll.size(); i++) {
				vector<int>p = TrianglesAll[i];
				sort(p.begin(), p.end());
				isTriangle[p] = true;
			}
		}
		map<vector<int>,bool>donetrisall;
		vector<pair<int, int>>doneedge;
		ofstream f3f;
		ofstream f2f;
		f3f.open(_partpos._folder+"/extras" + _partpos.input_connectivityfile);
		f2f.open(_partpos._folder+"/verts" + _partpos.input_connectivityfile);
		vector<int>alldone;
		for (int i = 0; i <_unfpt.t_size; i++) {
			vector<int>previouscompare;
			vector<int>unfnow = _unfpt.plottedfinalpts[i];
			vector<int>unfnoworg = { unfnow[0] % _partpos.n_part, unfnow[1] % _partpos.n_part, unfnow[2] % _partpos.n_part };
			TriPtsDS DSnow;
			if (!isScar)DSnow = TriPts(unfnoworg);
			else DSnow = ScarTriPts(unfnoworg, trals);
			vector<int>final_vert = DSnow.Vrt;
			sort(final_vert.begin(),final_vert.end());
			print_vector(final_vert);
			vector<Vector3d>plots = DSnow.Plts;
			int o1 = unfnow[0];
			int o2 = unfnow[1];
			int o3 = unfnow[2];
			//cout << i << endl;
			cout << o1 << " mps " << o2 << " mps " << o3 << endl;
			//cout << "final_vert: " << endl;
			nountripts.push_back(vector<int>());
			// Finding all the small triangles inside the triangular region bounded by o1,o2,o3
			vector<vector<int>>smalltris;
			for (int j = 0; j < final_vert.size(); j++) {
				//vector<int>badlist = { 183,184,185,186,187,188,189,190,191,192,2,13,28,137,153,154,156,157,158,23,24,25,26,27 };  
				vector<int>badlist = { 0,70,76,77,78,79,80,106,109,316,319,336,340,342,347,69,71,72,73,74,75,312,315,323,327,331,335,37,55,58,64,65,68,232,283,287,298,308,6,28,29,30,31,32,129,131,206,210,212,217,221,1,2,3,4,5,113,117,121,125 };
				//vector<int>badlist = { 1,2,3,4,5,6,7,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,171,178,180,184,188,192,196,198,282,284,291,293,300,302,309,311,313,315,317,321,325,329,333,337,341,345,349,353,357,361,365,369,373,377,381 };
				badlist.clear();
				int indn = final_vert[j];
				int sizenow = _partpos.pos_neighbor[indn].size();
				//cout << indn << "now" << endl;
				if (_partpos.input_smalltris != "NULL"  ) {
					for (int k = 0; k < sizenow; k++) {
						int nnow1 = _partpos.pos_neighbor[indn][k];
						for (int kk = 0; kk < sizenow; kk++) {
							int nnext1;
							if (nnow1 != _partpos.pos_neighbor[indn][kk]) {
								nnext1 = _partpos.pos_neighbor[indn][kk];
							}
							else continue;
							if ((vector_search(final_vert, nnow1) == 1) && (vector_search(final_vert, nnext1) == 1)) {
								vector<int>P;
								P.push_back(indn);
								P.push_back(nnow1);
								P.push_back(nnext1);
								sort(P.begin(), P.end());
								bool pushnow = 0;
								if ((vector_search(badlist, nnow1) == 1) && (vector_search(badlist, nnext1) == 1) && (vector_search(badlist, indn) == 1))pushnow = 1;
								for (int l = 0; l < smalltris.size(); l++) {
									if (smalltris[l] == P)pushnow = 1;
								}
								if (pushnow == 0) {
									//print_vector(P);
									if(isTriangle[P]==true)smalltris.push_back(P);
									totaltris++;
								}
							}
						}
					}
				}
				else {
					int sizefinal = sizenow;
					if (sizenow <= 4) {
						sizefinal = sizenow - 1;;
					}
					for (int k = 0; k < sizefinal; k++) {
						int nnow = _partpos.pos_neighbor[indn][k];
						int nnext = _partpos.pos_neighbor[indn][(k + 1) % sizenow];
						if ((vector_search(final_vert, nnow) == 1) && (vector_search(final_vert, nnext) == 1)) {
							vector<int>P;
							P.push_back(indn);
							P.push_back(nnow);
							P.push_back(nnext);
							Real l1 = dist(_partpos.pos[indn], _partpos.pos[nnow]);
							Real l2 = dist(_partpos.pos[indn], _partpos.pos[nnext]);
							Real l3 = dist(_partpos.pos[nnow], _partpos.pos[nnext]);
							Real AR = max(max(l1, l2), l3) / min(min(l1, l2), l3);
							
							sort(P.begin(), P.end());
							bool pushnow = 0;
							if (AR > 2.0)pushnow = 1;
							if ((vector_search(badlist, nnow) == 1) && (vector_search(badlist, nnext) == 1) && (vector_search(badlist, indn) == 1))pushnow = 1;
							for (int l = 0; l < smalltris.size(); l++) {
								if (smalltris[l] == P)pushnow = 1;
							}
							if (pushnow == 0) {
								print_vector(P);
								smalltris.push_back(P);
								totaltris++;
							}
						}
					}
				}
			}
			// Finding which particles to plot first based on previously plotted particles form earlier flattened large triangles
			queue<vector<int>>donetris;
			vector<int>smalldone;
			int firstnode, secondnode;
			int donesmalltris = 0;
			int startvrt = -5;
			vector<int>previousorgpts;
			map<int, int>previouspts;
			if (i == 0)startvrt = o1;
			midplots.clear();
			if (i != 0) {
				for (int j = 0; j < final_vert.size(); j++) {
					previouspts[final_vert[j]] = -5;
					for (int r = 0; r < untripts[_unfpt.finalparenttriangle[i]].size(); r++) {
						if (untripts[_unfpt.finalparenttriangle[i]][r] % _partpos.n_part == final_vert[j]) {
							previousorgpts.push_back(final_vert[j]);
							previouspts[final_vert[j]] = untripts[_unfpt.finalparenttriangle[i]][r];
							if (startvrt == -5 && final_vert[j]!=0) {
								for (int rk = 0; rk < smalltris.size(); rk++) {
									if (vector_search(smalltris[rk], final_vert[j] == 1)) {
										startvrt = final_vert[j];
										break;
									}
								}
							}
							Vector3d placeholder = { pos2D[ID2D[untripts[_unfpt.finalparenttriangle[i]][r]]][0],pos2D[ID2D[untripts[_unfpt.finalparenttriangle[i]][r]]][1],0.0 };
							midplots.push_back(make_pair(untripts[_unfpt.finalparenttriangle[i]][r], placeholder));
							smalldone.push_back(untripts[_unfpt.finalparenttriangle[i]][r] % _partpos.n_part);
						}
					}
				}
			}
			// This section finds the order of smaller triangles to plot sequentially so that for every triangle, two vertices are previously plotted already
			vector<vector<int>>odsmalltris;
			vector<bool>mappedtrism;
			mappedtrism.resize(smalltris.size(), false);
			queue<vector<int>>unftrsm;
			for (int j = 0; j < smalltris.size(); j++) {
				if (i == 0 && vector_search(smalltris[j],startvrt)==1) {
					unftrsm.push(smalltris[j]);
					mappedtrism[j] = true;
					//cout << "DONE" << endl;
					break;
				}
				if(i!=0) {
					int psum = 0;
					bool dont = false;
					for (int r = 0; r < 3; r++) {
						psum+= vector_search(previousorgpts, smalltris[j][r]);
						if (vector_search(nountripts[_unfpt.finalparenttriangle[i]], smalltris[j][r]) == 1)
							dont = true;

					}
					if (psum > 1 && dont==false) {
						unftrsm.push(smalltris[j]);
						mappedtrism[j] = true;
					}
				}
			}
			while (!unftrsm.empty()) {
				vector<int>tnow = unftrsm.front();
				pair<int, int>pair1 = make_pair(tnow[0], tnow[1]);
				pair<int, int>pair2 = make_pair(tnow[0], tnow[2]);
				pair<int, int>pair3 = make_pair(tnow[2], tnow[1]);
				odsmalltris.push_back(tnow);
				unftrsm.pop();
				for (int r = 0; r < smalltris.size(); r++) {
					if (vector_search(smalltris[r], pair1.first) == 1 && vector_search(smalltris[r], pair1.second) == 1 && tnow != smalltris[r]) {
						if (mappedtrism[r] == false) {
							unftrsm.push(smalltris[r]);
							mappedtrism[r] = true;
						}
					}
				}
				for (int r = 0; r < smalltris.size(); r++) {
					if (vector_search(smalltris[r], pair2.first) == 1 && vector_search(smalltris[r], pair2.second) == 1 && tnow != smalltris[r]) {
						if (mappedtrism[r] == false) {
							unftrsm.push(smalltris[r]);
							mappedtrism[r] = true;
						}
					}
				}
				for (int r = 0; r < smalltris.size(); r++) {
					if (vector_search(smalltris[r], pair3.first) == 1 && vector_search(smalltris[r], pair3.second) == 1 && tnow != smalltris[r]) {
						if (mappedtrism[r] == false) {
							unftrsm.push(smalltris[r]);
							mappedtrism[r] = true;
						}
					}
				}
			}
			for (int k = 0; k < 3; k++) {
				if (vector_search(_unfpt._tri._defectpos.defect_pos, odsmalltris[0][k]) == 1) {
					firstnode = odsmalltris[0][k];
					secondnode = odsmalltris[0][(k + 1) % 3];
				}
				else {
					firstnode = odsmalltris[0][k];
					secondnode = odsmalltris[0][(k + 1) % 3];
				}
			}
			if (i != 0) {
				untripts.push_back(vector<int>());
				for (int j = 0; j < previousorgpts.size(); j++) {
					untripts[i].push_back(previouspts[previousorgpts[j]]);
				}
			}
			//This section handles plotting each new particle and subsequently calling the distortion minimization function
			if (i == 0) {				// This condition is for the first large triangle where no previous reference vertices are given
				Real dis = dist(_partpos.pos[firstnode], _partpos.pos[secondnode]);
				Vector3d first = { 0,0,0 };
				Vector3d second = { dis,0,0 };
				pos2D[firstnode][0] = first(0);
				pos2D[firstnode][1] = first(1);
				pos2D[secondnode][0] = second(0);
				pos2D[secondnode][1] = second(1);
				f2f << firstnode << " ";
				f2f << secondnode << " ";
				ord[firstnode] = counter;
				ordrev[counter] = firstnode;
				posord.push_back({ first(0),first(1) });
				counter++;
				ord[secondnode] = counter;
				ordrev[counter] = secondnode;
				posord.push_back({ second(0),second(1) });
				counter++;
				mapped[firstnode] = 1;
				mapped[secondnode] = 1;
				midplots.push_back(make_pair(firstnode, first));
				midplots.push_back(make_pair(secondnode, second));
				smalldone.push_back(firstnode);
				smalldone.push_back(secondnode);
				alldone.push_back(firstnode);
				alldone.push_back(secondnode);
				for (int j = 0; j < odsmalltris.size(); j++) {
					donetrisglobal.push_back(odsmalltris[j]);
					donetrisall[odsmalltris[j]] = 1;
					tricounter++;
					int root, lever, image;
					if (vector_search(smalldone, odsmalltris[j][2]) == 0) {
						if (j == 0) {
							root = firstnode;
							lever = secondnode;
						}
						else {
							root = odsmalltris[j][0];
							lever = odsmalltris[j][1];
						}
						image = odsmalltris[j][2];
					}
					else if (vector_search(smalldone, odsmalltris[j][1]) == 0) {
						if (j == 0) {
							root = firstnode;
							lever = secondnode;
						}
						else {
							root = odsmalltris[j][0];
							lever = odsmalltris[j][2];
						}
						image = odsmalltris[j][1];
					}
					else if (vector_search(smalldone, odsmalltris[j][0]) == 0) {
						if (j == 0) {
							root = firstnode;
							lever = secondnode;
						}
						else {
							root = odsmalltris[j][1];
							lever = odsmalltris[j][2];
						}
						image = odsmalltris[j][0];
					}
					else {
						CalcDistortion(to_string(minimcounter));
						minimcounter++;
						cout << minimcounter << " " << tricounter << " in i=0 " << totaltris << endl;
						if (minimcounter == stopminim)return;
						continue;
					}
					int oldnode=-1;
					for (int k = 0; k < j; k++) {
						if (k != j && vector_search(odsmalltris[k], root) == 1 && vector_search(odsmalltris[k], lever) == 1)
						{
							if (root != odsmalltris[k][0] && lever != odsmalltris[k][0])oldnode = odsmalltris[k][0];
							else if (root != odsmalltris[k][1] && lever != odsmalltris[k][1])oldnode = odsmalltris[k][1];
							else oldnode = odsmalltris[k][2];
						}
					}
					Real dis1 = dist(_partpos.pos[root], _partpos.pos[image]);
					Real dis2 = dist(_partpos.pos[lever], _partpos.pos[image]);
					Real dism = dist(_partpos.pos[root], _partpos.pos[lever]);
					Vector3d rootpos, leverpos, oldpos;
					for (int k = 0; k < midplots.size(); k++) {
						if (root == midplots[k].first) {
							rootpos = midplots[k].second;
						}
						if (lever == midplots[k].first) {
							leverpos = midplots[k].second;
						}
						if (oldnode == midplots[k].first) {
							oldpos = midplots[k].second;
						}
					}
					Real a = (dis1 * dis1 - dis2 * dis2 + dism * dism) / (2 * dism);
					Real h = sqrt(dis1 * dis1 - a * a);
					Vector3d PP0 = { rootpos[0] + a * (leverpos[0] - rootpos[0]) / dism,rootpos[1] + a * (leverpos[1] - rootpos[1]) / dism,0 };
					Vector3d PP1 = { PP0[0] + h * (leverpos[1] - rootpos[1]) / dism,PP0[1] - h * (leverpos[0] - rootpos[0]) / dism ,0 };
					Vector3d PP2 = { PP0[0] - h * (leverpos[1] - rootpos[1]) / dism,PP0[1] + h * (leverpos[0] - rootpos[0]) / dism ,0 };
					Vector3d PP;
					if (j == 0) {
						if (PP1[1] < 0) {
							midplots.push_back(make_pair(image, PP2));
							PP = PP2;
						}
						else {
							midplots.push_back(make_pair(image, PP1));
							PP = PP1;
						}
					}
					else {
						Vector3d aa1, bb1, cc1, First1, Second1;
						aa1 = { leverpos[0] - rootpos[0],leverpos[1] - rootpos[1], 0 };
						bb1 = { oldpos[0] - rootpos[0],oldpos[1] - rootpos[1], 0 };
						cc1 = { PP1[0] - rootpos[0],PP1[1] - rootpos[1], 0 };
						Second1 = aa1.cross(bb1) / aa1.cross(bb1).norm();
						First1 = aa1.cross(cc1) / aa1.cross(cc1).norm();
						
						if (First1[2] == Second1[2]) {
							midplots.push_back(make_pair(image, PP2));
							PP = PP2;
						}
						else {
							midplots.push_back(make_pair(image, PP1));
							PP = PP1;
						}
					}
					ord[image] = counter;
					ordrev[counter] = image;
					posord.push_back({ PP(0),PP(1) });
					counter++;
					f2f << image << " ";
					smalldone.push_back(image);
					alldone.push_back(image);
					pos2D[image][0] = PP(0);
					pos2D[image][1] = PP(1);
					mapped[image] = 1;
					CalcDistortion(to_string(minimcounter));		// Calls distortion minimization function
					minimcounter++;
					cout << minimcounter << " " << tricounter << " in i=0 others " << totaltris << endl;
					if(minimcounter==stopminim)return;
				}
				untripts.push_back(final_vert);
			}
			else {			// This condition is for the next large triangles where previous reference vertices are given
				for (int j = 0; j < odsmalltris.size(); j++) {
					if (minimcounter >= startprint) {
						shouldprint = 1;			// If this condition is reached, then energy values are printed
					}
					int root, lever, image;
					image = -5;
					vector<int>nowtri = odsmalltris[j];
					bool donotprint = 0;
					int rchange;
					if (donetrisall[odsmalltris[j]] == 1)continue;
					for (int r = 0; r < 3; r++) {
						if (vector_search(previousorgpts, nowtri[r]) == 1 || vector_search(smalldone, nowtri[r]) == 1) {
							nowtri[r] = previouspts[nowtri[r]];
						}
						else {
							image = nowtri[r];
							rchange = r;
							if(previouspts[nowtri[(r + 1) % 3]]==-5)
							root = previouspts[nowtri[(r + 1) % 3]];
							else root = previouspts[nowtri[(r + 1) % 3]%_partpos.n_part];
							if (previouspts[nowtri[(r + 2) % 3]] == -5)
							lever = previouspts[nowtri[(r + 2) % 3]];
							else lever = previouspts[nowtri[(r + 2) % 3] % _partpos.n_part];
						}
					}
					// This section decides if the current particle is previously plotted already and no new plotting is necessary
					if (image == -5) {
						for (int r = 0; r < 3; r++) {
							smalldone.push_back(nowtri[r] % _partpos.n_part);
						}	
						donetrisglobal.push_back(nowtri);
						donetrisall[odsmalltris[j]] = 1;
						tricounter++;
						CalcDistortion(to_string(minimcounter));
						minimcounter++;
						cout << minimcounter << " " << tricounter << " in i!=0 not printing any new triangle " << totaltris << endl;
						if (minimcounter == stopminim)return;
						continue;
					}
					if (root == -5 || lever == -5)continue;
					// This section decides if a new instance of a previously plotted particle is necessary to be created
					int indexnow = image;
					if (mapped[image% _partpos.n_part] == 1) {
						int prevsize = rep_pts[indexnow].size();
						indexnow = rep_pts[indexnow][prevsize-1] + _partpos.n_part;
						rep_pts[indexnow% _partpos.n_part].push_back(indexnow);
						image = indexnow;
					}
					previouspts[image%_partpos.n_part] = image;
					nowtri[rchange] = image;
					untripts[i].push_back(image);
					// This section geometrically finds the location of the new particle based on previous two vertices of the current triangle
					int oldnode = -1;
					for (int k = 0; k < donetrisglobal.size(); k++) {
						if (donetrisglobal[k] != nowtri && vector_search(donetrisglobal[k], root) == 1 && vector_search(donetrisglobal[k], lever) == 1)
						{
							if (root != donetrisglobal[k][0] && lever != donetrisglobal[k][0])oldnode = donetrisglobal[k][0];
							else if (root != donetrisglobal[k][1] && lever != donetrisglobal[k][1])oldnode = donetrisglobal[k][1];
							else oldnode = donetrisglobal[k][2];
						}
					}
					Real dis1 = dist(_partpos.pos[root % _partpos.n_part], _partpos.pos[image% _partpos.n_part]);
					Real dis2 = dist(_partpos.pos[lever % _partpos.n_part], _partpos.pos[image%_partpos.n_part]);
					Real dism = dist(_partpos.pos[root % _partpos.n_part], _partpos.pos[lever % _partpos.n_part]);
					Vector3d rootpos, leverpos, oldpos;

					oldpos = { pos2D[ID2D[oldnode]][0],pos2D[ID2D[oldnode]][1],0.0 };
					for (int k = 0; k < midplots.size(); k++) {
						if (root == midplots[k].first) {
							rootpos = midplots[k].second;
						}
						if (lever == midplots[k].first) {
							leverpos = midplots[k].second;
						}
						if (oldnode == midplots[k].first) {
							oldpos = midplots[k].second;
						}
					}
					dism = dist(rootpos, leverpos);
					Real a = (dis1 * dis1 - dis2 * dis2 + dism * dism) / (2 * dism);
					Real h = sqrt(dis1 * dis1 - a * a);
					Vector3d PP0 = { rootpos[0] + a * (leverpos[0] - rootpos[0]) / dism,rootpos[1] + a * (leverpos[1] - rootpos[1]) / dism,0 };
					Vector3d PP1 = { PP0[0] + h * (leverpos[1] - rootpos[1]) / dism,PP0[1] - h * (leverpos[0] - rootpos[0]) / dism ,0 };
					Vector3d PP2 = { PP0[0] - h * (leverpos[1] - rootpos[1]) / dism,PP0[1] + h * (leverpos[0] - rootpos[0]) / dism ,0 };
					if (isnan(h)) {
						PP1 = PP2 = (rootpos + leverpos) / 2.0;
					}
					Vector3d aa1, bb1, cc1, First1, Second1;
					aa1 = { leverpos[0] - rootpos[0],leverpos[1] - rootpos[1], 0 };
					bb1 = { oldpos[0] - rootpos[0],oldpos[1] - rootpos[1], 0 };
					cc1 = { PP1[0] - rootpos[0],PP1[1] - rootpos[1], 0 };
					Second1 = aa1.cross(bb1) / aa1.cross(bb1).norm();
					First1 = aa1.cross(cc1) / aa1.cross(cc1).norm();
					Vector3d PP;
					if (First1[2] == Second1[2]) {
						midplots.push_back(make_pair(image, PP2));
						PP = PP2;
					}
					else {
						midplots.push_back(make_pair(image, PP1));
						PP = PP1;
					}
					if (image >= _partpos.n_part) {
						ID2D[image] = extid;
						pos2D.push_back({ 0,0 });
						pos2D[extid][0] = PP(0);
						pos2D[extid][1] = PP(1);
						extid++;
						f3f << extid - 1 << " " << image << endl;
					}
					else {
						pos2D[image][0] = PP(0);
						pos2D[image][1] = PP(1);
					}
					ord[image] = counter;
					ordrev[counter] = image;
					posord.push_back({PP(0),PP(1)});
					counter++;
					f2f << image << " ";
					mapped[image%_partpos.n_part] = 1;
					smalldone.push_back(image% _partpos.n_part);
					smalldone.push_back(lever% _partpos.n_part);
					smalldone.push_back(root% _partpos.n_part);
					alldone.push_back(image% _partpos.n_part);
					alldone.push_back(lever% _partpos.n_part);
					alldone.push_back(root% _partpos.n_part);
					previouspts[image] = image;
					donetrisglobal.push_back(nowtri);
					donetrisall[odsmalltris[j]] = 1;
					tricounter++;
					CalcDistortion(to_string(minimcounter));
					minimcounter++;
					cout << minimcounter << " " << tricounter << " in i!=0 printing new point " << totaltris << endl;
					if (minimcounter == stopminim)return;
				}
			}
			f2f << endl;
		}	
		f2f.close();
		f3f.close();
	}
	void UnfoldTriangles::CalcDistortion(string str) {
		map<string, int>midpts;
		set<string>midptstrs;
		map<int, int>neword;
		vector<Vector2d> posordnew;
		vector< pair<vector<int>, Real>> difs;
		vector< pair<vector<int>, Real>> trcs;
		vector< pair<vector<int>, Real>> js;
		if (minim==true) {
			posordnew = posord;
			neword = ord;
			// This section builds the midpoints for quadratic triangles
			for (int DTG = 0; DTG < donetrisglobal.size(); DTG++) {
				int p1 = donetrisglobal[DTG][0];
				int p2 = donetrisglobal[DTG][1];
				int p3 = donetrisglobal[DTG][2];
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
						int pm1 = ord[pt1];
						int pm2 = ord[pt2];
						Vector2d midpm = (posord[pm1] + posord[pm2]) / 2.0;
						posordnew.push_back(midpm);
						midpts[ptnow] = posordnew.size() - 1;
					}
				}
			}
			// Calls body function to build the model for running minimization
			Body bd(posordnew, donetrisglobal, ord, midpts, so, 0, _partpos.n_part);
			int totalcountersize = counter + midpts.size();
			VectorXd xo;
			xo.resize(totalcountersize * 2);
			vector<double>Pos2;
			Pos2.clear();
			Pos2 = bd.getPos();
			for (int l = 0; l < 2 * totalcountersize; l++) {
				xo[l] = Pos2[l];
			}
			Tris = donetrisglobal;
			ordr = ord;
			MIDP = midpts;
			s = _partpos.input_particlefile;
			nprt = _partpos.n_part;
			LBFGSBParam<double> param;
			param.epsilon = param_epsilon;
			param.epsilon_rel = param_epsilon_rel;
			param.delta = param_delta;
			//param.m = param_m;
			//param.max_linesearch = param_max_linesearch;
			//param.ftol = param_ftol;
			//param.wolfe = param_wolfe;
			//param.max_submin = param_max_submin;
			//param.max_iterations = param_max_iterations;
			LBFGSBSolver<double> solver(param);
			double fx;
			iter = 1;
			VectorXd lb = VectorXd::Constant(2 * totalcountersize, -1000.0);
			VectorXd ub = VectorXd::Constant(2 * totalcountersize, 1000.0);
			lb[0] = 0.0;
			ub[0] = 0.0;
			lb[1] = 0.0;
			ub[1] = 0.0;
			int niter = solver.minimize(foo, xo, fx, lb, ub);
			EnergyF.push_back(fx);
			//cout << fx << endl;
			for (int l = 0; l < counter; l++) {
				pos2D[ID2D[ordrev[l]]][0] = xo[2 * l + 0];
				pos2D[ID2D[ordrev[l]]][1] = xo[2 * l + 1];
				posord[l][0] = xo[2 * l + 0];
				posord[l][1] = xo[2 * l + 1];
				posordnew[l][0] = xo[2 * l + 0];
				posordnew[l][1] = xo[2 * l + 1];
			}
			for (int l = counter; l < (counter + midpts.size()); l++) {
				posordnew[l][0] = xo[2 * l + 0];
				posordnew[l][1] = xo[2 * l + 1];
			}
			for (int l = 0; l < midplots.size(); l++) {
				midplots[l].second = { pos2D[ID2D[midplots[l].first]][0],pos2D[ID2D[midplots[l].first]][1],0.0 };
			}
		}
		Body bd2(posordnew, donetrisglobal, ord, midpts, so, 1, _partpos.n_part);
		difs = bd2.getdiff();
		trcs = bd2.gettrcs();
		js = bd2.getjs();
		int totaltris = 0;
		vector<int>notprint;
		int problemco = 0;
		for (int DT = 0; DT < donetrisglobal.size(); DT++) {
			if (vector_search(notprint, DT) == 1) {
				continue;
			}
			vector<int>maint = { 0,0,0 };
			for (int j = 0; j < 3; j++) {
				maint[j] = donetrisglobal[DT][j] % _partpos.n_part;
			}
			sort(maint.begin(), maint.end());
			bool problem = 0;
			for (int j = 0; j < donetrisglobal.size(); j++) {
				if (DT != j) {
					vector<int>maing = { 0,0,0 };
					for (int r = 0; r < 3; r++) {
						maing[r] = donetrisglobal[j][r] % _partpos.n_part;
					}
					if (maint[0] == maing[0] && maint[1] == maing[1] && maint[2] == maing[2]) {
						problem = 1;
						problemco++;
						notprint.push_back(j);
					}
				}
			}
			totaltris++;
		}
		// This function handles printing of unfolded triangles with provided specifications about distortion
		if (stoi(str) > 0) {
			int init2dsize = pos2D.size();
			int initpossize = posord.size();
			int scheme = 1;
			ofstream f;
			ofstream ftr;
			ofstream fj;
			f.open(_partpos._folder+"/UnfoldedTrisAreaDistortion" + str + _partpos.input_particlefile);
			f << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS ";

			ftr.open(_partpos._folder + "/UnfoldedTrisTraceDistortion" + str + _partpos.input_particlefile);
			ftr << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS ";

			fj.open(_partpos._folder + "/UnfoldedTrisJacobianDistortion" + str + _partpos.input_particlefile);
			fj << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS ";

			if (scheme == 1) { 
				f << pos2D.size() + midpts.size();
				ftr << pos2D.size() + midpts.size();
				fj << pos2D.size() + midpts.size();
			}
			else {
				f << pos2D.size() << endl;
				ftr << pos2D.size() << endl;
				fj << pos2D.size() << endl;
			}
			f << " float\n";
			ftr << " float\n";
			fj << " float\n";
			for (int i = 0; i < pos2D.size()  ; i++) {
				f << pos2D[i][0] << " " << pos2D[i][1] << " " << 0 << endl;
				ftr << pos2D[i][0] << " " << pos2D[i][1] << " " << 0 << endl;
				fj << pos2D[i][0] << " " << pos2D[i][1] << " " << 0 << endl;
			}
			if (scheme == 1) {
				for (int l = counter; l < (counter + midpts.size()); l++) {
					f << posordnew[l][0] << " " << posordnew[l][1] << " " << 0 << endl;
					ftr << posordnew[l][0] << " " << posordnew[l][1] << " " << 0 << endl;
					fj << posordnew[l][0] << " " << posordnew[l][1] << " " << 0 << endl;
				}
			}
			f << endl;
			ftr << endl;
			fj << endl;
			f << "CELLS ";
			ftr << "CELLS ";
			fj << "CELLS ";
			if (scheme == 1) {
				f << totaltris << " " << 7 * totaltris << endl;
				ftr << totaltris << " " << 7 * totaltris << endl;
				fj << totaltris << " " << 7 * totaltris << endl;
			}
			else {
				f << totaltris << " " << 4 * totaltris << endl;
				ftr << totaltris << " " << 4 * totaltris << endl;
				fj << totaltris << " " << 4 * totaltris << endl;
			}
			if (scheme == 1) {
				for (int DTG = 0; DTG < donetrisglobal.size(); DTG++) {
					vector<int>pms;
					int p1 = donetrisglobal[DTG][0];
					int p2 = donetrisglobal[DTG][1];
					int p3 = donetrisglobal[DTG][2];
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
						pms.push_back(midpts[ptnow]);
					}
					if (vector_search(notprint, DTG) == 1)continue;
					f << 6 << " " << ID2D[donetrisglobal[DTG][0]] << " " << ID2D[donetrisglobal[DTG][1]] << " " << ID2D[donetrisglobal[DTG][2]] << " " << pms[0] + pos2D.size() - posord.size() << " " << pms[1] + pos2D.size() - posord.size() << " " << pms[2] + pos2D.size() - posord.size() << endl;
					ftr << 6 << " " << ID2D[donetrisglobal[DTG][0]] << " " << ID2D[donetrisglobal[DTG][1]] << " " << ID2D[donetrisglobal[DTG][2]] << " " << pms[0] + pos2D.size() - posord.size() << " " << pms[1] + pos2D.size() - posord.size() << " " << pms[2] + pos2D.size() - posord.size() << endl;
					fj << 6 << " " << ID2D[donetrisglobal[DTG][0]] << " " << ID2D[donetrisglobal[DTG][1]] << " " << ID2D[donetrisglobal[DTG][2]] << " " << pms[0] + pos2D.size() - posord.size() << " " << pms[1] + pos2D.size() - posord.size() << " " << pms[2] + pos2D.size() - posord.size() << endl;
				}
			}
			else {
				for (int DTG = 0; DTG < donetrisglobal.size(); DTG++) {
					if (vector_search(notprint, DTG) == 1)continue;
					f << 3 << " " << ID2D[donetrisglobal[DTG][0]] << " " << ID2D[donetrisglobal[DTG][1]] << " " << ID2D[donetrisglobal[DTG][2]] << endl;
					ftr << 3 << " " << ID2D[donetrisglobal[DTG][0]] << " " << ID2D[donetrisglobal[DTG][1]] << " " << ID2D[donetrisglobal[DTG][2]] << endl;
					fj << 3 << " " << ID2D[donetrisglobal[DTG][0]] << " " << ID2D[donetrisglobal[DTG][1]] << " " << ID2D[donetrisglobal[DTG][2]] << endl;
				}
			}
			f << endl;
			f << "CELL_TYPES ";
			f << totaltris << endl;

			ftr << endl;
			ftr << "CELL_TYPES ";
			ftr << totaltris << endl;

			fj << endl;
			fj << "CELL_TYPES ";
			fj << totaltris << endl;
			for (int i = 0; i < totaltris; i++) {
				if (scheme == 1) {
					f << 22 << endl;
					ftr << 22 << endl;
					fj << 22 << endl;
				}
				else {
					f << 5 << endl;
					ftr << 5 << endl;
					fj << 5 << endl;
				}
			}
			f << "\nCELL_DATA " << totaltris << "\n";
			f << "SCALARS AreaDistortion float\n";
			f << "LOOKUP_TABLE default\n";

			ftr << "\nCELL_DATA " << totaltris << "\n";
			ftr << "SCALARS AreaDistortion float\n";
			ftr << "LOOKUP_TABLE default\n";

			fj << "\nCELL_DATA " << totaltris << "\n";
			fj << "SCALARS AreaDistortion float\n";
			fj << "LOOKUP_TABLE default\n";
			Real totaldistortion = 0.0;
			for (int i = 0; i < donetrisglobal.size(); i++) {
				if (vector_search(notprint, i) == 1)continue;
				vector<int>maint = { 0,0,0 };
				for (int j = 0; j < 3; j++) {
					maint[j] = donetrisglobal[i][j] % _partpos.n_part;
				}
				sort(maint.begin(), maint.end());
				bool found = false;
				for (int j = 0; j < difs.size(); j++) {
					vector<int>ts = difs[j].first;
					sort(ts.begin(), ts.end());
					if (maint == ts) {
						f << difs[j].second << endl;
						ftr << trcs[j].second << endl;
						fj << js[j].second << endl;
						totaldistortion += difs[j].second;
						found = true;
						break;
					}
				}
			}
			f.close();
			ftr.close();
			fj.close();
		}		
	}
}