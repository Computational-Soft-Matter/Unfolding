#include "UnfoldTriangles.h"
#include <LBFGSB.h>
namespace unfolding {
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
		Body bdcomp(Pos, Tris, ordr, s, 0, 0, 0,nprt);
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
		double enrg = bdcomp.getenrg();
		return enrg;
	}
	double foo(const VectorXd& x, VectorXd& grad)
	{
		//cout << "iter" << iter << endl;
		iter++;
		double f;
		f = computeEnrg(x);
		//cout << "Energy: " << f << endl;
		EnergyF.push_back(f);
		grad = computeGrd(x);
		//cout << grad << endl;
		//cout << "grd" << endl;
		return f;
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
		int nn, int un)
	{
		// Base Case
		if (un == -1) {
			pathsn.push_back(pathn);
			return;
		}

		// Loop for all the parents
		// of the given vertex
		for (int par : parentn[un]) {

			// Insert the current
			// vertex in path
			pathn.push_back(un);

			// Recursive call for its parent
			find_paths(pathsn, pathn, parentn,
				nn, par);

			// Remove the current vertex
			pathn.pop_back();
		}
	}
	TriPtsDS UnfoldTriangles::TriPts(vector<int>T) {
		//print_vector(T);
		vector<int>final_vert;
		vector<vector<int>>trails;
		vector<Vector3d>plots;
		queue<int>path;
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
		for (int j = 0; j < 3; j++) {
			visitedpath.clear();
			visitedpath.resize(_partpos.n_part, 0);
			int source = tr_edge[j].first;
			int target = tr_edge[j].second;
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
			vector<int>graphdis(_partpos.n_part, 999);
			graphdis[source] = 0;
			for (int k = 0; k < allpaths.size(); k++) {
				if ((allpaths[k].first.first == source && allpaths[k].first.second == target) || (allpaths[k].first.first == source && allpaths[k].first.second == target)) {
					possiblepath = allpaths[k].second;
					//print_vector(possiblepath);
					break;
				}
			}
			vector<int>trail;
			if (possiblepath.size() > 0) {
				for (int k = 0; k < possiblepath.size(); k++) {
					if (visitedpath[possiblepath[k]] == false) {
						final_vert.push_back(possiblepath[k]);
						path.push(possiblepath[k]);
						visitedpath[possiblepath[k]] = true;
					}
				}
				trail = possiblepath;
				pathsnow.push_back(possiblepath);
			}
			else {
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
						if (!visited[adj]) {
							visited[adj] = true;
							bfs.push(adj);
							pathfinder[adj] = source;
						}
						if (graphdis[adj] > graphdis[source] + 1) {
							graphdis[adj] = graphdis[source] + 1;
							bfs.push(adj);
							parents[adj].clear();
							parents[adj].push_back(source);
						}
						else if (graphdis[adj] == graphdis[source] + 1) {
							parents[adj].push_back(source);
						}
						if (adj == target)stop = true;
					}
				}
				vector<int>targetparent = { target };
				vector<vector<int>>mulpaths;
				vector<int>singlepath;
				find_paths(mulpaths, singlepath, parents, _partpos.n_part, target);
				int biggestdist = -1;
				int pathdorkar;
				for (int k = 0; k < mulpaths.size(); k++) {
					Real totdis = 0;
					for (int kk = 0; kk < mulpaths[k].size(); kk++) {
						totdis += _partpos.MinDisGrp[othernode][mulpaths[k][kk]];
					}
					if (totdis > biggestdist) {
						biggestdist = totdis;
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
						//cout << target << " skyline " << endl;
						if(vector_search(final_vert,target)==0)
						final_vert.push_back(target);
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
		
		for (int rr = 0; rr < trails.size(); rr++) {
			int tester;
			print_vector(trails[rr]);
			for (int r = 0; r < 3; r++) {
				if (T[r] != tr_edge[rr].first && T[r] != tr_edge[rr].second) {
					tester = T[r];
				}
			}
			//cout << tester<<"tester" << endl;
			for (int r = 0; r < trails[rr].size(); r++) {
				int trailnow = trails[rr][r];
				cout << trailnow <<" "<< vector_search(_unfpt._tri._defectpos.defect_pos, trailnow) << endl;
				if (tester == 0) {
					//if (vector_search(_unfpt._tri._defectpos.defect_pos, trailnow) == 1)continue;
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
				if (totalfns >= 3 && tester!=0)continue;
				//cout << trailnow <<"ankhiyan "<< endl;
				vector<Real>distdone;
				for (int l = 0; l < 2; l++) {
					//print_vector(donesets[l]);
					//cout << l << endl;
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
				if (tester == 0) {
					cout << trailnow << endl;
						print_vector(donesets[0]);
						print_vector(donesets[1]);
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
					//print_vector(distdone);
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
			//cout << final_vert[j] << endl;
			int vn_ngh = final_vert[j];
			if (vn_ngh > _partpos.n_part)continue;
			Real x = _partpos.pos[vn_ngh][0];
			Real y = _partpos.pos[vn_ngh][1];
			Real z = _partpos.pos[vn_ngh][2];
			Real t = (n(0) * d - n(0) * x + n(1) * e - n(1) * y + n(2) * f - n(2) * z) / (n(0) * n(0) + n(1) * n(1) + n(2) * n(2));
			Vector3d proj = { x + t * n(0),y + t * n(1),z + t * n(2) };
			plots.push_back(_partpos.pos[vn_ngh]);
		}
		//path.push(n1);
		//path.push(n2);
		//path.push(n3);
		
		vector<int>testvec;
		vector<Vector3d>plotsbackup;
		visitedpath.clear();
		visitedpath.resize(_partpos.n_part, false);
		visitedpath[n1] = 1;
		visitedpath[n2] = 1;
		visitedpath[n3] = 1;
		while (!path.empty()) {
			int verticenow = path.front();
			//cout << verticenow << endl;
			path.pop();
			if (verticenow > _partpos.n_part)continue;
			int vnsize = _partpos.pos_neighbor[verticenow].size();
			for (int k = 0; k < vnsize; k++) {
				int vn_ngh = _partpos.pos_neighbor[verticenow][k];
				//cout << vn_ngh << endl;
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
						if(vector_search(final_vert,vn_ngh)==0)final_vert.push_back(vn_ngh);
						plots.push_back(_partpos.pos[vn_ngh]);
					}
					/*else {
						if (verticenow != n1 && verticenow != n2 && verticenow != n3)
						{
							//visitedpath[vn_ngh] = true;
							//testvec.push_back(vn_ngh);
							//plotsbackup.push_back(_partpos.pos[vn_ngh]);
						}
					}*/
				}
			}
		}
		//for (int k = 0; k < testvec.size(); k++) {
			//final_vert.push_back(testvec[k]);
			//plots.push_back(plotsbackup[k]);
		//}
		TriPtsDS tempDS;
		tempDS.Vrt = final_vert;
		tempDS.Plts = plots;
		return tempDS;
	}
	void UnfoldTriangles::Unfold() {
		vector<pair<int, int>>doneedge;
		for (int i = 0; i < 5/*_unfpt.t_size*/; i++) {
			vector<int>previouscompare;
			vector<int>unfnow = _unfpt.plottedfinalpts[i];
			vector<int>unfnoworg = { unfnow[0] % _partpos.n_part, unfnow[1] % _partpos.n_part, unfnow[2] % _partpos.n_part };
			TriPtsDS DSnow = TriPts(unfnoworg);
			vector<int>final_vert = DSnow.Vrt;
			vector<Vector3d>plots = DSnow.Plts;
			int o1 = unfnow[0];
			int o2 = unfnow[1];
			int o3 = unfnow[2];
			cout << i << endl;
			cout << o1 << " " << o2 << " " << o3 << endl;
			cout << "final_vert: " << endl;
			print_vector(final_vert);
			nountripts.push_back(vector<int>());
			// Finding all the small triangles inside the triangular region bounded by o1,o2,o3
			vector<vector<int>>smalltris;

			for (int j = 0; j < final_vert.size(); j++) {
				int indn = final_vert[j];
				//cout << indn << "now" << endl;
				int sizenow = _partpos.pos_neighbor[indn].size();
				//print_vector(_partpos.pos_neighbor[indn]);
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
						sort(P.begin(), P.end());
						bool pushnow = 0;
						for (int l = 0; l < smalltris.size(); l++) {
							if (smalltris[l] == P)pushnow = 1;
						}
						if (pushnow == 0) {
							smalltris.push_back(P);
							//print_vector(P);
						}
					}
				}
			}
			queue<vector<int>>donetris;
			vector<int>smalldone;
			int firstnode, secondnode;
			int donesmalltris = 0;
			int startvrt = -5;
			vector<int>previousorgpts;
			map<int, int>previouspts;
			if (i == 0)startvrt = o1;
			vector<pair<int, Vector3d>>midplots;
			if (i != 0) {
				for (int j = 0; j < final_vert.size(); j++) {
					previouspts[final_vert[j]] = -5;
					for (int r = 0; r < untripts[_unfpt.finalparenttriangle[i]].size(); r++) {
						if (untripts[_unfpt.finalparenttriangle[i]][r] % _partpos.n_part == final_vert[j]) {
							previousorgpts.push_back(final_vert[j]);
							previouspts[final_vert[j]] = untripts[_unfpt.finalparenttriangle[i]][r];
							if (startvrt == -5) {
								for (int rk = 0; rk < smalltris.size(); rk++) {
									if (vector_search(smalltris[rk], final_vert[j] == 1)) {
										startvrt = final_vert[j];
										break;
									}
								}
							}
							//if (startvrt == -5) startvrt = final_vert[j];
							Vector3d placeholder = { pos2D[ID2D[untripts[_unfpt.finalparenttriangle[i]][r]]][0],pos2D[ID2D[untripts[_unfpt.finalparenttriangle[i]][r]]][1],0.0 };
							midplots.push_back(make_pair(untripts[_unfpt.finalparenttriangle[i]][r], placeholder));
							smalldone.push_back(untripts[_unfpt.finalparenttriangle[i]][r] % _partpos.n_part);
						}
					}
				}
			}
			vector<vector<int>>odsmalltris;
			vector<bool>mappedtrism;
			mappedtrism.resize(smalltris.size(), false);
			queue<vector<int>>unftrsm;
			for (int j = 0; j < smalltris.size(); j++) {
				if (i == 0 && vector_search(smalltris[j],startvrt)==1) {
					unftrsm.push(smalltris[j]);
					mappedtrism[j] = true;
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
			}
			if (i != 0) {
				untripts.push_back(vector<int>());
				for (int j = 0; j < previousorgpts.size(); j++) {
					untripts[i].push_back(previouspts[previousorgpts[j]]);
				}
			}
			if (i == 0) {
				Real dis = dist(_partpos.pos[firstnode], _partpos.pos[secondnode]);
				Vector3d first = { 0,0,0 };
				Vector3d second = { dis,0,0 };
				pos2D[firstnode][0] = first(0);
				pos2D[firstnode][1] = first(1);
				pos2D[secondnode][0] = second(0);
				pos2D[secondnode][1] = second(1);
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
				for (int j = 0; j < odsmalltris.size(); j++) {
					donetrisglobal.push_back(odsmalltris[j]);
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
						cout << "didnotplot " << minimcounter << endl;
						if (minim == true) {
							Body bd(posord, donetrisglobal, ord, so, 0, 0, 0, _partpos.n_part);
							VectorXd xo;
							xo.resize(counter * 2);
							vector<double>Pos2;
							Pos2.clear();
							Pos2 = bd.getPos();
							for (int l = 0; l < 2 * counter; l++) {
								xo[l] = Pos2[l];
							}
							Tris = donetrisglobal;
							ordr = ord;
							s = _partpos.input3Dcon;
							nprt = _partpos.n_part;
							LBFGSBParam<double> param;
							param.epsilon = 1e-10;
							param.delta = 1e-15;
							LBFGSBSolver<double> solver(param);
							double fx;
							iter = 1;
							VectorXd lb = VectorXd::Constant(2 * counter, -1000.0);
							VectorXd ub = VectorXd::Constant(2 * counter, 1000.0);
							lb[0] = 0.0;
							ub[0] = 0.0;
							lb[1] = 0.0;
							ub[1] = 0.0;
							lb[3] = 0.0;
							ub[3] = 0.0;
							int niter = solver.minimize(foo, xo, fx, lb, ub);
							EnergyF.push_back(fx);
							for (int l = 0; l < counter; l++) {
								pos2D[ID2D[ordrev[l]]][0] = xo[2 * l + 0];
								pos2D[ID2D[ordrev[l]]][1] = xo[2 * l + 1];
								posord[l][0] = xo[2 * l + 0];
								posord[l][1] = xo[2 * l + 1];
							}
							for (int l = 0; l < midplots.size(); l++) {
								midplots[l].second = { pos2D[ID2D[midplots[l].first]][0],pos2D[ID2D[midplots[l].first]][1],0.0 };
							}
						}
						CalcDistortion(to_string(minimcounter));
						minimcounter++;
						
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
					//cout << oldnode << endl;
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
					smalldone.push_back(image);
					pos2D[image][0] = PP(0);
					pos2D[image][1] = PP(1);
					mapped[image] = 1;
					if (minim == true) {
						cout << "didplot " << minimcounter << endl;
						Body bd(posord, donetrisglobal, ord, so, 0, 0, 0, _partpos.n_part);
						VectorXd xo;
						xo.resize(counter * 2);
						vector<double>Pos2;
						Pos2.clear();
						Pos2 = bd.getPos();
						for (int l = 0; l < 2 * counter; l++) {
							xo[l] = Pos2[l];
						}
						Tris = donetrisglobal;
						ordr = ord;
						s = _partpos.input3Dcon;
						nprt = _partpos.n_part;
						LBFGSBParam<double> param;
						param.epsilon = 1e-10;
						param.delta = 1e-15;
						LBFGSBSolver<double> solver(param);
						double fx;
						iter = 1;
						VectorXd lb = VectorXd::Constant(2 * counter, -1000.0);
						VectorXd ub = VectorXd::Constant(2 * counter, 1000.0);
						lb[0] = 0.0;
						ub[0] = 0.0;
						lb[1] = 0.0;
						ub[1] = 0.0;
						lb[3] = 0.0;
						ub[3] = 0.0;
						int niter = solver.minimize(foo, xo, fx, lb, ub);
						EnergyF.push_back(fx);
						for (int l = 0; l < counter; l++) {
							pos2D[ID2D[ordrev[l]]][0] = xo[2 * l + 0];
							pos2D[ID2D[ordrev[l]]][1] = xo[2 * l + 1];
							posord[l][0] = xo[2 * l + 0];
							posord[l][1] = xo[2 * l + 1];
						}
						for (int l = 0; l < midplots.size(); l++) {
							midplots[l].second = { pos2D[ID2D[midplots[l].first]][0],pos2D[ID2D[midplots[l].first]][1],0.0 };
						}
					}
					CalcDistortion(to_string(minimcounter));
					minimcounter++;
				}
				untripts.push_back(final_vert);
			}
			else {
				for (int j = 0; j < odsmalltris.size(); j++) {
					int root, lever, image;
					image = -5;
					vector<int>nowtri = odsmalltris[j];
					print_vector(nowtri);
					int rchange;
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
					if (image == -5) {
						for (int r = 0; r < 3; r++) {
							smalldone.push_back(nowtri[r] % _partpos.n_part);
						}
						donetrisglobal.push_back(nowtri);
						cout << "didnotplot " << minimcounter << endl;
						if (minim == true) {
							Body bd(posord, donetrisglobal, ord, so, 0, 0, 0, _partpos.n_part);
							VectorXd xo;
							xo.resize(counter * 2);
							vector<double>Pos2;
							Pos2.clear();
							Pos2 = bd.getPos();
							for (int l = 0; l < 2 * counter; l++) {
								xo[l] = Pos2[l];
							}
							Tris = donetrisglobal;
							ordr = ord;
							s = _partpos.input3Dcon;
							nprt = _partpos.n_part;
							LBFGSBParam<double> param;
							param.epsilon = 1e-10;
							param.delta = 1e-15;
							LBFGSBSolver<double> solver(param);
							double fx;
							iter = 1;
							VectorXd lb = VectorXd::Constant(2 * counter, -1000.0);
							VectorXd ub = VectorXd::Constant(2 * counter, 1000.0);
							lb[0] = 0.0;
							ub[0] = 0.0;
							lb[1] = 0.0;
							ub[1] = 0.0;
							lb[3] = 0.0;
							ub[3] = 0.0;
							int niter = solver.minimize(foo, xo, fx, lb, ub);
							EnergyF.push_back(fx);
							for (int l = 0; l < counter; l++) {
								pos2D[ID2D[ordrev[l]]][0] = xo[2 * l + 0];
								pos2D[ID2D[ordrev[l]]][1] = xo[2 * l + 1];
								posord[l][0] = xo[2 * l + 0];
								posord[l][1] = xo[2 * l + 1];
							}
							for (int l = 0; l < midplots.size(); l++) {
								midplots[l].second = { pos2D[ID2D[midplots[l].first]][0],pos2D[ID2D[midplots[l].first]][1],0.0 };
							}
						}
						CalcDistortion(to_string(minimcounter));
						minimcounter++;
						continue;
					}
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

					//cout << image << " image " <<root<< " "<<lever<< endl;
					int oldnode = -1;

					for (int k = 0; k < donetrisglobal.size(); k++) {
						if (donetrisglobal[k] != nowtri && vector_search(donetrisglobal[k], root) == 1 && vector_search(donetrisglobal[k], lever) == 1)
						{
							if (root != donetrisglobal[k][0] && lever != donetrisglobal[k][0])oldnode = donetrisglobal[k][0];
							else if (root != donetrisglobal[k][1] && lever != donetrisglobal[k][1])oldnode = donetrisglobal[k][1];
							else oldnode = donetrisglobal[k][2];
						}
					}
					//cout << oldnode << endl;
					Real dis1 = dist(_partpos.pos[root % _partpos.n_part], _partpos.pos[image% _partpos.n_part]);
					Real dis2 = dist(_partpos.pos[lever % _partpos.n_part], _partpos.pos[image%_partpos.n_part]);
					Real dism = dist(_partpos.pos[root % _partpos.n_part], _partpos.pos[lever % _partpos.n_part]);
					//cout<<dis1<<" "<<dis2<<" "<<dism<<endl;
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
						//cout << extid - 1 <<"extra "<< endl;
					}
					else {
						pos2D[image][0] = PP(0);
						pos2D[image][1] = PP(1);
					}
					//cout << PP << endl;
					ord[image] = counter;
					ordrev[counter] = image;
					posord.push_back({PP(0),PP(1)});
					counter++;
					mapped[image%_partpos.n_part] = 1;
					smalldone.push_back(image% _partpos.n_part);
					smalldone.push_back(lever% _partpos.n_part);
					smalldone.push_back(root% _partpos.n_part);
					previouspts[image] = image;
					donetrisglobal.push_back(nowtri);
					if (minim == true) {
						cout << "didplot " << minimcounter << endl;
						Body bd(posord, donetrisglobal, ord, so, 0, 0, 0, _partpos.n_part);
						VectorXd xo;
						xo.resize(counter * 2);
						vector<double>Pos2;
						Pos2.clear();
						Pos2 = bd.getPos();
						for (int l = 0; l < 2 * counter; l++) {
							xo[l] = Pos2[l];
						}
						Tris = donetrisglobal;
						ordr = ord;
						s = _partpos.input3Dcon;
						nprt = _partpos.n_part;
						LBFGSBParam<double> param;
						param.epsilon = 1e-10;
						param.delta = 1e-15;
						LBFGSBSolver<double> solver(param);
						double fx;
						iter = 1;
						VectorXd lb = VectorXd::Constant(2 * counter, -1000.0);
						VectorXd ub = VectorXd::Constant(2 * counter, 1000.0);
						lb[0] = 0.0;
						ub[0] = 0.0;
						lb[1] = 0.0;
						ub[1] = 0.0;
						lb[3] = 0.0;
						ub[3] = 0.0;
						int niter = solver.minimize(foo, xo, fx, lb, ub);
						EnergyF.push_back(fx);
						for (int l = 0; l < counter; l++) {
							pos2D[ID2D[ordrev[l]]][0] = xo[2 * l + 0];
							pos2D[ID2D[ordrev[l]]][1] = xo[2 * l + 1];
							posord[l][0]= xo[2 * l + 0];
							posord[l][1] = xo[2 * l + 1];
						}
						for (int l = 0; l < midplots.size(); l++) {
							midplots[l].second = { pos2D[ID2D[midplots[l].first]][0],pos2D[ID2D[midplots[l].first]][1],0.0 };
						}
					}
					CalcDistortion(to_string(minimcounter));
					minimcounter++;
				}
			}
			for (int j = 0; j < final_vert.size(); j++) {
				if (vector_search(smalldone, final_vert[j]) == 0) {
					//cout <<"missing"<< final_vert[j] << endl;
					int image = final_vert[j];
					int indexnow = image;
					if (mapped[image % _partpos.n_part] == 1) {
						int prevsize = rep_pts[indexnow].size();
						indexnow = rep_pts[indexnow][prevsize - 1] + _partpos.n_part;
						rep_pts[indexnow % _partpos.n_part].push_back(indexnow);
						image = indexnow;
					}
					//cout << image << endl;
					previouspts[image % _partpos.n_part] = image;
					untripts[i].push_back(image);
					nountripts[i].push_back(image% _partpos.n_part);
					int root = midplots[0].first;
					Vector3d rootpos = midplots[0].second;
					int lever = midplots[1].first;
					Vector3d leverpos = midplots[1].second;
					int lever2 = midplots[midplots.size()-1].first;
					Vector3d lever2pos = midplots[midplots.size() - 1].second;
					//cout << root << " " << lever << " " << lever2 << " " << endl;
					//cout << rootpos << " r " << leverpos << " l " << lever2pos << endl;
					Real dis1 = dist(_partpos.pos[root % _partpos.n_part], _partpos.pos[image % _partpos.n_part]);
					Real dis2 = dist(_partpos.pos[lever % _partpos.n_part], _partpos.pos[image % _partpos.n_part]);
					Real dism = dist(_partpos.pos[root % _partpos.n_part], _partpos.pos[lever % _partpos.n_part]);
					dism = dist(rootpos, leverpos);
					Real a = (dis1 * dis1 - dis2 * dis2 + dism * dism) / (2 * dism);
					Real h = sqrt(dis1 * dis1 - a * a);
					Vector3d PP0 = { rootpos[0] + a * (leverpos[0] - rootpos[0]) / dism,rootpos[1] + a * (leverpos[1] - rootpos[1]) / dism,0 };
					Vector3d PP1 = { PP0[0] + h * (leverpos[1] - rootpos[1]) / dism,PP0[1] - h * (leverpos[0] - rootpos[0]) / dism ,0 };
					Vector3d PP2 = { PP0[0] - h * (leverpos[1] - rootpos[1]) / dism,PP0[1] + h * (leverpos[0] - rootpos[0]) / dism ,0 };
					if (isnan(h)) {
						PP1 = PP2 = (rootpos + leverpos) / 2.0;
					}
					Real dislever2 = ceil(dist(_partpos.pos[lever2 % _partpos.n_part], _partpos.pos[image % _partpos.n_part]) * 1000.0) / 1000.0;
					Real dispp1 = dist(PP1, lever2pos);
					Real dispp2 = dist(PP2, lever2pos);
					Vector3d PP;
					if (abs(dispp1 - dislever2) < abs(dispp2 - dislever2)) {
						PP = PP1;
					}
					else PP = PP2;
					if (image >= _partpos.n_part) {
						ID2D[image] = extid;
						pos2D.push_back({ 0,0 });
						pos2D[extid][0] = PP(0);
						pos2D[extid][1] = PP(1);
						extid++;
						//cout << extid - 1 << "extra " << endl;
					}
					else {
						pos2D[image][0] = PP(0);
						pos2D[image][1] = PP(1);
					}
					ord[image] = counter;
					ordrev[counter] = image;
					posord.push_back({PP(0),PP(1)});
					counter++;
					mapped[image % _partpos.n_part] = 1;
					smalldone.push_back(image% _partpos.n_part);
					previouspts[image] = image;
				}

			}
		}	
	}
	void UnfoldTriangles::CalcDistortion(string str) {
		if (0) {
			Body bd(posord, donetrisglobal, ord, so, 0, 0, 0, _partpos.n_part);
			VectorXd xo;
			xo.resize(counter * 2);
			vector<double>Pos2;
			Pos2.clear();
			Pos2 = bd.getPos();
			for (int l = 0; l < 2 * counter; l++) {
				xo[l] = Pos2[l];
			}
			Tris = donetrisglobal;
			ordr = ord;
			s = _partpos.input3Dcon;
			nprt = _partpos.n_part;
			LBFGSBParam<double> param;
			param.epsilon = 1e-10;
			param.delta = 1e-15;
			LBFGSBSolver<double> solver(param);
			double fx;
			iter = 1;
			VectorXd lb = VectorXd::Constant(2 * counter, -1000.0);
			VectorXd ub = VectorXd::Constant(2 * counter, 1000.0);
			lb[0] = 0.0;
			ub[0] = 0.0;
			lb[1] = 0.0;
			ub[1] = 0.0;
			lb[3] = 0.0;
			ub[3] = 0.0;
			int niter = solver.minimize(foo, xo, fx, lb, ub);
			EnergyF.push_back(fx);
			for (int l = 0; l < counter; l++) {
				pos2D[ID2D[ordrev[l]]][0] = xo[2 * l + 0];
				pos2D[ID2D[ordrev[l]]][1] = xo[2 * l + 1];
				posord[l][0] = xo[2 * l + 0];
				posord[l][1] = xo[2 * l + 1];
			}
		}
		Body bd2(posord, donetrisglobal, ord, so, 0, 0, 0, _partpos.n_part);
		vector< pair<vector<int>, Real>> difs;
		difs = bd2.checkDiffTri();
		pair<vector< pair<vector<int>, Real>>, vector< pair<vector<int>, Real>>> resis;
		resis = bd2.plotNeoHookeanGrd2D();
		Real total = bd2.checkTotalSurfaceTri();
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
			bool problem = 0;
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
		f.open("UnfoldedTrisDistortion"+str + _partpos.input_particlefile);
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
			if (vector_search(notprint, i) == 1)continue;
			f << 3 << " " << ID2D[donetrisglobal[i][0]] << " " << ID2D[donetrisglobal[i][1]] << " " << ID2D[donetrisglobal[i][2]] << endl;
		}
		f << endl;
		f << "CELL_TYPES ";
		f << totaltris << endl;
		for (int i = 0; i < totaltris; i++) {
			f << 5 << endl;
		}
		f << "\nCELL_DATA " << totaltris << "\n";
		f << "SCALARS neighbour float\n";
		f << "LOOKUP_TABLE default\n";
		Real totaldistortion = 0.0;
		for (int i = 0; i < donetrisglobal.size(); i++) {
			if (vector_search(notprint, i) == 1)continue;
			vector<int>maint = { 0,0,0 };
			for (int j = 0; j < 3; j++) {
				maint[j] = donetrisglobal[i][j] % _partpos.n_part;
			}
			sort(maint.begin(), maint.end());
			//print_vector(maint);
			bool found = false;
			for (int j = 0; j < difs.size(); j++) {
				vector<int>ts = difs[j].first;
				sort(ts.begin(), ts.end());
				if (maint == ts) {
					f << difs[j].second << endl;
					totaldistortion += difs[j].second;
					//print_vector(ts);
					found = true;
					break;
				}
			}
			//if (found == 0)print_vector(maint);
		}
		cout << "Total distortion is: " << totaldistortion/total << endl;
		f.close();


		ofstream f1;
		f1.open("UnfoldedTrisI1bar" + str + _partpos.input_particlefile);
		f1 << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS ";
		f1 << pos2D.size();
		f1 << " float\n";
		for (int i = 0; i < pos2D.size(); i++) {
			f1 << pos2D[i][0] << " " << pos2D[i][1] << " " << 0 << endl;
		}
		f1 << endl;
		f1 << "CELLS ";
		f1 << totaltris << " " << 4 * totaltris << endl;
		for (int i = 0; i < donetrisglobal.size(); i++) {
			if (vector_search(notprint, i) == 1)continue;
			f1 << 3 << " " << ID2D[donetrisglobal[i][0]] << " " << ID2D[donetrisglobal[i][1]] << " " << ID2D[donetrisglobal[i][2]] << endl;
		}
		f1 << endl;
		f1 << "CELL_TYPES ";
		f1 << totaltris << endl;
		for (int i = 0; i < totaltris; i++) {
			f1 << 5 << endl;
		}
		f1 << "\nCELL_DATA " << totaltris << "\n";
		f1 << "SCALARS neighbour float\n";
		f1 << "LOOKUP_TABLE default\n";
		for (int i = 0; i < donetrisglobal.size(); i++) {
			if (vector_search(notprint, i) == 1)continue;
			vector<int>maint = { 0,0,0 };
			for (int j = 0; j < 3; j++) {
				maint[j] = donetrisglobal[i][j] % _partpos.n_part;
			}
			sort(maint.begin(), maint.end());
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
		f3.open("UnfoldedTrisI3trace" + str + _partpos.input_particlefile);
		f3 << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS ";
		f3 << pos2D.size();
		f3 << " float\n";
		for (int i = 0; i < pos2D.size(); i++) {
			f3 << pos2D[i][0] << " " << pos2D[i][1] << " " << 0 << endl;
		}
		f3 << endl;
		f3 << "CELLS ";
		f3 << totaltris << " " << 4 * totaltris << endl;
		for (int i = 0; i < donetrisglobal.size(); i++) {
			if (vector_search(notprint, i) == 1)continue;
			f3 << 3 << " " << ID2D[donetrisglobal[i][0]] << " " << ID2D[donetrisglobal[i][1]] << " " << ID2D[donetrisglobal[i][2]] << endl;
		}
		f3 << endl;
		f3 << "CELL_TYPES ";
		f3 << totaltris << endl;
		for (int i = 0; i < totaltris; i++) {
			f3 << 5 << endl;
		}
		f3 << "\nCELL_DATA " << totaltris << "\n";
		f3 << "SCALARS neighbour float\n";
		f3 << "LOOKUP_TABLE default\n";
		for (int i = 0; i < donetrisglobal.size(); i++) {
			if (vector_search(notprint, i) == 1)continue;
			vector<int>maint = { 0,0,0 };
			for (int j = 0; j < 3; j++) {
				maint[j] = donetrisglobal[i][j] % _partpos.n_part;
			}
			sort(maint.begin(), maint.end());
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

		
	}
}