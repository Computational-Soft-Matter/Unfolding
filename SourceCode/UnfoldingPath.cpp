#include "UnfoldingPath.h"
namespace unfolding {
	// Comparator function to sort edges class
	bool CompareDataedges(edges a, edges b)
	{
		return a.cost < b.cost;
	}
	//UnionFind algorithm - helper function for MST algorithm
	class UnionFind {
	public:
		UnionFind(int sz) : root(sz), rank(sz) {
			for (int i = 0; i < sz; i++) {
				root[i] = i;
				rank[i] = 1;
			}
		}

		int find(int x) {
			if (x == root[x]) {
				return x;
			}
			return root[x] = find(root[x]);
		}

		void unionSet(int x, int y) {
			int rootX = find(x);
			int rootY = find(y);
			if (rootX != rootY) {
				if (rank[rootX] > rank[rootY]) {
					root[rootY] = rootX;
				}
				else if (rank[rootX] < rank[rootY]) {
					root[rootX] = rootY;
				}
				else {
					root[rootY] = rootX;
					rank[rootX] += 1;
				}
			}
		}

		bool connected(int x, int y) {
			return find(x) == find(y);
		}

	private:
		vector<int> root;
		vector<int> rank;
	};
	void UnfoldingPath::buildAdjacency() {
		triconn.resize(t_size, vector<int>());
		for (int j = 0; j < t_size; j++) {
			vector<int>tnow = _tri.triangles[j];
			pair<int, int>pair1 = make_pair(tnow[0], tnow[1]);
			pair<int, int>pair2 = make_pair(tnow[0], tnow[2]);
			pair<int, int>pair3 = make_pair(tnow[2], tnow[1]);
			for (int i = 0; i < t_size; i++) {
				if (vector_search(_tri.triangles[i], pair1.first) == 1 && vector_search(_tri.triangles[i], pair1.second) == 1 && i != j) {
					triconn[j].push_back(i);
				}
				if (vector_search(_tri.triangles[i], pair2.first) == 1 && vector_search(_tri.triangles[i], pair2.second) == 1 && i != j) {
					triconn[j].push_back(i);
				}
				if (vector_search(_tri.triangles[i], pair3.first) == 1 && vector_search(_tri.triangles[i], pair3.second) == 1 && i != j) {
					triconn[j].push_back(i);
				}
			}
		}
	}
	void UnfoldingPath::buildEdges() {
		map<pair<int, int>, bool>edgedone;
		for (int i = 0; i < t_size; i++) {
			for (int j = 0; j < triconn[i].size(); j++) {
				pair<int, int>p;
				if (i < triconn[i][j]) p = make_pair(i, triconn[i][j]);
				else p = make_pair(triconn[i][j], i);
				if (edgedone[p] == 0) {
					ed.push_back(edges(p, (float)rand() / RAND_MAX , 1));
					edfnl.push_back(edges(p, 0.0, 1));
					edgedone[p] = 1;
				}
			}
		}
		edsize = ed.size();
	}
	void UnfoldingPath::printAdjacency() {
		for (int i = 0; i < t_size; i++) {
			cout << " Triangle " << i << " is: " << endl;
			print_vector(_tri.triangles[i]);
			cout << "Triangle " << i << " is connected to triangles: ";
			print_vector(triconn[i]);
		}
	}
	void UnfoldingPath::sortEdges() {
		for (int si = 0; si < edsize; si++) {
			for (int sj = 0; sj < edsize - 1; sj++) {
				if (ed[sj].cost > ed[sj + 1].cost) {
					Real tempcost = ed[sj + 1].cost;
					pair<int, int>temppair = ed[sj + 1].edge;
					int tempvisited = ed[sj + 1].visited;
					ed[sj + 1] = ed[sj];
					ed[sj].cost = tempcost;
					ed[sj].visited = tempvisited;
					ed[sj].edge = temppair;
				}
			}
		}
	}
	void UnfoldingPath::sortEdgesfnl() {
		for (int si = 0; si < edsize; si++) {
			for (int sj = 0; sj < edsize - 1; sj++) {
				if (edfnl[sj].cost > edfnl[sj + 1].cost) {
					Real tempcost = edfnl[sj + 1].cost;
					pair<int, int>temppair = edfnl[sj + 1].edge;
					int tempvisited = edfnl[sj + 1].visited;
					edfnl[sj + 1] = edfnl[sj];
					edfnl[sj].cost = tempcost;
					edfnl[sj].visited = tempvisited;
					edfnl[sj].edge = temppair;
				}
			}
		}
	}
	void UnfoldingPath::runMST(int steps) {
		for (int mst_n = 0; mst_n < steps; mst_n++) {					// Runs MST algorithm for 'steps' amount of iteration
			vector<pair<int, int>>mstedges;
			if (mst_n<randlimit) {										// Builds an unfolding path using "ed" list. This is the exploration stage. Runs for randlimit number of iterations.
				sortEdges();
				UnionFind uf(t_size);
				for (int e = 0; e < ed.size(); e++) {
					int node1 = ed[e].edge.first;
					int node2 = ed[e].edge.second;
					if (uf.connected(node1, node2) != 1) {
						uf.unionSet(node1, node2);
						mstedges.push_back(make_pair(ed[e].edge.first, ed[e].edge.second));
					}
				}
			}
			else {														// Builds an unfolding path using "edfnl" list. This uses the list of edges with weights calculated from exploration and searches for slightly better conditions. (exploitation) 
				sortEdgesfnl();
				UnionFind uf(t_size);
				for (int e = 0; e < edfnl.size(); e++) {
					int node1 = edfnl[e].edge.first;
					int node2 = edfnl[e].edge.second;
					if (uf.connected(node1, node2) != 1) {
						uf.unionSet(node1, node2);
						mstedges.push_back(make_pair(edfnl[e].edge.first, edfnl[e].edge.second));
					}
				}
			}
			
			vector<int>orderingfromedges = edgestoTris(mstedges);		// Indices of large triangles in the order of unfolding path
			// Builds the same list but with the triangles expressed as list of vertices
			for (int i = 0; i < t_size; i++) {
				orderedtriangles.push_back(vector<int>());
				orderedtriangles[i] = _tri.triangles[orderingfromedges[i]];
			}
			int msz = mstedges.size();
			
			// This section finds the number of  nodes in the graph that has more than 2 neighbors
			vector<vector<int>>n_nodes;
			n_nodes.resize(t_size, vector<int>());
			for (int i = 0; i < msz; i++) {
				int edge1 = mstedges[i].first;
				int edge2 = mstedges[i].second;
				for (int j = 0; j < edfnl.size(); j++) {
					int ed0 = edfnl[j].edge.first;
					int ed1 = edfnl[j].edge.second;
					if ((ed0 == edge1 && ed1 == edge2) || (ed1 == edge2 && ed0 == edge1)) {
						n_nodes[ed0].push_back(ed1);
						n_nodes[ed1].push_back(ed0);
					}
				}
			}
			int badtriangles = 0;
			for (int i = 0; i < n_nodes.size(); i++) {
				if (n_nodes[i].size() == 3) {
					badtriangles += 1;
				}
			}
			
			Real refoldingcost = Real(badtriangles) / Real(n_nodes.size());				// This is the weight that expresses if the unfolding path is strip
			Real costfinal = buildfinalTriangles(orderingfromedges,refoldingcost);		// Calls the build triangles function that prints the triangle vertices in the given unfolding order
																						// Subsequently calculates the weight (costfinal) depending on the overall unfolded vertices positions. 
			
			
			for (int i = 0; i < msz; i++) {								// At every iteration the 'edfnl' list is update with new weights
				int edge1 = mstedges[i].first;
				int edge2 = mstedges[i].second;
				for (int j = 0; j < edfnl.size(); j++) {
					int ed0 = edfnl[j].edge.first;
					int ed1 = edfnl[j].edge.second;
					if ((ed0 == edge1 && ed1 == edge2) || (ed1 == edge2 && ed0 == edge1)) {
						edfnl[j].cost = (edfnl[j].visited * edfnl[j].cost + costfinal) / (edfnl[j].visited + 1);
						edfnl[j].visited++;
					}
				}
			}

			if (mst_n < randlimit) {									// If exploring, this function populates 'ed' list with random weights again
				for (int j = 0; j < ed.size(); j++) {
					ed[j].cost = (float)rand() / RAND_MAX;
					ed[j].visited++;
				}
			}
			else {														
				if (mst_n % 50 == 0) {									// At every 50 steps in the exploitation phase, declusters the unfolding path in search for better paths
					vector<edges>decluster;
					Real maxcost = -1.0;
					int maxvisit=0;
					// This section finds the clustering weight of each edge in the graph based on how many neighbors and neighbors of neighbors their nodes have
					for (int i = 0; i < msz; i++) {
						int edge1 = mstedges[i].first;
						int edge2 = mstedges[i].second;
						for (int j = 0; j < edfnl.size(); j++) {
							int ed0 = edfnl[j].edge.first;
							int ed1 = edfnl[j].edge.second;
							if (edfnl[j].cost >= maxcost) {
								maxcost = edfnl[j].cost;
							}
							if (edfnl[j].visited >= maxvisit) {
								maxvisit = edfnl[j].visited;
							}
							if ((ed0 == edge1 && ed1 == edge2) || (ed1 == edge2 && ed0 == edge1)) {
								Real totalpenalty = 0;
								totalpenalty += ((n_nodes[ed0].size()) * (n_nodes[ed1].size()));
								for (int edt = 0; edt < n_nodes[ed0].size(); edt++) {
									totalpenalty += (n_nodes[n_nodes[ed0][edt]].size());
								}
								for (int edt = 0; edt < n_nodes[ed1].size(); edt++) {
									totalpenalty += (n_nodes[n_nodes[ed1][edt]].size());
								}
								decluster.push_back(edges(make_pair(ed0, ed1), totalpenalty, 0));
							}
						}
					}
					sort(decluster.begin(), decluster.end(),CompareDataedges);
					// Ths section sets the edge with highest clustering to the worst possible cost
					// So that in the next iteration of MST
					// This edge is avoided to build the unfolding path, creates a declustering
					for (int j = 0; j < edfnl.size(); j++) {
						int ed0 = edfnl[j].edge.first;
						int ed1 = edfnl[j].edge.second;
						if ((ed0 == decluster[decluster.size()-1].edge.first && ed1 == decluster[decluster.size() - 1].edge.second) || (ed1 == decluster[decluster.size() - 1].edge.second && ed0 == decluster[decluster.size() - 1].edge.first)) {
							edfnl[j].cost = maxcost;
							edfnl[j].visited = maxvisit;
						}
					}
				}
			}
		}
	}
	vector<int> UnfoldingPath::edgestoTris(vector<pair<int, int>>& ME) {
		vector<int>orderedtrngls;
		queue<int>odnodes;
		map<int, bool>nodedone;
		bool found = 0;
		plotableedges.clear();
		int insertedges=ME.size();

		// Finds which triangle to start the unfolding from
		vector<int>Tally;
		Tally.resize(insertedges + 1, 0);
		for (int i = 0; i < insertedges; i++) {
			Tally[ME[i].first]++;
			Tally[ME[i].second]++;
		}
		plotableedges[0] = -1;
		for (int i = 0; i < insertedges + 1; i++) {
			if (Tally[i] == 1) {
				odnodes.push(i);
				nodedone[i] = 1;
				found = true;
				break;
			}
		}
		if (found == 0) {
			odnodes.push(ME[0].first);
			nodedone[ME[0].first] = 1;
		}

		// Creates the order of unfolding using adjacent edges from the first triangle and keeps building the path continuously in the same way
		parenttriangle.resize(insertedges + 1, -1); 
		int counttri = 1;
		while (!odnodes.empty()) {
			int now = odnodes.front();
			orderedtrngls.push_back(now);
			odnodes.pop();
			for (int in = 0; in < insertedges; in++) {
				if (ME[in].first == now) {
					if (nodedone[ME[in].second] != 1)
					{
						parenttriangle[counttri] = find_ind(orderedtrngls, now);
						plotableedges[ME[in].second] = now;
						odnodes.push(ME[in].second);
						counttri++;
						nodedone[ME[in].second] = 1;
					}
				}
				else if (ME[in].second == now) {
					if (nodedone[ME[in].first] != 1)
					{
						parenttriangle[counttri] = find_ind(orderedtrngls, now);
						plotableedges[ME[in].first] = now;
						odnodes.push(ME[in].first);
						counttri++;
						nodedone[ME[in].first] = 1;
					}
				}
			}
		}
		return orderedtrngls;
	}
	Real UnfoldingPath::buildfinalTriangles(vector<int>&OE,Real metric3) {
		int npt = _tri._partpos.n_part;
		int overlap = 0;
		vector<vector<int>>plottedpts;
		vector<Vector3d>points;
		vector<Vector3d>centroids;
		vector < vector<int>>con2D;
		vector<Vector3d>pos2D_inner;
		vector<vector<int>>refoldingdata;
		vector<int>doneplot;
		doneplot.resize(npt, false);
		con2D.resize(npt, vector<int>());
		pos2D_inner.resize(npt, { 0,0,0 });
		
		// This section keeps track of repeated vertices
		int extct = npt;
		map<int, int>extrapts;
		for (int i = 0; i < npt; i++) {
			extrapts[i] = i;
		}
		map<int, vector<int>>plotedposes;
		vector<pair<int, int>>doneedge;
		map<int, int>plottedtemp;
		vector<vector<int>>rep_pts_inner;
		for (int i = 0; i < npt; i++) {
			rep_pts_inner.push_back(vector<int>());
		}
		vector<vector<Vector3d>>tripoints;

		// Based on the order, this section prints the triangle vertices
		for (int i = 0; i < t_size; i++) {
			tripoints.push_back(vector<Vector3d>());
			int n1 = orderedtriangles[i][0];
			int n2 = orderedtriangles[i][1];
			int n3 = orderedtriangles[i][2];
			pair<int, int>path1 = make_pair(n1, n2);
			pair<int, int>path2 = make_pair(n1, n3);
			pair<int, int>path3 = make_pair(n2, n3);
			int o1, o2, o3, tempo3, older;
			Vector3d origin;
			Vector3d origin2;
			int oldertemp = -5;
			vector<int>testnew;
			if (i == 0) {
				testnew.push_back(n1);
				testnew.push_back(n2);
				testnew.push_back(n3);
			}
			else {
				vector<int>tempnow;
				int prevtemp = plotableedges[OE[i]];
				tempnow = plotedposes[prevtemp];
				tempnow[0] = tempnow[0] % npt;
				tempnow[1] = tempnow[1] % npt;
				tempnow[2] = tempnow[2] % npt;
				if (vector_search(tempnow, n1) == 0) {
					int tempsearch = n1;
					while (plottedtemp[tempsearch] == 5) {
						tempsearch += npt;
					}
					if (vector_search(orderedtriangles[i], tempnow[0]) == 0) {
						testnew.push_back(plotedposes[prevtemp][1]);
						testnew.push_back(plotedposes[prevtemp][2]);
						oldertemp = plotedposes[prevtemp][0];
					}
					else if (vector_search(orderedtriangles[i], tempnow[1]) == 0) {
						testnew.push_back(plotedposes[prevtemp][0]);
						testnew.push_back(plotedposes[prevtemp][2]);
						oldertemp = plotedposes[prevtemp][1];
					}
					else if (vector_search(orderedtriangles[i], tempnow[2]) == 0) {
						testnew.push_back(plotedposes[prevtemp][0]);
						testnew.push_back(plotedposes[prevtemp][1]);
						oldertemp = plotedposes[prevtemp][2];
					}
					testnew.push_back(tempsearch);
				}
				if (vector_search(tempnow, n2) == 0) {
					int tempsearch = n2;
					while (plottedtemp[tempsearch] == 5) {
						tempsearch += npt;
					}
					if (vector_search(orderedtriangles[i], tempnow[0]) == 0) {
						testnew.push_back(plotedposes[prevtemp][1]);
						testnew.push_back(plotedposes[prevtemp][2]);
						oldertemp = plotedposes[prevtemp][0];
					}
					else if (vector_search(orderedtriangles[i], tempnow[1]) == 0) {
						testnew.push_back(plotedposes[prevtemp][0]);
						testnew.push_back(plotedposes[prevtemp][2]);
						oldertemp = plotedposes[prevtemp][1];
					}
					else if (vector_search(orderedtriangles[i], tempnow[2]) == 0) {
						testnew.push_back(plotedposes[prevtemp][0]);
						testnew.push_back(plotedposes[prevtemp][1]);
						oldertemp = plotedposes[prevtemp][2];
					}
					testnew.push_back(tempsearch);
				}
				if (vector_search(tempnow, n3) == 0) {
					int tempsearch = n3;
					while (plottedtemp[tempsearch] == 5) {
						tempsearch += npt;
					}
					if (vector_search(orderedtriangles[i], tempnow[0]) == 0) {
						testnew.push_back(plotedposes[prevtemp][1]);
						testnew.push_back(plotedposes[prevtemp][2]);
						oldertemp = plotedposes[prevtemp][0];
					}
					else if (vector_search(orderedtriangles[i], tempnow[1]) == 0) {
						testnew.push_back(plotedposes[prevtemp][0]);
						testnew.push_back(plotedposes[prevtemp][2]);
						oldertemp = plotedposes[prevtemp][1];
					}
					else if (vector_search(orderedtriangles[i], tempnow[2]) == 0) {
						testnew.push_back(plotedposes[prevtemp][0]);
						testnew.push_back(plotedposes[prevtemp][1]);
						oldertemp = plotedposes[prevtemp][2];
					}
					testnew.push_back(tempsearch);
				}
			}
			plotedposes[OE[i]] = testnew;
			plottedtemp[testnew[0]] = 5;
			plottedtemp[testnew[1]] = 5;
			plottedtemp[testnew[2]] = 5;
			if (i == 0) {
				origin = { 0,0,0 };
				Real d = dist(_tri._partpos.pos[n1], _tri._partpos.pos[n2]);
				origin2 = { d,0,0 };
				o1 = n1, o2=n2, o3=n3;
			}
			else {
				bool found = 0;
				for (int repi = 0; repi < rep_pts_inner[n1].size(); repi++) {
					int now = rep_pts_inner[n1][repi];
					for (int repi2 = 0; repi2 < rep_pts_inner[n2].size(); repi2++) {
						int now2 = rep_pts_inner[n2][repi2];
						if (vector_search(con2D[extrapts[now]], now2) == 1) {
							o3 = n3;
							o1 = now;
							o2 = now2;
							found = 1;
							break;
						}
					}
					if (found != 1) {
						for (int repi2 = 0; repi2 < rep_pts_inner[n3].size(); repi2++) {
							int now2 = rep_pts_inner[n3][repi2];
							if (vector_search(con2D[extrapts[now]], now2) == 1) {
								o3 = n2;
								o1 = now;
								o2 = now2;
								found = 1;
								break;
							}
						}
					}
				}
				if (found != 1) {
					for (int repi = 0; repi < rep_pts_inner[n2].size(); repi++) {
						int now = rep_pts_inner[n2][repi];
						for (int repi2 = 0; repi2 < rep_pts_inner[n3].size(); repi2++) {
							int now2 = rep_pts_inner[n3][repi2];
							if (vector_search(con2D[extrapts[now]], now2) == 1) {
								o3 = n1;
								o1 = now;
								o2 = now2;
								found = 1;
								break;
							}
						}
					}
				}
				int consize = con2D[extrapts[o1]].size();
				for (int j = 0; j < consize; j++) {
					int nowcon = con2D[extrapts[o1]][j];
					if (vector_search(con2D[extrapts[o2]], nowcon) == 1) {
						older = nowcon;
					}
				}
				tempo3 = o3;
				if (rep_pts_inner[o3 % npt].size() > 0)o3 = rep_pts_inner[tempo3][rep_pts_inner[tempo3].size() - 1] + npt;
			}
			if (i == 0)older = -5;
			o1 = testnew[0];
			o2 = testnew[1];
			o3 = testnew[2];
			older = oldertemp;
			olderverts.push_back(older);
			refoldingdata.push_back({ o1,o2,o3,older });
			Vector3d oldpos;
			if (i == 0) {
				pos2D_inner[o1] = origin;
				pos2D_inner[o2] = origin2;
				rep_pts_inner[o1 % npt].push_back(o1);
				rep_pts_inner[o2 % npt].push_back(o2);
				points.push_back(origin);
				points.push_back(origin2);
			}
			Vector3d rootpos = pos2D_inner[extrapts[o1]];
			Vector3d leverpos = pos2D_inner[extrapts[o2]];
			tripoints[i].push_back(rootpos);
			tripoints[i].push_back(leverpos);
			if (i != 0)oldpos = pos2D_inner[extrapts[older]];
			Real dis1 = dist(_tri._partpos.pos[o1 % npt], _tri._partpos.pos[o3 % npt]);
			Real dis2 = dist(_tri._partpos.pos[o2 % npt], _tri._partpos.pos[o3 % npt]);
			Real dism = dist(_tri._partpos.pos[o1 % npt], _tri._partpos.pos[o2 % npt]);
			Real a = (dis1 * dis1 - dis2 * dis2 + dism * dism) / (2 * dism);
			Real h = sqrt(dis1 * dis1 - a * a);
			Vector3d PP0 = { rootpos[0] + a * (leverpos[0] - rootpos[0]) / dism,rootpos[1] + a * (leverpos[1] - rootpos[1]) / dism,0 };
			Vector3d PP1 = { PP0[0] + h * (leverpos[1] - rootpos[1]) / dism,PP0[1] - h * (leverpos[0] - rootpos[0]) / dism ,0 };
			Vector3d PP2 = { PP0[0] - h * (leverpos[1] - rootpos[1]) / dism,PP0[1] + h * (leverpos[0] - rootpos[0]) / dism ,0 };
			int indexnow = o3;
			int prevsize = rep_pts_inner[o3 % npt].size();
			if (doneplot[o3 % npt] == 1) {
				indexnow = rep_pts_inner[indexnow % npt][prevsize - 1] + npt;
			}
			rep_pts_inner[indexnow % npt].push_back(indexnow);
			int extctadd = 0;
			if (i == 0) {
				if (PP1[1] < 0) {
					if (indexnow >= npt) {
						pos2D_inner.push_back(PP2);
						extrapts[indexnow] = extct;
						extctadd++;
						points.push_back(PP2);
						centroids.push_back((origin + origin2 + PP2) / 3);
					}
					else {
						pos2D_inner[indexnow] = PP2;
						points.push_back(PP2);
						centroids.push_back((origin + origin2 + PP2) / 3);
					}
				}
				else {
					if (indexnow >= npt) {
						pos2D_inner.push_back(PP1);
						extrapts[indexnow] = extct;
						extctadd++;
						points.push_back(PP1);
						centroids.push_back((origin + origin2 + PP1) / 3);
					}
					else {
						pos2D_inner[indexnow] = PP1;
						points.push_back(PP1);
						centroids.push_back((origin + origin2 + PP1) / 3);
					}
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
					if (indexnow >= npt) {
						pos2D_inner.push_back(PP2);
						extrapts[indexnow] = extct;
						extctadd++;
						points.push_back(PP2);
						centroids.push_back((origin + origin2 + PP2) / 3);
					}
					else {
						pos2D_inner[indexnow] = PP2;
						points.push_back(PP2);
						centroids.push_back((origin + origin2 + PP2) / 3);
					}
				}
				else {
					if (indexnow >= npt) {
						pos2D_inner.push_back(PP1);
						extrapts[indexnow] = extct;
						extctadd++;
						points.push_back(PP1);
						centroids.push_back((origin + origin2 + PP1) / 3);
					}
					else {
						pos2D_inner[indexnow] = PP1;
						points.push_back(PP1);
						centroids.push_back((origin + origin2 + PP1) / 3);
					}
				}
			}
			doneplot[n1]++;
			doneplot[n2]++;
			doneplot[n3]++;
			con2D[extrapts[o1]].push_back(o2);
			con2D[extrapts[o1]].push_back(o3);
			con2D[extrapts[o2]].push_back(o1);
			con2D[extrapts[o2]].push_back(o3);
			if (o3 >= npt) {
				extrapts[o3] = extct;
				con2D.push_back(vector<int>());
				con2D[extct].push_back(o1);
				con2D[extct].push_back(o2);
				extctadd++;
			}
			else {
				con2D[o3].push_back(o1);
				con2D[o3].push_back(o2);
			}
			if (extctadd != 0) {
				extct++;
			}
			plottedpts.push_back({ o1,o2,indexnow });
			doneedge.push_back(make_pair(extrapts[o1], extrapts[o2]));
			doneedge.push_back(make_pair(extrapts[o1], extrapts[indexnow]));
			doneedge.push_back(make_pair(extrapts[o2], extrapts[indexnow]));
			Vector3d pointnow = points[points.size() - 1];
			tripoints[i].push_back(pointnow);
			// This section finds if there is any overlap in the printed set of triangles 
			if (i != 0) {
				for (int k = 0; k < i; k++) {
					float MainT = trianglearea(tripoints[k][0], tripoints[k][1], tripoints[k][2]);
					float T1 = trianglearea(pointnow, tripoints[k][0], tripoints[k][1]);
					float T2 = trianglearea(pointnow, tripoints[k][1], tripoints[k][2]);
					float T3 = trianglearea(pointnow, tripoints[k][2], tripoints[k][0]);
					if (MainT == T1 + T2 + T3) {
						overlap = 1;
						//cout << "THERE IS OVERLAP" << endl;
					}
				}
			}
		}

		// This section keep tracks of which lines to join in terms of nodes to visualize the unfolding path
		vector<vector<int>>lines;
		for (int i = 0; i < doneedge.size(); i++) {
			vector<int>linenow;
			linenow.push_back(doneedge[i].first);
			linenow.push_back(doneedge[i].second);
			sort(linenow.begin(), linenow.end());
			bool isthere = 0;
			for (int j = 0; j < lines.size(); j++) {
				if (lines[j] == linenow) {
					isthere = 1;
				}
			}
			if (isthere == 0) {
				lines.push_back(linenow);
			}
		}

		// This section transforms the set of points to x and y direction and subsequently calculates the new weight
		vector<Vector3d>pointscopy;
		pointscopy = points;
		points.clear();
		points = centroids;
		int ps = points.size();
		Real avgx, avgy;
		Real totalx = 0, totaly = 0;
		for (int i = 0; i < ps; i++) {
			totalx += points[i][0];
			totaly += points[i][1];
		}
		avgx = totalx / ps;
		avgy = totaly / ps;
		for (int i = 0; i < ps; i++) {
			points[i][0] -= avgx;
			points[i][1] -= avgy;
		}
		Real avgx1, avgy1;
		Real totalx1 = 0, totaly1 = 0;
		for (int i = 0; i < ps; i++) {
			totalx1 += points[i][0];
			totaly1 += points[i][1];
		}
		avgx1 = totalx1 / ps;
		avgy1 = totaly1 / ps;
		Real xx = 0, xy = 0, yy = 0;
		for (int i = 0; i < ps; i++) {
			xx += pow(points[i][0] - avgx1, 2);
			yy += pow(points[i][1] - avgy1, 2);
			xy += (points[i][0] - avgx1) * (points[i][1] - avgy1);
		}
		Matrix2d M;
		M(0, 0) = xx / (ps - 1);
		M(0, 1) = xy / (ps - 1);
		M(1, 0) = xy / (ps - 1);
		M(1, 1) = yy / (ps - 1);
		Real a = M(0, 0), b = M(0, 1), c = b, d = M(1, 1);
		Real Tr = a + d;
		Real Det = a * d - b * c;
		Real L1 = Tr / (2.0) + sqrt((Tr * Tr / (4.0)) - Det);
		Real L2 = Tr / (2.0) - sqrt((Tr * Tr / (4.0)) - Det);
		Vector2d A;
		EigenSolver<Matrix2d>s(M);
		VectorXcd v = s.eigenvectors().col(1);
		Vector3d side2 = { v[0].real(),v[1].real(),0 };
		Vector3d side = { 1,0,0 };
		Vector3d signcheck = side2.cross(side) / side2.cross(side).norm();
		Real t1 = signcheck[2] * acos(side2.dot(side) / (side2.norm() * side.norm()));
		Real x1 = sin(t1 / 2) * 0;
		Real y1 = sin(t1 / 2) * 0;
		Real z1 = sin(t1 / 2) * 1;
		Real s1 = cos(t1 / 2);
		Matrix3d Q1;
		Q1 << 1 - 2 * y1 * y1 - 2 * z1 * z1, 2 * x1 * y1 - 2 * s1 * z1, 2 * x1 * z1 + 2 * s1 * y1,
			2 * x1 * y1 + 2 * s1 * z1, 1 - 2 * x1 * x1 - 2 * z1 * z1, 2 * y1 * z1 - 2 * s1 * x1,
			2 * x1 * z1 - 2 * s1 * y1, 2 * y1 * z1 + 2 * s1 * x1, 1 - 2 * x1 * x1 - 2 * y1 * y1;
		int pspo = pointscopy.size();
		for (int i = 0; i < ps; i++) {
			points[i] = Q1 * points[i];
		}
		for (int i = 0; i < pspo; i++) {
			pointscopy[i] = Q1 * pointscopy[i];
		}
		ps = points.size();
		points.clear();
		points = pointscopy;
		Real areanow;
		Real maxy = -999.0;
		Real maxx = -999.0;
		Real miny = 999.0;
		Real minx = 999.0;
		for (int i = 0; i < ps; i++) {
			if (points[i][0] < minx)minx = points[i][0];
			if (points[i][1] < miny)miny = points[i][1];
			if (points[i][0] > maxx)maxx = points[i][0];
			if (points[i][1] > maxy)maxy = points[i][1];
		}
		areanow = (maxx - minx) * (maxy - miny);
		//Real areaconst = 34.461;					// Surface area of icosahedron
		Real areaconst = 25.13;					// Surface area of dome
		//Real areaconst = 61.27;						// Surface area of parabola
		Real metric1 = 1.0-(areaconst / areanow);	// Metric for minimum bounding box
		//Real metric2 = abs((maxx - minx) - (maxy - miny)) / (20 * 2.0 - sqrt(3.0));	// Metric for minimum x-y span (normalized)
		Real metric2 = abs((maxx - minx) - (maxy - miny));
		Real finalmetric = metric3 + 50.0 * overlap;	// Assign which weights to use to find final sorted order
		//cout << metric3 << " " << metric1 << endl;
		//This sections stores the best unfolding path and prints path for visualization 
		if (finalmetric <= minmetric) {
			plottedfinalpts.clear();
			plottedfinalpts = plottedpts;
			finalparenttriangle = parenttriangle;
			minmetric = finalmetric;
			ofstream ff;
			ff.open(_partpos._folder + "/Best"+_partpos.input_particlefile);
			ff << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS ";
			ff << pos2D_inner.size();
			ff << " float\n";
			for (int r = 0; r < pos2D_inner.size(); r++) {
				ff << pos2D_inner[r][0] << " " << pos2D_inner[r][1] << " " << 0 << endl;
			}
			ff << endl;
			ff << "LINES ";
			ff << lines.size() << " " << 3 * lines.size() << endl;
			for (int i = 0; i < lines.size(); i++) {
				ff << 2 << " " << lines[i][0] << " " << lines[i][1] << endl;
			}
			ff << endl;
			ff.close();

			ofstream ft;
			ft.open(_partpos._folder + "/Besttris"+_partpos.input_connectivityfile);
			for (int i = 0; i < t_size; i++) {
				ft << orderedtriangles[i][0] << " " << orderedtriangles[i][1] << " " << orderedtriangles[i][2] << endl;
			}
			ft.close();
			ofstream fff;
			fff.open(_partpos._folder + "/BestTrisOlder" + _partpos.input_particlefile);
			for (int pp = 0; pp < refoldingdata.size(); pp++) {
				for (int ppp = 0; ppp < 4; ppp++) {
					fff << refoldingdata[pp][ppp] << " ";
				}
				fff << endl;
			}
			fff.close();

		}
		Real costfinal = finalmetric;
		return costfinal;
	}
}