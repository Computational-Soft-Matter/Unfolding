#include "UnfoldTriangles.h"
namespace unfolding {
	TriPtsDS UnfoldTriangles::TriPts(vector<int>T) {
		vector<int>final_vert;
		vector<Vector3d>plots;
		queue<int>path;
		int n1 = T[0], n2 = T[1], n3 = T[2];
		final_vert.push_back(n1);
		final_vert.push_back(n2);
		final_vert.push_back(n3);
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
			Real x = _partpos.pos[vn_ngh][0];
			Real y = _partpos.pos[vn_ngh][1];
			Real z = _partpos.pos[vn_ngh][2];
			Real t = (n(0) * d - n(0) * x + n(1) * e - n(1) * y + n(2) * f - n(2) * z) / (n(0) * n(0) + n(1) * n(1) + n(2) * n(2));
			Vector3d proj = { x + t * n(0),y + t * n(1),z + t * n(2) };
			plots.push_back(_partpos.pos[vn_ngh]);
		}
		path.push(n1);
		path.push(n2);
		path.push(n3);
		vector<int>testvec;
		vector<Vector3d>plotsbackup;
		vector<bool>visitedpath;
		visitedpath.resize(_partpos.n_part, false);
		visitedpath[n1] = 1;
		visitedpath[n2] = 1;
		visitedpath[n3] = 1;
		while (!path.empty()) {
			int verticenow = path.front();
			path.pop();
			int vnsize = _partpos.pos_neighbor[verticenow].size();
			for (int k = 0; k < vnsize; k++) {
				int vn_ngh = _partpos.pos_neighbor[verticenow][k];
				//cout << verticenow << " colorado " << vn_ngh << endl;
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
						final_vert.push_back(vn_ngh);
						plots.push_back(_partpos.pos[vn_ngh]);
					}
					else {
						if (verticenow != n1 && verticenow != n2 && verticenow != n3)
						{
							visitedpath[vn_ngh] = true;
							testvec.push_back(vn_ngh);
							plotsbackup.push_back(_partpos.pos[vn_ngh]);
						}
					}
				}
			}
		}
		for (int k = 0; k < testvec.size(); k++) {
			final_vert.push_back(testvec[k]);
			plots.push_back(plotsbackup[k]);
		}
		TriPtsDS tempDS;
		tempDS.Vrt = final_vert;
		tempDS.Plts = plots;
		return tempDS;
	}
	void UnfoldTriangles::Unfold() {
		vector<pair<int, int>>doneedge;
		for (int i = 0; i < _unfpt.t_size; i++) {
			vector<int>previouscompare;
			vector<int>unfnow = _unfpt.plottedfinalpts[i];
			vector<int>unfnoworg = {unfnow[0]%_partpos.n_part, unfnow[1] % _partpos.n_part, unfnow[2] % _partpos.n_part };
			TriPtsDS DSnow = TriPts(unfnoworg);
			vector<int>final_vert = DSnow.Vrt;
			vector<Vector3d>plots = DSnow.Plts;
			int o1= unfnow[0];
			int o2 = unfnow[1];
			int o3 = unfnow[2];
			cout << o1 << " " << o2 << " " << o3 << endl;
			Vector3d origin,origin2,side,side2;
			origin = { pos2D[ID2D[o1]][0],pos2D[ID2D[o1]][1],0.0 };
			origin2 = { pos2D[ID2D[o2]][0],pos2D[ID2D[o2]][1],0.0 };
			side = { pos2D[ID2D[o2]][0] - pos2D[ID2D[o1]][0],pos2D[ID2D[o2]][1] - pos2D[ID2D[o1]][1],0 };
			vector<vector<int>>smalltris;
			for (int j = 0; j < final_vert.size(); j++) {
				int indn = final_vert[j];
				int sizenow = _partpos.pos_neighbor[indn].size();
				for (int k = 0; k < sizenow; k++) {
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
						if (pushnow == 0)smalltris.push_back(P);
					}
				}
			}
			int startvrt;
			bool foundstartvrt = 0;
			for (int j = 0; j < smalltris.size(); j++) {
				if (vector_search(smalltris[j], o1 % _partpos.n_part) == 1) {
					startvrt = j;
					foundstartvrt = 1;
				}
			}
			int closeo1;
			Real disto1 = 9999;
			if (foundstartvrt == 0) {
				for (int k = 0; k < final_vert.size(); k++) {
					if (o1 % _partpos.n_part != final_vert[k]) {
						Real disno11 = dist(_partpos.pos[o1 % _partpos.n_part], _partpos.pos[final_vert[k]]);
						if (disno11 < disto1) {
							disto1 = disno11;
							closeo1 = final_vert[k];
						}
					}
				}
				for (int j = 0; j < smalltris.size(); j++) {
					if (vector_search(smalltris[j], closeo1) == 1) {
						startvrt = j;
					}
				}
			}
			vector<vector<int>>odsmalltris;
			vector<bool>mappedtrism;
			mappedtrism.resize(smalltris.size(), false);
			queue<vector<int>>unftrsm;
			unftrsm.push(smalltris[startvrt]);
			mappedtrism[startvrt] = true;
			int totalsm = 1;
			while (!unftrsm.empty()) {
				vector<int>tnow = unftrsm.front();
				pair<int, int>pair1 = make_pair(tnow[0], tnow[1]);
				pair<int, int>pair2 = make_pair(tnow[0], tnow[2]);
				pair<int, int>pair3 = make_pair(tnow[2], tnow[1]);
				odsmalltris.push_back(tnow);
				unftrsm.pop();
				for (int i = 0; i < smalltris.size(); i++) {
					if (vector_search(smalltris[i], pair1.first) == 1 && vector_search(smalltris[i], pair1.second) == 1 && tnow != smalltris[i]) {
						if (mappedtrism[i] == false) {
							unftrsm.push(smalltris[i]);
							mappedtrism[i] = true;
							totalsm++;
						}
					}
				}
				for (int i = 0; i < smalltris.size(); i++) {
					if (vector_search(smalltris[i], pair2.first) == 1 && vector_search(smalltris[i], pair2.second) == 1 && tnow != smalltris[i]) {
						if (mappedtrism[i] == false) {
							unftrsm.push(smalltris[i]);
							mappedtrism[i] = true;
							totalsm++;
						}
					}
				}
				for (int i = 0; i < smalltris.size(); i++) {
					if (vector_search(smalltris[i], pair3.first) == 1 && vector_search(smalltris[i], pair3.second) == 1 && tnow != smalltris[i]) {
						if (mappedtrism[i] == false) {
							unftrsm.push(smalltris[i]);
							mappedtrism[i] = true;
							totalsm++;
						}
					}
				}
			}
			vector<pair<int, Vector3d>>midplots;
			vector<int>smalldone;
			Vector3d P1 = _partpos.pos[odsmalltris[0][0]];
			Vector3d P2 = _partpos.pos[odsmalltris[0][1]];
			Vector3d P3 = _partpos.pos[odsmalltris[0][2]];
			Vector3d An, Bn, Cn, Dn;
			Vector3d P12 = P2 - P1;
			Vector3d P23 = P3 - P2;
			Vector3d N0 = P12.cross(P23);
			Vector3d Nn0 = N0 / N0.norm();
			int firstnode, secondnode;
			for (int k = 0; k < 3; k++) {
				if (vector_search(_unfpt._tri._defectpos.defect_pos, odsmalltris[0][k]) == 1) {
					firstnode = odsmalltris[0][k];
					secondnode = odsmalltris[0][(k + 1) % 3];
				}
			}
			Real dis = dist(_partpos.pos[firstnode], _partpos.pos[secondnode]);
			Vector3d first = { 0,0,0 };
			Vector3d second = { dis,0,0 };
			midplots.push_back(make_pair(firstnode, first));
			midplots.push_back(make_pair(secondnode, second));
			smalldone.push_back(firstnode);
			smalldone.push_back(secondnode);
			for (int j = 0; j < odsmalltris.size(); j++) {
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
				else continue;
				int oldnode;
				for (int k = 0; k < j; k++) {
					if (k != j && vector_search(odsmalltris[k], root) == 1 && vector_search(odsmalltris[k], lever) == 1)
					{
						if (root != odsmalltris[k][0] && lever != odsmalltris[k][0])oldnode = odsmalltris[k][0];
						else if (root != odsmalltris[k][1] && lever != odsmalltris[k][1])oldnode = odsmalltris[k][1];
						else oldnode = odsmalltris[k][2];
					}
				}
				if (j == 0)oldnode = -1;
				Real dis1 = ceil(dist(_partpos.pos[root], _partpos.pos[image]) * 1000.0) / 1000.0;
				Real dis2 = ceil(dist(_partpos.pos[lever], _partpos.pos[image]) * 1000.0) / 1000.0;
				Real dism = ceil(dist(_partpos.pos[root], _partpos.pos[lever]) * 1000.0) / 1000.0;
				Vector3d rootpos,leverpos,oldpos;
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
				if (j == 0) {
					if (PP1[1] < 0)midplots.push_back(make_pair(image, PP2));
					else midplots.push_back(make_pair(image, PP1));
				}
				else if (foundstartvrt == 0 && j == 0) {
					if (PP1[1] < 0) midplots.push_back(make_pair(image, PP2));
					else midplots.push_back(make_pair(image, PP1));
				}
				else {
					Vector3d aa1, bb1, cc1, First1, Second1;
					aa1 = { leverpos[0] - rootpos[0],leverpos[1] - rootpos[1], 0 };
					bb1 = { oldpos[0] - rootpos[0],oldpos[1] - rootpos[1], 0 };
					cc1 = { PP1[0] - rootpos[0],PP1[1] - rootpos[1], 0 };
					Second1 = aa1.cross(bb1) / aa1.cross(bb1).norm();
					First1 = aa1.cross(cc1) / aa1.cross(cc1).norm();
					if (First1[2] == Second1[2]) midplots.push_back(make_pair(image, PP2));
					else midplots.push_back(make_pair(image, PP1));
				}
				smalldone.push_back(image);
			}
			bool Bnf = 0;
			for (int j = 0; j < final_vert.size(); j++) {
				int indno = final_vert[j];
				if (foundstartvrt == 0 && indno == (o1 % _partpos.n_part)) {
					Vector3d tempvec = { pos2D[indno][0],pos2D[indno][1],0 };
					plots[j] = tempvec;
					An = tempvec;
					continue;
				}
				for (int k = 0; k < midplots.size(); k++) {
					if (indno == midplots[k].first) {
						plots[j] = midplots[k].second;
						if (indno == (o1 % _partpos.n_part))An = plots[j];
						if (indno == (o2 % _partpos.n_part)) {
							Bn = plots[j];
							Bnf = 1;
						}
						if (indno == (o3 % _partpos.n_part))Cn = plots[j];
						if (foundstartvrt == 0) {
							if (indno == closeo1)Dn = plots[j];
						}
					}

				}
			}
			if (i != 0) {
				for (int k = 0; k < triptsset.size(); k++) {
					if (vector_search(triptsset[k].first, (o1 % _partpos.n_part)) == 1 && vector_search(triptsset[k].first, (o2 % _partpos.n_part)) == 1) {
						previouscompare = triptsset[k].second;
					}
				}
			}
			Vector3d A = _partpos.pos[o1 % _partpos.n_part];
			Vector3d B = _partpos.pos[o2 % _partpos.n_part];
			Vector3d C = _partpos.pos[o3 % _partpos.n_part];
			Vector3d AB = B - A;
			Vector3d AC = C - A;
			Vector3d N = AB.cross(AC);
			Vector3d Nn = N / N.norm();
			Vector3d Z = { 0,0,1 };
			Vector3d u = Nn.cross(Z);
			Vector3d U = u / u.norm();
			Real t = acos(Nn.dot(Z) / (Nn.norm() * Z.norm()));
			Matrix3d M;
			M << cos(t) + pow(U(0), 2)* (1 - cos(t)), U(0)* U(1)* (1 - cos(t)) - U(2)* sin(t), U(0)* U(1)* (1 - cos(t)) + U(2)* sin(t),
				U(1)* U(2)* (1 - cos(t)) + U(2)* sin(t), cos(t) + pow(U(1), 2)* (1 - cos(t)), U(2)* U(1)* (1 - cos(t)) - U(0)* sin(t),
				U(0)* U(2)* (1 - cos(t)) - U(1)* sin(t), U(2)* U(1)* (1 - cos(t)) + U(0)* sin(t), cos(t) + pow(U(2), 2)* (1 - cos(t));
			Real x = sin(t / 2) * U(0);
			Real y = sin(t / 2) * U(1);
			Real z = sin(t / 2) * U(2);
			Real s = cos(t / 2);
			Matrix3d Q;
			Q << 1 - 2 * y * y - 2 * z * z, 2 * x * y - 2 * s * z, 2 * x * z + 2 * s * y,
				2 * x * y + 2 * s * z, 1 - 2 * x * x - 2 * z * z, 2 * y * z - 2 * s * x,
				2 * x * z - 2 * s * y, 2 * y * z + 2 * s * x, 1 - 2 * x * x - 2 * y * y;
			Vector3d diff;
			diff = origin - An;
			An = An + diff;
			Bn = Bn + diff;
			Cn = Cn + diff;
			side2 = { Bn(0) - An(0),Bn(1) - An(1),0 };
			Vector3d side3 = -side2;
			Vector3d diff2;
			Vector3d signcheck = side2.cross(side) / side2.cross(side).norm();
			Real t1 = signcheck[2] * acos(side2.dot(side) / (side2.norm() * side.norm()));
			if (float(side[0]) == float(side3[0]) && float(side[1]) == float(side3[1])) {
				t1 = 3.1415;
			}
			Real x1 = sin(t1 / 2) * 0;
			Real y1 = sin(t1 / 2) * 0;
			Real z1 = sin(t1 / 2) * 1;
			Real s1 = cos(t1 / 2);
			Matrix3d Q1;
			Q1 << 1 - 2 * y1 * y1 - 2 * z1 * z1, 2 * x1 * y1 - 2 * s1 * z1, 2 * x1 * z1 + 2 * s1 * y1,
				2 * x1 * y1 + 2 * s1 * z1, 1 - 2 * x1 * x1 - 2 * z1 * z1, 2 * y1 * z1 - 2 * s1 * x1,
				2 * x1 * z1 - 2 * s1 * y1, 2 * y1 * z1 + 2 * s1 * x1, 1 - 2 * x1 * x1 - 2 * y1 * y1;
			vector<int>printing;
			int extnow = extid;
			for (int j = 0; j < final_vert.size(); j++) {
				if ((i == 0 || (vector_search(previouscompare, final_vert[j]) == 0))) {
					int indexnow = final_vert[j];
					int prevsize = rep_pts[indexnow].size();
					if (mapped[ID2D[indexnow]] == 1) {
						indexnow = rep_pts[indexnow][prevsize - 1] + _partpos.n_part;
					}
					rep_pts[final_vert[j]].push_back(indexnow);
					Vector3d mp = plots[j];
					if (indexnow > _partpos.n_part) {
						ID2D[indexnow] = extid;
						pos2D.push_back({ 0,0 });
						pos2D[extid][0] = mp(0) + diff(0);
						pos2D[extid][1] = mp(1) + diff(1);
						extid++;
					}
					else {
						pos2D[indexnow][0] = mp(0) + diff(0);
						pos2D[indexnow][1] = mp(1) + diff(1);
					}
					printing.push_back(indexnow);
				}
			}
			extnow = extid - extnow;
			for (int j = 0; j < extnow; j++) {
				mapped.push_back(false);
			}
			for (int j = 0; j < final_vert.size(); j++) {
				if ((i == 0 || (vector_search(previouscompare, final_vert[j]) == 0))) {
					int indexnow = final_vert[j];
					int prevsize = rep_pts[indexnow].size();
					indexnow = rep_pts[indexnow][prevsize - 1];
					Vector3d ver = { pos2D[ID2D[indexnow]][0],pos2D[ID2D[indexnow]][1],0 };
					ver[0] = ver[0] - pos2D[ID2D[o1]][0];
					ver[1] = ver[1] - pos2D[ID2D[o1]][1];
					if (side2.cross(side).norm() != 0)ver = Q1 * ver;
					pos2D[ID2D[indexnow]][0] = ver[0] + pos2D[ID2D[o1]][0];
					pos2D[ID2D[indexnow]][1] = ver[1] + pos2D[ID2D[o1]][1];
					if (i == 0)mapped[ID2D[final_vert[j]]] = true;
				}
			}
			bool changedirection = false;
			if (i != 0) {
				Vector3d aa, bb, cc, First, Second;
				aa = { pos2D[ID2D[o2]][0] - pos2D[ID2D[o1]][0],pos2D[ID2D[o2]][1] - pos2D[ID2D[o1]][1], 0 };
				bb = { pos2D[ID2D[_unfpt.olderverts[i]]][0] - pos2D[ID2D[o1]][0],pos2D[ID2D[_unfpt.olderverts[i]]][1] - pos2D[ID2D[o1]][1], 0 };
				cc = { pos2D[ID2D[o3]][0] - pos2D[ID2D[o1]][0],pos2D[ID2D[o3]][1] - pos2D[ID2D[o1]][1], 0 };
				Second = aa.cross(bb) / aa.cross(bb).norm();
				First = aa.cross(cc) / aa.cross(cc).norm();
				if (First[2] == Second[2])changedirection = true;
			}
			if (changedirection == true && i != 0) {
				for (int j = 0; j < final_vert.size(); j++) {
					if ((i == 0 || (vector_search(previouscompare, final_vert[j]) == 0))) {
						int indexnow = final_vert[j];
						int prevsize = rep_pts[indexnow].size();
						indexnow = rep_pts[indexnow][prevsize - 1];
						Vector3d ver = { pos2D[ID2D[indexnow]][0],pos2D[ID2D[indexnow]][1],0 };
						Real p = ver(0);
						Real q = ver(1);
						Real a, b, c, d;
						a = pos2D[ID2D[o1]][0];
						b = pos2D[ID2D[o1]][1];
						c = pos2D[ID2D[o2]][0];
						d = pos2D[ID2D[o2]][1];
						Matrix2d D, Dx, Dy;
						D << (a - c) / (d - b), -1,
							(d - b) / (c - a), -1;
						Dx << (a - c) * p / (d - b) - q, -1,
							(d - b) * a / (c - a) - b, -1;
						Dy << (a - c) / (d - b), (a - c)* p / (d - b) - q,
							(d - b) / (c - a), (d - b)* a / (c - a) - b;
						if (mapped[ID2D[indexnow]] == false) {					
							if ((d - b) == 0 || (c - a) == 0) {
								pos2D[ID2D[indexnow]][0] = p;
								pos2D[ID2D[indexnow]][1] = -q;
								mapped[ID2D[indexnow]] = true;
							}
							else {
								pos2D[ID2D[indexnow]][0] = 2 * Dx.determinant() / D.determinant() - p;
								pos2D[ID2D[indexnow]][1] = 2 * Dy.determinant() / D.determinant() - q;
								mapped[ID2D[indexnow]] = true;

							}
						}

					}
				}
			}
			else {
				for (int j = 0; j < printing.size(); j++) {
					mapped[ID2D[printing[j]]] = true;
				}
			}
			vector<int>ptest;
			ptest.push_back(o1 % _partpos.n_part);
			ptest.push_back(o2% _partpos.n_part);
			ptest.push_back(o3% _partpos.n_part);
			triptsset.push_back(make_pair(ptest, final_vert));
			cout << ID2D[o3] << endl;
		}
		cout << pos2D.size() << endl;
		ofstream f;
		f.open("Unfolded.vtk");
		f << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS ";
		f << pos2D.size();
		f << " float\n";
		for (int i = 0; i < pos2D.size(); i++) {
			f << pos2D[i][0] << " " << pos2D[i][1] << " " << 0 << endl;
		}
		f << endl;
		f.close();
	}
}