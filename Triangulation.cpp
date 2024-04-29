#include "Triangulation.h"
namespace unfolding {
	void Triangulation::buildTriangles() {
		if (_inputtri != "NULL") {
			string _text;
			ifstream readFile(_inputtri);
			int indexnow = 0;
			while (getline(readFile, _text)) {
				triangles.push_back(vector<int>());
				string word;
				istringstream ss(_text);
				while (ss >> word)
				{
					triangles[indexnow].push_back(stod(word));
				}
				indexnow++;
			}
			readFile.close();
		}
		else {
			vector<vector<int>>ordered_neighbor;    // ordered list of neighboring defects
			for (int i = 0; i < _defectpos.n_def; i++) {
				ordered_neighbor.push_back(vector<int>());
				int d_size = _defectpos.defect_neighbor[i].size();
				if (d_size == 3) {

				}
				vector<vector<int>>temp;
				for (int j = 0; j < d_size; j++) {
					int def_j = _defectpos.defect_neighbor[i][j];
					temp.push_back(vector<int>());
					priority_queue<pair<int, int>>q;
					for (int k = 0; k < d_size; k++) {
						if (k != j) {
							int def_k = _defectpos.defect_neighbor[i][k];
							pair<int, int>p;
							p = make_pair(_partpos.MinDisGrp[def_j][def_k], def_k);
							q.push(p);
						}
					}
					while (q.size() > 2) {
						q.pop();
					}

					pair<int, int>p = q.top();
					temp[j].push_back(p.second);
					q.pop();
					temp[j].push_back(def_j);
					pair<int, int>p1 = q.top();
					temp[j].push_back(p1.second);
				}
				Vector3d AB = _partpos.pos[temp[0][0]] - _partpos.pos[temp[0][1]];
				Vector3d BC = _partpos.pos[temp[0][2]] - _partpos.pos[temp[0][1]];
				Vector3d Cr = AB.cross(BC);
				Vector3d n = Cr / Cr.norm();
				int countside = 0;
				for (int k = 0; k < _partpos.n_part; k++) {
					Vector3d P = _partpos.pos[k];
					Vector3d x = P - _partpos.pos[temp[0][0]];
					Real Crr = n.dot(x);
					if (Crr > 0)countside++;
				}
				if (countside < ceil(_partpos.n_part / 2.0)) {
					ordered_neighbor[i].push_back(temp[0][2]);
					ordered_neighbor[i].push_back(temp[0][1]);
					ordered_neighbor[i].push_back(temp[0][0]);
				}
				else {
					ordered_neighbor[i].push_back(temp[0][0]);
					ordered_neighbor[i].push_back(temp[0][1]);
					ordered_neighbor[i].push_back(temp[0][2]);
				}
				int current = ordered_neighbor[i][2];
				int previous = ordered_neighbor[i][1];
				while (ordered_neighbor[i].size() < d_size) {
					for (int l = 0; l < d_size; l++) {
						if (current == temp[l][1]) {
							if (previous == temp[l][0]) {
								ordered_neighbor[i].push_back(temp[l][2]);
								previous = current;
								current = temp[l][2];
								break;
							}
							else {
								ordered_neighbor[i].push_back(temp[l][0]);
								previous = current;
								current = temp[l][0];
								break;
							}
						}
					}
				}
			}
			for (int i = 0; i < _defectpos.n_def; i++) {
				int def_i = _defectpos.defect_pos[i];
				for (int j = 0; j < ordered_neighbor[i].size(); j++) {
					int t_pos1 = j;
					int t_pos2 = (j + 1) % ordered_neighbor[i].size();
					vector<int>trianglenow = { def_i,ordered_neighbor[i][t_pos1],ordered_neighbor[i][t_pos2] };
					sort(trianglenow.begin(), trianglenow.end());
					if (triangles.size() == 0)triangles.push_back(trianglenow);
					bool found = 0;
					for (int k = 0; k < triangles.size(); k++) {
						if (triangles[k] == trianglenow) {
							found = 1;
						}
					}
					if (found == 0) {
						triangles.push_back(trianglenow);
					}
				}
			}
		}
	}
	void Triangulation::getTriangles() {
		int t_size = triangles.size();
		cout << "Printing unfoldable triangles' position of vertices: " << endl;
		for (int i = 0; i < t_size; i++) {
			cout << i << ": ";
			print_vector(triangles[i]);
		}
	}
	void Triangulation::shuffleTriangles() {
		random_shuffle(triangles.begin(), triangles.end());
	}
}