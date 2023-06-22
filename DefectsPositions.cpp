#include "DefectsPositions.h"
namespace unfolding {
	void DefectsPositions::findDefects() {
		if (_inputdefectfile != "NULL") {
			string _text;
			ifstream readFile(_inputdefectfile);
			while (getline(readFile, _text)) {
				string word;
				istringstream s(_text);
				while (s >> word)
				{
					defect_pos.push_back(stod(word));
				}
			}
			readFile.close();
		}
		else {
			for (int i = 0; i < _partpos.n_part; i++) {
				int isize = _partpos.pos_neighbor[i].size();
				if (_partpos.pos_neighbor[i].size() != 6) {
					defect_pos.push_back(i);
				}
			}
		}
		n_def = defect_pos.size();
	}
	void DefectsPositions::findConnectivity() {
		vector<int>def_neighbor;
		for (int i = 0; i < n_def; i++) {
			int def_i = defect_pos[i];
			int isstop = 0, isconst = 0, init_r = _initialr, lastcount_n = 0, maxn = -10, maxr;
			vector<int>list_p;
			vector<int>list_n;
			map<Real, Real>constcount;
			map<Real, Real>trackcount;	
			while (isstop != 2) {
				Real count_n = 0;
				for (int j = 0; j < n_def; j++) {
					int def_j = defect_pos[j];
					Real r = _partpos.MinDisGrp[def_i][def_j];
					if (r <= init_r + _delr / 2 && r >= init_r - _delr / 2 && i != j) {
						count_n++;
					}
				}
				int countdist = count_n - maxn;
				lastcount_n = count_n;
				if (maxn < count_n) {
					maxn = count_n;
					maxr = init_r;
				}
				list_p.push_back(count_n);
				list_n.push_back(init_r);
				if (init_r > _rc || count_n > 15 || countdist < -2) {
					isstop = 2;
				}
				init_r += _dr;
			}
			def_neighbor.push_back(maxr);
		}
		defect_neighbor.resize(n_def, vector<int>());
		for (int i = 0; i < n_def; i++) {
			int def_i = defect_pos[i];
			int range = def_neighbor.at(i);
			for (int j = 0; j < n_def; j++) {
				int def_j = defect_pos[j];
				int r = _partpos.MinDisGrp[def_i][def_j];
				if (r <= range + _delr / 2 && i != j) {
					int isize = defect_neighbor[i].size();
					int jsize = defect_neighbor[j].size();
					// If 'i' contains 'j' then the vice versa should be true, made sure in the next part of the code
					int ci = 0, cj = 0;
					for (int k = 0; k < isize; k++) {
						if (defect_neighbor[i].at(k) == def_j)ci++;
					}
					for (int k = 0; k < jsize; k++) {
						if (defect_neighbor[j].at(k) == def_i)cj++;
					}
					if (ci == 0)defect_neighbor[i].push_back(def_j);
					if (cj == 0)defect_neighbor[j].push_back(def_i);
				}
			}
		}
	}
	void DefectsPositions::mapdefectid() {
		for (int i = 0; i < n_def; i++) {
			id_defects[defect_pos[i]] = i;
		}
	}
	int DefectsPositions::getNofDefects() {
		return n_def;
	}
	void DefectsPositions::getDefectsLocations() {
		cout << "Locations of defects are: " << endl;
		print_vector(defect_pos);
	}
	void DefectsPositions::getAllDefectsNeighbors() {
		for (int i = 0; i < n_def; i++) {
			cout << "Neighbors of defect at " << defect_pos[i] << " are: " << endl;
			print_vector(defect_neighbor[i]);
		}
	}
	void DefectsPositions::getDefectsNeighborsID(int i) {
			cout << "Neighbors of defect at " << i << " are: " << endl;
			print_vector(defect_neighbor[id_defects[i]]);
	}
}