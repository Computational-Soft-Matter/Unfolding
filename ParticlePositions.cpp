#include "ParticlePositions.h"

namespace unfolding {
	void ParticlePositions::readParticlePositions() {
		string _text;
		ifstream _readFile(input_particlefile);
		while (getline(_readFile, _text)) {
			string word;
			istringstream s(_text);
			Real temp_pos[3] = {};
			int indexnow = 0;
			while (s >> word)
			{
				temp_pos[indexnow] = stof(word);
				indexnow++;
				if (indexnow == 3) {
					pos.push_back({ temp_pos[0], temp_pos[1], temp_pos[2] });
					indexnow = 0;
				}
			}
		}
		_readFile.close();
	}

	void ParticlePositions::readParticleConnectivity() {
		string _text;
		ifstream _readFile(input_connectivityfile);
		int indexnow = 0;
		while (getline(_readFile, _text)) {
			pos_neighbor.push_back(vector<int>());
			string word;
			istringstream s(_text);
			while (s >> word)
			{
				pos_neighbor[indexnow].push_back(stod(word));
			}
			indexnow++;
		}
		_readFile.close();
	}

	void ParticlePositions::getAllParticles() {
		for (int i = 0; i < n_part; i++) {
			cout << "Printing position of particle: " << i <<endl;
			cout << pos[i] << endl;
		}
	}

	void ParticlePositions::getParticleID(int i) {
		cout << "Printing position of particle: " << i <<endl;
		cout << pos[i] << endl;
	}

	void ParticlePositions::getParticleNeighbors(int i) {
		cout << "Printing neighboring particles of particle: " << i <<endl;
		print_vector(pos_neighbor[i]);
	}

	void ParticlePositions::calcFloyd() {
		for (int i = 0; i < n_part; i++) {
			for (int j = 0; j < n_part; j++) {
				MinDisGrp[i][j] = INF;
			}
		}
		for (int i = 0; i < n_part; i++) {
			MinDisGrp[i][i] = 0;
			for (int j = 0; j < pos_neighbor[i].size(); j++) {
				MinDisGrp[i][pos_neighbor[i][j]] = 1;
			}
		}
		for (int k = 0; k < n_part; k++) {
			for (int i = 0; i < n_part; i++) {
				for (int j = 0; j < n_part; j++) {
					if (MinDisGrp[i][k] + MinDisGrp[k][j] < MinDisGrp[i][j])
						MinDisGrp[i][j] = MinDisGrp[i][k] + MinDisGrp[k][j];
				}
			}
		}
	}
}