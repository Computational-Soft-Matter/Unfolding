#include "Unfolding.h"
namespace unfolding {
	void print_vector(vector<Real>& A) {
		int Asize = A.size();
		for (int i = 0; i < Asize; i++) {
			cout << A.at(i) << ", ";
		}
		cout << "Reals" << endl;
	}

	void print_vector(vector<int>& A) {
		int Asize = A.size();
		for (int i = 0; i < Asize; i++) {
			cout << A.at(i) << ", ";
		}
		cout << "ints" << endl;
	}

	bool vector_search(vector<int>&A, int a) {
		int Asize = A.size();
		int counta = 0;
		for (int i = 0; i < Asize; i++) {
			if (A.at(i) == a) {
				counta++;
			}
		}
		if (counta != 0) {
			return 1;
		}
		else return 0;
	}
	int find_ind(vector<int>& A, int a) {
		int Asize = A.size();
		int counta = 0;
		for (int i = 0; i < Asize; i++) {
			if (A.at(i) == a) {
				return i;
			}
		}
		return 0;
	}
	Real dist(Vector3d A, Vector3d B) {
		return sqrt(pow(A(0) - B(0), 2) + pow(A(1) - B(1), 2) + pow(A(2) - B(2), 2));
	}

	Real dist(Vector2d A, Vector2d B) {
		return sqrt(pow(A(0) - B(0), 2) + pow(A(1) - B(1), 2));
	}

	Real trianglearea(Vector3d A, Vector3d B, Vector3d C) {
		Vector3d AB = B - A;
		Vector3d AC = C - A;
		Vector3d Cr = AB.cross(AC);
		return 0.5 * (Cr.norm());
	}
}