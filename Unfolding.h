#ifndef _UNFOLDING_H_
#define _UNFOLDING_H_

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<algorithm>
#include<array>
#include<vector>
#include<Eigen>
#include<math.h>
#include<stdlib.h>
#include<queue>
#include<map>
#include<utility>
#include<time.h>
#include<ctime>
#define INF 999
#define PI 3.14159265

using namespace std;
using namespace Eigen;

namespace unfolding {

	typedef double Real;
	void print_vector(vector<Real>& A);
	void print_vector(vector<int>& A);
	bool vector_search(vector<int>&A, int a);
	class edges
	{
	public:
		pair<int, int>edge;
		Real cost;
		int visited;
		edges(pair<int, int>edge, Real cost, int visited) {
			this->edge = edge;
			this->cost = cost;
			this->visited = visited;
		}
	};
	Real dist(Vector3d A, Vector3d B);
	Real dist(Vector2d A, Vector2d B);
	Real trianglearea(Vector3d A, Vector3d B, Vector3d C);
	class TriPtsDS {
	public:
		vector<int>Vrt;
		vector<Vector3d>Plts;
	};
}

#endif // _UNDOLFING_H_
