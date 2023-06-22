#ifndef __UNFOLDINGPATH_H__
#define __UNFOLDINGPATH_H__
#include "Triangulation.h"
namespace unfolding {
	class UnfoldingPath {
	public:
		Triangulation _tri;
		int t_size;
		vector<vector<int>>triconn;					//AdjacencyList of connectivty among triangles
		vector<edges>ed;							//List of unfoldable edges
		int _mststeps;
		int edsize;									//Number of edges
		map<int, int>plotableedges;					//Need to make it local inside cpp file
		vector<vector<int>>plottedfinalpts;			//Plotted triangle locations
		vector<vector<int>>orderedtriangles;		//Triangle correct order to plot
		vector<int>olderverts;
		UnfoldingPath(Triangulation TR, int MSTSteps) :_tri(TR), _mststeps(MSTSteps) {
			t_size = _tri.triangles.size();
			buildAdjacency();
			buildEdges();
			runMST(_mststeps);
		}

		void buildAdjacency();
		void buildEdges();
		void printAdjacency();
		void runMST(int steps);
		void sortEdges();
		vector<int> edgestoTris(vector<pair<int, int>>& ME);
		void buildfinalTriangles(vector<int>&OE);
	};
}
#endif 