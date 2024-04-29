#ifndef __UNFOLDINGPATH_H__
#define __UNFOLDINGPATH_H__
#include "Triangulation.h"
namespace unfolding {
	class UnfoldingPath {
	public:
		Triangulation _tri;
		ParticlePositions _partpos;
		Real minarea;
		int t_size;
		vector<vector<int>>triconn;					//AdjacencyList of connectivty among triangles
		vector<edges>edfnl;
		vector<edges>ed;							//List of unfoldable edges
		int _mststeps;
		int edsize;									//Number of edges
		map<int, int>plotableedges;					//Need to make it local inside cpp file
		vector<vector<int>>plottedpts;
		vector<vector<int>>plottedfinalpts;			//Plotted triangle locations
		vector<int>finalparenttriangle;
		vector<vector<int>>orderedtriangles;		//Triangle correct order to plot
		map<int, int>previous_triangle;
		vector<int>parenttriangle;
		vector<int>olderverts;
		int randlimit;
		int bestcount;
		UnfoldingPath(Triangulation TR, ParticlePositions PartPos, int MSTSteps) :_tri(TR), _partpos(PartPos), _mststeps(MSTSteps) {
			t_size = _tri.triangles.size();
			buildAdjacency();
			//srand(time(NULL));
			buildEdges();
			minarea = 99999.0;
			randlimit = 1500;
			bestcount = 0;
			runMST(_mststeps);
		}

		void buildAdjacency();
		void buildEdges();
		void printAdjacency();
		void runMST(int steps);
		void sortEdges();
		void sortEdgesfnl();
		vector<int> edgestoTris(vector<pair<int, int>>& ME);
		Real buildfinalTriangles(vector<int>&OE);
	};
}
#endif 