#ifndef __UNFOLDINGPATH_H__
#define __UNFOLDINGPATH_H__
#include "Triangulation.h"
namespace unfolding {
	class UnfoldingPath {
	public:
		Triangulation _tri;
		ParticlePositions _partpos;

		Real minmetric;							// This is used to compare the current minimum metric after every MST iterations
		int t_size;								// Numer of large triangles
		int randlimit;							// Number of iterations for random unfolding path search 
		int _mststeps;							// Number of iterations where MST algorithm finds an unfolding path
		int edsize;								// Number of edges in the large triangles graph

		vector<vector<int>>triconn;				// AdjacencyList of connectivty among triangles
		vector<edges>ed;						// Initial List of unfoldable edges
		vector<edges>edfnl;						// Final list of unfoldable edges after completing the path
		map<int, int>plotableedges;				// Maps the parent triangle node of each unfolded triangle node, this is updated with each MST iteration	
		vector<vector<int>>orderedtriangles;	// List of triangles in the order of unfolding, expressed as a set of vertices modded by the number of total points, this is useful initially to locate the unfolding path
		vector<vector<int>>plottedfinalpts;		// List of triangles in the order of unfolding, expressed as a set of vertices, every repeated vertex has a new index in this list, this list is used later in "UnfoldTriangles.cpp"
		vector<int>parenttriangle;				// List of parent triangle node index of each unfolded triangle node, this is updated with each MST iteration
		vector<int>finalparenttriangle;			// List of parent triangle node index of each unfolded triangle node for the best unfolding path, and will be later used for unfolding smaller triangles
		vector<int>olderverts;
		
		UnfoldingPath(Triangulation TR, ParticlePositions PartPos, int MSTSteps) :_tri(TR), _partpos(PartPos), _mststeps(MSTSteps) {
			t_size = _tri.triangles.size();
			buildAdjacency();
			srand(time(NULL));
			buildEdges();
			minmetric = 99999.0;
			randlimit = 8500;
			runMST(_mststeps);
		}

		void buildAdjacency();					// Function to build adjacancy list of large triangles where each large triangle is a node of the graph
		void buildEdges();						// Populates the edges class with each edge in the large triangles graph, and assigns a random weight for that edge
		void printAdjacency();					// Self explanatory
		void sortEdges();						// Helper function to sort edges class named "ed"
		void sortEdgesfnl();					// Helper function to sort edges class named "edfnl"
		void runMST(int steps);					// Runs MST algorithm on given list of edges class
		vector<int> edgestoTris(vector<pair<int, int>>& ME);		// From the unfolding path given by MST algorithm, this function builds a list of triangles in the order which they will be unfolded
		Real buildfinalTriangles(vector<int>&OE, Real metric3);		// Prints the position of vertices based on the given unfolding path, subsequently calculates the new weight based on the unfolded particle positions
	};
}
#endif 