#ifndef GRAPH_H
#define GRAPH_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

class Biclique {
public:
  vector<int> left_vertices;
  vector<int> right_vertices;
};

class Graph {
public:
  Graph(int vertices);
  void addEdge(int u, int v);
  vector<Biclique> getBicliquePartition(int num_threads = 1);
  vector<Biclique> getDensityAwareBicliquePartition(double gamma = -1.0);

  bool isEdge(int u, int v) const;
  int getNumVertices() const;
  int getNumEdges() const;

private:
  int vertices;
  int edges = 0;
  vector<vector<int>> adjList;
  vector<uint8_t> adjMatrix; // Flattened 1D array: row * vertices + col
};

Graph getCompleteGraph(int n);
Graph getRandomGraph(int n, double p);

#endif
