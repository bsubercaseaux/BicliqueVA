#include "cnf_reader.h"
#include "graph.h"
#include "transform_partition.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char *argv[]) {
  ios::sync_with_stdio(false);
  cin.tie(NULL);

  if (argc < 3 || argc > 4) {
    cerr << "Usage: " << argv[0] << " <input_cnf> <output_cnf> [num_threads]"
         << endl;
    return 1;
  }

  string input_cnf = argv[1];
  string output_cnf = argv[2];
  int num_threads = 1;
  if (argc == 4) {
    num_threads = atoi(argv[3]);
  }

  // 1. Read Graph from CNF
  Graph G = fromCNF(input_cnf);
  if (G.getNumVertices() == 0) {
    cerr << "Error reading CNF or empty graph." << endl;
    return 1;
  }
  cout << "Read graph with " << G.getNumVertices() << " vertices and "
       << G.getNumEdges() << " edges." << endl;

  // 2. Compute Biclique Partition
  vector<Biclique> partition = G.getBicliquePartition(num_threads);
  // vector<Biclique> partition = G.getDensityAwareBicliquePartition();
  int total_weight = 0;
  for (const auto &b : partition) {
    auto l = b.left_vertices.size();
    auto r = b.right_vertices.size();
    total_weight += min(l + r, l * r);
    if (l + r < l * r) {
      // cout << "saved" << (l * r) - (l + r) << endl;
    }
  }
  cout << "Computed partition with " << partition.size()
       << " bicliques, total weight: " << total_weight << endl;

  // 3. Transform to new CNF
  transformCNF(input_cnf, partition, output_cnf);

  return 0;
}
