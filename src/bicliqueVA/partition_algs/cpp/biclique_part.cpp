#include "cnf_reader.h"
#include "graph.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

#ifndef BICLIQUE_PARTITIONS_NO_MAIN
int main(int argc, char *argv[]) {
  if (argc < 2 || argc > 3) {
    cerr << "Usage: " << argv[0] << " <input_cnf> [output_partition_file]"
         << endl;
    return 1;
  }
  string filename = argv[1];
  Graph G = fromCNF(filename);
  // cout << "Read graph with " << G.getNumVertices() << " vertices and "
  //      << G.getNumEdges() << " edges." << endl;

  // vector<Biclique> partition = G.getDensityAwareBicliquePartition();
  vector<Biclique> partition = G.getBicliquePartition();
  int total_weight = 0;
  for (int t = 0; t < (int)partition.size(); ++t) {
    auto &b = partition[t];
    auto l = b.left_vertices.size();
    auto r = b.right_vertices.size();
    total_weight += min(l + r, l * r);
  }

  // cout << "Total bicliques: " << partition.size()
  //      << ", total weight: " << total_weight << endl;

  if (argc == 3) {
    string out_filename = argv[2];
    ofstream outfile(out_filename);
    if (!outfile.is_open()) {
      cerr << "Error: Could not open output file " << out_filename << endl;
      return 1;
    }
    outfile << partition.size() << endl;
    for (const auto &b : partition) {
      outfile << b.left_vertices.size() << " " << b.right_vertices.size()
              << endl;
      for (size_t i = 0; i < b.left_vertices.size(); ++i) {
        outfile << b.left_vertices[i]
                << (i == b.left_vertices.size() - 1 ? "" : " ");
      }
      outfile << endl;
      for (size_t i = 0; i < b.right_vertices.size(); ++i) {
        outfile << b.right_vertices[i]
                << (i == b.right_vertices.size() - 1 ? "" : " ");
      }
      outfile << endl;
    }
    // cout << "Partition written to " << out_filename << endl;
  }

  return 0;
}
#endif
