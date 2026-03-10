#include "transform_partition.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

// Helper to read partition
vector<Biclique> readPartition(const string &filename) {
  vector<Biclique> partition;
  ifstream file(filename);
  if (!file.is_open()) {
    cerr << "Error: Could not open partition file " << filename << endl;
    return partition;
  }

  int num_bicliques;
  file >> num_bicliques;
  for (int i = 0; i < num_bicliques; ++i) {
    int l_size, r_size;
    file >> l_size >> r_size;
    Biclique b;
    b.left_vertices.resize(l_size);
    b.right_vertices.resize(r_size);
    for (int j = 0; j < l_size; ++j)
      file >> b.left_vertices[j];
    for (int j = 0; j < r_size; ++j)
      file >> b.right_vertices[j];
    partition.push_back(b);
  }
  return partition;
}

void transformCNF(const string &original_cnf, const vector<Biclique> &partition,
                  const string &output_cnf) {
  // 1. Read original CNF to get num_vars and non-biclique clauses
  ifstream cnf_in(original_cnf);
  if (!cnf_in.is_open()) {
    cerr << "Error: Could not open original CNF " << original_cnf << endl;
    return;
  }

  int num_vars = 0;
  int num_clauses = 0;
  // Store kept clauses in a flat vector to avoid overhead of
  // vector<vector<int>> Format: lit1 lit2 ... 0 lit1 lit2 ... 0
  vector<int> other_clauses_flat;
  int other_clauses_count = 0;
  string line;

  while (getline(cnf_in, line)) {
    if (line.empty() || line[0] == 'c')
      continue;
    if (line[0] == 'p') {
      stringstream ss(line);
      string p, cnf;
      ss >> p >> cnf >> num_vars >> num_clauses;
      continue;
    }
    stringstream ss(line);
    int lit;
    vector<int> clause;
    while (ss >> lit && lit != 0) {
      clause.push_back(lit);
    }
    // Logic from bicliques_to_formula.py:
    // if len(clause) != 2: keep
    // elif clause[0] > 0 or clause[1] > 0: keep
    // (negative 2-clauses are the edges, which we replace)
    if (clause.size() != 2) {
      for (int l : clause)
        other_clauses_flat.push_back(l);
      other_clauses_flat.push_back(0);
      other_clauses_count++;
    } else if (clause[0] > 0 || clause[1] > 0) {
      for (int l : clause)
        other_clauses_flat.push_back(l);
      other_clauses_flat.push_back(0);
      other_clauses_count++;
    }
  }
  cnf_in.close();

  // 3. Generate new clauses
  // We need new variables for bicliques if we use the auxiliary variable
  // encoding Logic from bicliques_to_formula.py: if len(L)*len(R) - len(L) -
  // len(R) < 1:
  //   use direct encoding: -u -v for u in L, v in R
  // else:
  //   use auxiliary variable y_b: -u y_b, -v -y_b

  // We need to assign IDs to y_b. They should start after num_vars.
  int next_var = num_vars + 1;
  long long new_clauses_count = 0;

  // First pass: count new clauses and assign variables
  for (const auto &b : partition) {
    long long l_sz = b.left_vertices.size();
    long long r_sz = b.right_vertices.size();

    if (l_sz * r_sz - l_sz - r_sz < 1) {
      new_clauses_count += l_sz * r_sz;
    } else {
      next_var++; // Allocate y_b
      new_clauses_count += l_sz + r_sz;
    }
  }

  // 4. Write output CNF
  // Use C-style IO for performance
  FILE *cnf_out = fopen(output_cnf.c_str(), "w");
  if (!cnf_out) {
    cerr << "Error: Could not open output CNF " << output_cnf << endl;
    return;
  }

  // Set buffer size to 64KB
  char buffer[65536];
  setvbuf(cnf_out, buffer, _IOFBF, 65536);

  fprintf(cnf_out, "p cnf %d %lld\n", (next_var - 1),
          (other_clauses_count + new_clauses_count));

  // Write kept clauses
  for (size_t i = 0; i < other_clauses_flat.size(); ++i) {
    fprintf(cnf_out, "%d", other_clauses_flat[i]);
    if (other_clauses_flat[i] == 0)
      fprintf(cnf_out, "\n");
    else
      fprintf(cnf_out, " ");
  }

  // Generate and write new clauses on the fly
  int current_y_b_var = num_vars + 1;
  for (const auto &b : partition) {
    long long l_sz = b.left_vertices.size();
    long long r_sz = b.right_vertices.size();

    if (l_sz * r_sz - l_sz - r_sz < 1) {
      for (int u : b.left_vertices) {
        for (int v : b.right_vertices) {
          fprintf(cnf_out, "%d %d 0\n", -(u + 1), -(v + 1));
        }
      }
    } else {
      int y_b = current_y_b_var++;
      for (int u : b.left_vertices) {
        fprintf(cnf_out, "%d %d 0\n", -(u + 1), y_b);
      }
      for (int v : b.right_vertices) {
        fprintf(cnf_out, "%d %d 0\n", -(v + 1), -y_b);
      }
    }
  }
  fclose(cnf_out);

  // cout << "Transformed CNF written to " << output_cnf << endl;
  // cout << "Original vars: " << num_vars << ", New vars: " << (next_var - 1)
  //      << endl;
  // cout << "Original clauses (kept): " << other_clauses_count
  //      << ", New clauses: " << new_clauses_count << endl;
}

#ifndef BICLIQUE_PARTITIONS_NO_MAIN
int main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << "Usage: " << argv[0]
         << " <original_cnf> <partition_file> <output_cnf>" << endl;
    return 1;
  }

  string original_cnf = argv[1];
  string partition_file = argv[2];
  string output_cnf = argv[3];

  vector<Biclique> partition = readPartition(partition_file);
  transformCNF(original_cnf, partition, output_cnf);

  return 0;
}
#endif
