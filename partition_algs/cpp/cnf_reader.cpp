#include "cnf_reader.h"
#include <fstream>
#include <iostream>
#include <sstream>

Graph fromCNF(string filename) {
  FILE *file = fopen(filename.c_str(), "r");
  if (!file) {
    cerr << "Error: Could not open file " << filename << endl;
    return Graph(0); // Return an empty graph
  }

  int num_vars = 0;
  int num_clauses = 0;
  char buffer[1024];

  // Parse header
  while (fgets(buffer, sizeof(buffer), file)) {
    if (buffer[0] == 'c') {
      continue;
    } else if (buffer[0] == 'p') {
      if (sscanf(buffer, "p cnf %d %d", &num_vars, &num_clauses) != 2) {
        cerr << "Error: Invalid DIMACS header format." << endl;
        fclose(file);
        return Graph(0);
      }
      break;
    }
  }

  if (num_vars == 0) {
    cerr << "Error: Could not determine number of variables from CNF file or "
            "file is empty/malformed."
         << endl;
    fclose(file);
    return Graph(0);
  }

  Graph G(num_vars);

  // Process clauses
  int literal;
  vector<int> current_clause_literals;
  current_clause_literals.reserve(10); // Reserve some space

  while (fscanf(file, "%d", &literal) == 1) {
    if (literal == 0) {
      // End of clause
      // Check if it's a negative 2-CNF clause
      if (current_clause_literals.size() == 2 &&
          current_clause_literals[0] < 0 && current_clause_literals[1] < 0) {

        int u_dimacs = abs(current_clause_literals[0]);
        int v_dimacs = abs(current_clause_literals[1]);

        // Convert to 0-indexed graph vertices
        // DIMACS variables are 1-indexed, graph vertices are 0-indexed.
        int u_graph = u_dimacs - 1;
        int v_graph = v_dimacs - 1;

        // Validate vertex indices
        if (u_graph < 0 || u_graph >= num_vars || v_graph < 0 ||
            v_graph >= num_vars) {
          // Warning suppressed for performance? Or keep it?
          // cerr << "Warning: Literal " ...
        } else if (u_graph == v_graph) {
          // Skip self-loops
        } else {
          // Ensure u < v for addEdge, as per Graph::addEdge assertion
          G.addEdge(std::min(u_graph, v_graph), std::max(u_graph, v_graph));
        }
      }
      current_clause_literals.clear();
    } else {
      current_clause_literals.push_back(literal);
    }
  }

  fclose(file);
  return G;
}
