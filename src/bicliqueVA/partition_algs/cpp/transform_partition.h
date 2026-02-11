#ifndef TRANSFORM_PARTITION_H
#define TRANSFORM_PARTITION_H

#include "graph.h"
#include <string>
#include <vector>

using namespace std;

void transformCNF(const string &original_cnf, const vector<Biclique> &partition,
                  const string &output_cnf);

#endif
