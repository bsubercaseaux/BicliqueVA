#include "graph.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <thread>
#include <tuple>
#include <vector>

using namespace std;

// ---------- utilities for Density Aware Algorithm ----------

static double h2(double gamma) {
  if (gamma <= 0 || gamma >= 1)
    return 0.0;
  return -(gamma * log2(gamma) + (1 - gamma) * log2(1 - gamma));
}

static bool choose_leq_threshold(int L, int y, long long T) {
  if (y < 0 || y > L)
    return y == 0;
  if (y == 0 || y == L)
    return 1 <= T;
  if (y == 1 || y == L - 1)
    return L <= T;

  y = min(y, L - y);
  long long num = 1;
  for (int i = 1; i <= y; ++i) {
    num = (num * (L - y + i)) / i;
    if (num > T)
      return false;
  }
  return num <= T;
}

static vector<int> precompute_max_len_per_y(int r, long long T) {
  vector<int> X(r + 1);
  for (int y = 1; y <= r; ++y) {
    int lo = y, hi = r, ans = y - 1;
    while (lo <= hi) {
      int mid = (lo + hi) / 2;
      if (choose_leq_threshold(mid, y, T - 1)) {
        ans = mid;
        lo = mid + 1;
      } else {
        hi = mid - 1;
      }
    }
    X[y] = ans;
  }
  return X;
}

static int circ_dist(int i, int j, int t) {
  int d = abs(i - j);
  return min(d, t - d);
}

static bool R_rel_density(int i, int j, int t) {
  if (i == j)
    return false;
  if (i < j)
    return circ_dist(i, j, t) % 2 == 0;
  else
    return circ_dist(i, j, t) % 2 == 1;
}

Graph::Graph(int vertices) {
  this->vertices = vertices;
  this->adjList = vector<vector<int>>(vertices, vector<int>());
  // Flattened adjacency matrix initialized to 0 (false)
  // Use uint8_t for byte-addressability and speed over vector<bool>
  this->adjMatrix.resize((size_t)vertices * vertices, 0);
}

void Graph::addEdge(int u, int v) {
  // assert that 0 <= u < v < this->vertices.
  if (u < 0 || v <= u || v >= this->vertices) {
    cerr << "Error in addEdge, invalid indices " << u << ", " << v << endl;
    return;
  }
  if (this->adjMatrix[(size_t)u * vertices + v]) {
    cerr << "edge already present for " << u << ", " << v << endl;
    return;
  }
  this->adjList[u].push_back(v);
  this->adjList[v].push_back(u);
  this->adjMatrix[(size_t)u * vertices + v] = 1;
  this->adjMatrix[(size_t)v * vertices + u] = 1;
  this->edges++;
}

bool Graph::isEdge(int u, int v) const {
  if (u < 0 || v < 0 || u >= vertices || v >= vertices)
    return false;
  return this->adjMatrix[(size_t)u * vertices + v];
}

int Graph::getNumVertices() const { return vertices; }

int Graph::getNumEdges() const { return edges; }

/* --------- helpers for Algorithm 1 --------- */

// tournament orientation R(i, j)
static inline bool R_orient(int i, int j, int r) {
  if (i == j)
    return false;
  // Use linear distance to match Python implementation
  int d = abs(i - j);
  if (i < j)
    return (d % 2 == 0);
  else
    return (d % 2 == 1);
}

vector<Biclique> Graph::getBicliquePartition(int num_threads) {
  vector<Biclique> partition;
  int n = this->vertices;
  if (n == 0)
    return partition;

  // M := floor(lg n - 1.1 lg lg n), but keep it >= 1 and handle tiny n
  // robustly.
  // double lg_n = (n > 1) ? log2(static_cast<double>(n)) : 1.0;
  // double lg_lg_n = (lg_n > 0) ? log2(lg_n) : 0.0;
  // cout << "lg_n = " << lg_n << ", lg_lg_n = " << lg_lg_n << endl;
  // cout << "lg_n - 1.1 lg_lg_n = " << lg_n - 1.1 * lg_lg_n << endl;
  // int M = static_cast<int>(floor(lg_n - 1.1 * lg_lg_n));
  // if (M < 1)
  //   M = 1;

  int M = 2;
  int best_est = 1e9;
  for (int i = 3; i <= 15; ++i) {
    int estS = n * (1 << (i - 1));
    int estA = n * (n / (2 * i));
    int est = estS + estA;
    cout << "i = " << i << ", est = " << est << endl;
    if (est < best_est) {
      best_est = est;
      M = i;
    }
  }

  cout << "Using M = " << M << " for n = " << n << endl;

  // Partition V into r parts of size at most M: P0, P1, ..., P_{r-1}
  int r = (n + M - 1) / M;
  vector<vector<int>> parts(r);
  vector<int> partOf(n, -1);
  for (int i = 0; i < n; ++i) {
    int p = i / M;
    parts[p].push_back(i);
    partOf[i] = p;
  }

  if (num_threads < 1)
    num_threads = 1;
  vector<std::thread> threads;
  vector<vector<Biclique>> thread_partitions(num_threads);

  auto worker = [&](int thread_id, int start_r, int end_r) {
    // Thread-local buckets
    size_t max_subset_count = 1ULL << M;
    vector<vector<int>> buckets(max_subset_count);
    vector<size_t> masks(n, 0); // Thread-local masks buffer

    for (int i = start_r; i < end_r; ++i) {
      const vector<int> &Pi = parts[i];
      int partSize = static_cast<int>(Pi.size());
      if (partSize == 0)
        continue;
      size_t subsetCount = 1ULL << partSize;

      if (buckets.size() < subsetCount)
        buckets.resize(subsetCount);

      // --- Step 1: Bicliques (S, A(S)) ---

      // Build A(S): scan all vertices v that point to part i in the tournament
      // R Reset masks for this part iteration Optimization: We only need to
      // reset masks that we touch? Actually, we iterate all v, so we must reset
      // all masks or just overwrite. But we accumulate bits `|=`, so we must
      // reset. `std::fill` is fast enough.
      std::fill(masks.begin(), masks.end(), 0);

      for (int bit = 0; bit < partSize; ++bit) {
        int u = Pi[bit];
        size_t bitVal = (size_t(1) << bit);
        size_t rowOffset = (size_t)u * vertices;
        for (int v = 0; v < n; ++v) {
          if (adjMatrix[rowOffset + v]) {
            masks[v] |= bitVal;
          }
        }
      }

      for (int v = 0; v < n; ++v) {
        if (masks[v] == 0)
          continue;
        int pv = partOf[v];
        if (!R_orient(pv, i, r))
          continue;
        buckets[masks[v]].push_back(v);
      }

      for (size_t mask = 1; mask < subsetCount; ++mask) {
        vector<int> &Aset = buckets[mask];
        if (Aset.empty())
          continue;
        vector<int> S;
        S.reserve(partSize);
        for (int bit = 0; bit < partSize; ++bit) {
          if (mask & (size_t(1) << bit))
            S.push_back(Pi[bit]);
        }
        Biclique b;
        b.left_vertices = std::move(S);
        b.right_vertices = std::move(Aset);
        thread_partitions[thread_id].push_back(std::move(b));
      }

      // --- Step 2: Within-part edges ---
      int sz = (int)Pi.size();
      for (int a = 0; a < sz; ++a) {
        for (int b = a + 1; b < sz; ++b) {
          int u = Pi[a], v = Pi[b];
          if (adjMatrix[(size_t)u * vertices + v]) {
            Biclique bq;
            bq.left_vertices = {u};
            bq.right_vertices = {v};
            thread_partitions[thread_id].push_back(std::move(bq));
          }
        }
      }
    }
  };

  int chunk_size = (r + num_threads - 1) / num_threads;
  for (int t = 0; t < num_threads; ++t) {
    int start = t * chunk_size;
    int end = min(start + chunk_size, r);
    if (start < end) {
      threads.emplace_back(worker, t, start, end);
    }
  }

  for (auto &t : threads) {
    if (t.joinable())
      t.join();
  }

  size_t total_bicliques = 0;
  for (const auto &tp : thread_partitions)
    total_bicliques += tp.size();
  partition.reserve(total_bicliques);

  for (auto &tp : thread_partitions) {
    partition.insert(partition.end(), std::make_move_iterator(tp.begin()),
                     std::make_move_iterator(tp.end()));
  }

  return partition;
}

vector<Biclique> Graph::getDensityAwareBicliquePartition(double gamma) {
  int n = this->vertices;
  if (n == 0)
    return {};

  if (gamma < 0) {
    long long denom = (long long)n * (n - 1) / 2;
    gamma = (denom > 0) ? (double)this->edges / denom : 0.0;
  }

  auto lg = [](double x) { return log2(max(x, 2.0)); };
  double h = max(h2(gamma), 1e-9);

  // r = ceil( lg^2 n / h2(gamma) )
  int r = max(1, (int)ceil(pow(lg(n), 2) / h));

  // cout << "Using r = " << r << " for n = " << n << ", gamma = " << gamma <<
  // endl;

  vector<vector<int>> parts;
  vector<pair<int, int>> pos_in_part(n, {-1, -1});
  vector<int> order(n);
  for (int i = 0; i < n; ++i)
    order[i] = i;

  for (int i = 0; i < n; i += r) {
    vector<int> block;
    for (int k = 0; k < r && i + k < n; ++k) {
      block.push_back(order[i + k]);
    }
    parts.push_back(block);
    int bidx = parts.size() - 1;
    for (int j = 0; j < block.size(); ++j) {
      pos_in_part[block[j]] = {bidx, j + 1}; // 1-based index
    }
  }

  int t = parts.size();
  long long T = max(6LL, (long long)(n / pow(r, 4)));
  vector<int> X = precompute_max_len_per_y(r, T);

  // neigh_pos: key (v, i) -> list of positions
  // Index = v * t + i.
  vector<vector<int>> neigh_pos(n * t);
  vector<pair<int, int>> within_part_edges;

  for (int u = 0; u < n; ++u) {
    for (int v : adjList[u]) {
      if (u < v) { // Process each edge once
        pair<int, int> pu = pos_in_part[u];
        pair<int, int> pv = pos_in_part[v];
        int iu = pu.first, ju = pu.second;
        int iv = pv.first, jv = pv.second;

        if (iu == iv) {
          within_part_edges.push_back({u, v});
        } else {
          if (R_rel_density(iu, iv, t)) {
            neigh_pos[u * t + iv].push_back(jv);
          }
          if (R_rel_density(iv, iu, t)) {
            neigh_pos[v * t + iu].push_back(ju);
          }
        }
      }
    }
  }

  for (auto &vec : neigh_pos) {
    sort(vec.begin(), vec.end());
  }

  // Ai buckets
  // Key: (i, a, b, S)
  // S is vector<int>
  map<tuple<int, int, int, vector<int>>, vector<int>> Ai;

  for (int v = 0; v < n; ++v) {
    for (int i = 0; i < t; ++i) {
      const vector<int> &P = neigh_pos[v * t + i];
      if (P.empty())
        continue;

      int mlen = parts[i].size();
      int idx = 0;
      while (idx < P.size()) {
        int a = P[idx];
        int tmax = 1;
        while (idx + tmax < P.size()) {
          int cand_t = tmax + 1;
          if (cand_t > r)
            break;
          int needed_len = P[idx + cand_t - 1] - a + 1;
          if (needed_len <= X[cand_t]) {
            tmax = cand_t;
          } else {
            break;
          }
        }

        int next_break = (idx + tmax < P.size()) ? P[idx + tmax] - 1 : mlen;
        int b_cap = a + X[tmax] - 1;
        int b = min(b_cap, next_break);

        vector<int> positions;
        for (int k = 0; k < tmax; ++k) {
          positions.push_back(P[idx + k]);
        }

        Ai[{i, a, b, positions}].push_back(v);

        idx += tmax;
      }
    }
  }

  vector<Biclique> bicliques;
  for (auto const &[key, Averts] : Ai) {
    int i = get<0>(key);
    // int a = get<1>(key);
    // int b = get<2>(key);
    const vector<int> &positions = get<3>(key);

    const vector<int> &P_i = parts[i];
    vector<int> Lverts;
    for (int p : positions) {
      // p is 1-based index in P_i
      Lverts.push_back(P_i[p - 1]);
    }

    if (!Lverts.empty() && !Averts.empty()) {
      Biclique bq;
      bq.left_vertices = Lverts;
      bq.right_vertices = Averts;
      bicliques.push_back(bq);
    }
  }

  for (auto &edge : within_part_edges) {
    Biclique bq;
    bq.left_vertices = {edge.first};
    bq.right_vertices = {edge.second};
    bicliques.push_back(bq);
  }

  return bicliques;
}

Graph getCompleteGraph(int n) {
  Graph G = Graph(n);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      G.addEdge(i, j);
    }
  }
  return G;
}

Graph getRandomGraph(int n, double p) {
  Graph G = Graph(n);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double rnd = static_cast<double>(rand()) / RAND_MAX;
      if (rnd < p) {
        G.addEdge(i, j);
      }
    }
  }
  return G;
}
