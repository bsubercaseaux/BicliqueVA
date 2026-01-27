CXX = g++
CXXFLAGS = -O3 -Wall -std=c++17 -pthread

all: biclique_part transform_partition bp

biclique_part: biclique_part.o cnf_reader.o graph.o
	$(CXX) $(CXXFLAGS) -o biclique_part biclique_part.o cnf_reader.o graph.o

transform_partition: transform_partition.o
	$(CXX) $(CXXFLAGS) -o transform_partition transform_partition.o

bp: bp.o cnf_reader.o graph.o transform_partition_lib.o
	$(CXX) $(CXXFLAGS) -o bp bp.o cnf_reader.o graph.o transform_partition_lib.o

biclique_part.o: biclique_part.cpp cnf_reader.h graph.h
	$(CXX) $(CXXFLAGS) -c biclique_part.cpp

cnf_reader.o: cnf_reader.cpp cnf_reader.h graph.h
	$(CXX) $(CXXFLAGS) -c cnf_reader.cpp

graph.o: graph.cpp graph.h
	$(CXX) $(CXXFLAGS) -c graph.cpp

transform_partition.o: transform_partition.cpp transform_partition.h graph.h
	$(CXX) $(CXXFLAGS) -c transform_partition.cpp

transform_partition_lib.o: transform_partition.cpp transform_partition.h graph.h
	$(CXX) $(CXXFLAGS) -DBICLIQUE_PARTITIONS_NO_MAIN -c transform_partition.cpp -o transform_partition_lib.o

bp.o: bp.cpp cnf_reader.h graph.h transform_partition.h
	$(CXX) $(CXXFLAGS) -c bp.cpp

clean:
	rm -f *.o biclique_part transform_partition bp
