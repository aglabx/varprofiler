CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -march=native -pthread

TARGETS = kmer_profiler cenpb_finder motif_discovery

all: $(TARGETS)

kmer_profiler: kmer_profiler.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

cenpb_finder: cenpb_finder.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

motif_discovery: motif_discovery.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

install: $(TARGETS)
	mkdir -p ../../bin
	cp $(TARGETS) ../../bin/

clean:
	rm -f $(TARGETS)

.PHONY: all install clean