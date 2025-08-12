CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -pthread
TARGET = kmer_profiler
SOURCE = kmer_profiler.cpp

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCE)

clean:
	rm -f $(TARGET)

install: $(TARGET)
	cp $(TARGET) /usr/local/bin/

uninstall:
	rm -f /usr/local/bin/$(TARGET)

.PHONY: all clean install uninstall