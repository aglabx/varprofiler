# VarProfiler Makefile

CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -pthread -Wall -Wextra
TARGET = kmer_profiler
SOURCE = kmer_profiler.cpp

# Main target
all: $(TARGET)

# Compile k-mer profiler
$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCE)

# Clean build artifacts
clean:
	rm -f $(TARGET)
	rm -rf __pycache__
	rm -f *.pyc

# Install to /usr/local/bin
install: $(TARGET)
	cp $(TARGET) /usr/local/bin/

# Uninstall from /usr/local/bin
uninstall:
	rm -f /usr/local/bin/$(TARGET)

# Install Python dependencies
install-deps:
	pip install pandas matplotlib numpy

# Run tests (if we have test data)
test: $(TARGET)
	@echo "Running test pipeline..."
	@if [ -f "test_data/test_genome.fa" ]; then \
		python varprofiler_pipeline.py test_data/test_genome.fa -o test_output --profile quick; \
	else \
		echo "No test data found. Please add test_data/test_genome.fa"; \
	fi

# Run quick analysis
quick-run: $(TARGET)
	@if [ -z "$(GENOME)" ]; then \
		echo "Usage: make quick-run GENOME=path/to/genome.fa"; \
		exit 1; \
	fi
	python varprofiler_pipeline.py $(GENOME) -o quick_results --profile quick

# Run standard analysis  
standard-run: $(TARGET)
	@if [ -z "$(GENOME)" ]; then \
		echo "Usage: make standard-run GENOME=path/to/genome.fa"; \
		exit 1; \
	fi
	python varprofiler_pipeline.py $(GENOME) -o standard_results --profile standard

# Run detailed analysis
detailed-run: $(TARGET)
	@if [ -z "$(GENOME)" ]; then \
		echo "Usage: make detailed-run GENOME=path/to/genome.fa"; \
		exit 1; \
	fi
	python varprofiler_pipeline.py $(GENOME) -o detailed_results --profile detailed

# Run with config file
config-run: $(TARGET)
	@if [ -z "$(CONFIG)" ]; then \
		echo "Usage: make config-run CONFIG=path/to/config.json"; \
		exit 1; \
	fi
	@if [ -z "$(GENOME)" ]; then \
		echo "Usage: make config-run GENOME=path/to/genome.fa CONFIG=path/to/config.json"; \
		exit 1; \
	fi
	python varprofiler_pipeline.py $(GENOME) -c $(CONFIG)

# Help message
help:
	@echo "VarProfiler Makefile"
	@echo "===================="
	@echo ""
	@echo "Available targets:"
	@echo "  make              - Build kmer_profiler"
	@echo "  make clean        - Remove build artifacts"
	@echo "  make install      - Install to /usr/local/bin"
	@echo "  make uninstall    - Remove from /usr/local/bin"
	@echo "  make install-deps - Install Python dependencies"
	@echo "  make test         - Run test pipeline (requires test_data/test_genome.fa)"
	@echo ""
	@echo "Analysis targets (require GENOME=path/to/genome.fa):"
	@echo "  make quick-run GENOME=file.fa    - Run quick analysis (k=15, large windows)"
	@echo "  make standard-run GENOME=file.fa - Run standard analysis (k=23, default)"
	@echo "  make detailed-run GENOME=file.fa - Run detailed analysis (k=31, small windows)"
	@echo ""
	@echo "Config-based run:"
	@echo "  make config-run GENOME=file.fa CONFIG=config.json - Run with configuration file"
	@echo ""
	@echo "Example:"
	@echo "  make standard-run GENOME=human_genome.fa"

.PHONY: all clean install uninstall install-deps test quick-run standard-run detailed-run config-run help