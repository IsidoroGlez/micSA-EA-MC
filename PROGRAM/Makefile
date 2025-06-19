# Global defines
SHELL=bash
L := 160
NBETAS := 1
NBITSPREBUSQUEDAS := 10
GLOBAL_DEFINES := -DL=$(L) -DNBETAS=$(NBETAS) -DNUMBITSPREBUSQUEDAS=$(NBITSPREBUSQUEDAS) -DJANUS -DWALK

# Directory paths
SRC_DIR := src
INCLUDE_DIR := include
OBJ_DIR := obj
BIN_DIR := bin

# Compiler flags
C_FLAGS := -c -O $(GLOBAL_DEFINES) -m64 -mcmodel=large
CUDADIR := /usr/local/cuda
CUDAINCLUDEDIR := $(CUDADIR)/include
NVCC_FLAGS := -c -O3 --ptxas-options=-v -DHOSTREDU -arch=sm_90 $(GLOBAL_DEFINES) -I$(CUDAINCLUDEDIR) -L$(CUDADIR)/lib -lcudart -Xcompiler -mcmodel=large
LDFLAGS := -lm -lquadmath

# Compiler commands
CC := /usr/bin/g++
NVCC := nvcc

# Target executable name
TARGET := $(BIN_DIR)/EA_CUBE_L${L}_NB${NBETAS}_NBPRE${NBITSPREBUSQUEDAS}

# Source files
C_SOURCES := $(wildcard $(SRC_DIR)/*.c)
CU_SOURCES := $(wildcard $(SRC_DIR)/*.cu)
HEADERS := $(wildcard $(INCLUDE_DIR)/*.h)

# Object files
C_OBJECTS := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(C_SOURCES))
CU_OBJECTS := $(patsubst $(SRC_DIR)/%.cu, $(OBJ_DIR)/%.o, $(CU_SOURCES))

# All object files
OBJECTS := $(C_OBJECTS) $(CU_OBJECTS)

# Default target
all: $(TARGET)

# Link the final executable
$(TARGET): $(OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(NVCC) $^ -o $@ $(LDFLAGS) -lcudart -Xcompiler -mcmodel=large

# Compile C source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS) Makefile
	@mkdir -p $(OBJ_DIR)
	$(CC) $(C_FLAGS) -c $< -o $@ $(LDFLAGS)

# Compile CUDA source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu $(HEADERS) Makefile
	@mkdir -p $(OBJ_DIR)
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@ $(LDFLAGS)

# Clean up
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)/*~ $(SRC_DIR)/*~ $(INCLUDE_DIR)/*~ *~

# Phony targets
.PHONY: all clean
