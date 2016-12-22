CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(patsubst src/%.cpp,bin/%.o,$(CPP_FILES))
CPP_FILES += lib/libmrio/src/MRIOIndexSet.cpp lib/libmrio/src/MRIOTable.cpp
OBJ_FILES += bin/MRIOIndexSet.o bin/MRIOTable.o
.SECONDARY: OBJ_FILES
LD_FLAGS := -Wl,--gc-sections -lnetcdf_c++4
CC_FLAGS := -std=c++0x -I lib/cpp-library -I lib/libmrio/include -fdata-sections -ffunction-sections -Wshadow

all: fast

fast: CC_FLAGS += -fopenmp -O3
fast: LD_FLAGS += -fopenmp
fast: gcc

debug: CC_FLAGS += -g -DDEBUG
debug: gcc

gcc: CXX = g++
gcc: makedir
gcc: binaries

clean:
	@rm -f $(OBJ_FILES) bin/libmrio psi

makedir:
	@mkdir -p bin

binaries: psi

psi: bin/psi.o $(OBJ_FILES)
	@echo Linking to $@
	@$(CXX) -o $@ $^ $(LD_FLAGS)

bin/psi.o: src/psi.cpp lib/cpp-library/nvector.h
	@mkdir -p bin
	@echo Compiling $<
	@$(CXX) $(CC_FLAGS) $(ACC_OPTIONS) -c -o $@ $<

bin/MRIO%.o: lib/libmrio/src/MRIO%.cpp lib/libmrio/include/MRIO%.h
	@echo Compiling $<
	@$(CXX) $(CC_FLAGS) $(ACC_OPTIONS) -c -o $@ $<

bin/%.o: src/%.cpp include/%.h
	@echo Compiling $<
	@$(CXX) $(CC_FLAGS) $(ACC_OPTIONS) -c -o $@ $<
