EXE=main.prog

SRC_DIR = src
OBJ_DIR = obj
INC_DIR = include

SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
HEADERS = $(wildcard $(INC_DIR)/*.hpp)

CXX ?= g++

CXXFLAGS += -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion \
			-Wcast-align -Wunused -Wlogical-op -Wnull-dereference \
			-std=c++17

CPPFLAGS += -DHAVE_INLINE -I/space/ge52sir/local/include/

all:   CXXFLAGS += -DDEBUG=0 -O3
debug: CPPFLAGS += -DDEBUG=1
debug: CXXFLAGS += -O0 -g

LDFLAGS += -L/space/ge52sir/local/lib/
LDLIBS += -lgsl -lgslcblas

.PHONY: all run debug force

all: $(EXE)
debug: $(EXE)

run: all
	./$(EXE)

$(EXE): main.o $(OBJ)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

main.o: main.cpp $(HEADERS) compiler_flags
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS) compiler_flags
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

compiler_flags: force
	echo '$(CPPFLAGS) $(CXXFLAGS)' | cmp -s - $@ || echo '$(CPPFLAGS) $(CXXFLAGS)' > $@

clean:
	$(RM) $(OBJ) main.o
