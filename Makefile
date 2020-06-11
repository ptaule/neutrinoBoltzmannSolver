EXE=main.prog

SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

CFLAGS += -Wall -Wextra -Wpedantic # -O3
CPPFLAGS += -I/scratch/gsl-2.5/ $(OPTIONS)

LDFLAGS_GSL  = -L/scratch/gsl-2.5/.libs/ -L/scratch/gsl-2.5/cblas/.libs
LDLIBS_GSL   = -lgsl -lgslcblas

LDFLAGS += $(LDFLAGS_GSL)
LDLIBS += $(LDLIBS_GSL) -lm

all: CFLAGS += -O3
debug: CFLAGS += -O0 -g

.PHONY: all, run, debug, force

all: $(EXE)
debug: $(EXE)

run: all
	./$(EXE)

$(EXE): main.o $(OBJ)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

%.o : %.cpp $(HEADERS) compiler_flags
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

compiler_flags: force
	echo '$(CPPFLAGS) $(CFLAGS)' | cmp -s - $@ || echo '$(CPPFLAGS) $(CFLAGS)' > $@

clean:
	$(RM) $(OBJ) main.o
