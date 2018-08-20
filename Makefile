CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix ,$(notdir $(CPP_FILES:.cpp=.o)))
CC = g++-8
CC_FLAGS := -fopenmp -std=c++11

main: $(OBJ_FILES)
	$(CC) $(CC_FLAGS) -o $@ $^

%.o: src/%.cpp src/*.h
	$(CC) $(CC_FLAGS) -c -o $@ $<
