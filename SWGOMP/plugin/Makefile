SRCS     := $(wildcard *.cpp)
OBJS     := $(patsubst %.cpp, %.o, $(SRCS))
DEPS     := $(patsubst %.cpp, %.dep, $(SRCS))
CXXFLAGS := -fno-rtti -O3 -fPIC -std=gnu++11 -g
INCLUDES := -I$(shell sw9gfortran -print-file-name=plugin)/include -I$(shell dirname $(shell which sw9gfortran))/../lib_for_gcc/include/ -I../include/cross
all: swgomp.so
swgomp.so: $(OBJS)
	g++ $^ -shared -g -o $@
%.o: %.cpp
	g++ $< $(INCLUDES) $(CXXFLAGS) -c -o $@
%.dep: %.cpp
	g++ -MM -D__sw_host__ -D__sw_slave__ $(INCLUDES) -MT $<.o $< -o $@
-include $(DEPS)
clean:
	-rm swgomp.so *.o *.dep -f
