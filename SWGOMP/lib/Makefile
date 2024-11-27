C_SRCS  := $(wildcard *.c)
S_SRCS  := $(wildcard *.S)
F_SRCS  := $(wildcard *.F90)
MODULES := omnicopy_mod.mod
IMODS   := $(patsubst %, ../include/fortran/%, $(MODULES))
OBJS    := $(patsubst %.c, %.c.o, $(C_SRCS)) $(patsubst %.S, %.o, $(S_SRCS)) $(patsubst %.F90, %.F90.o, $(F_SRCS)) 
DEPS    := $(patsubst %.c, %.dep, $(C_SRCS))
CPPFLAGS:= -I../include/cross
CFLAGS  := -g -O3 -msimd -faddress_align=64 $(CPPFLAGS)
FFLAGS  := $(CFLAGS) -cpp
all: libswgomp.a $(IMODS)
libswgomp.a: $(OBJS)
	ar cr $@ $^
%.c.o: %.c
	sw9gcc -mhost $(CFLAGS) $(CFLAGS_HOST) $< -c -o $<.host.o
	sw9gcc -mslave $(CFLAGS) $(CFLAGS_SLAVE) $< -c -o $<.slave.o
	swld -r $<.host.o $<.slave.o -o $@
%.F90.o: %.F90
	sw9gfortran -mhost $(FFLAGS) $(CFLAGS_HOST) $< -c -o $<.host.o
	sw9gfortran -mslave $(FFLAGS) $(CFLAGS_SLAVE) $< -c -o $<.slave.o
	swld -r $<.host.o $<.slave.o -o $@
../include/fortran/%.mod: %.mod
	cp $< -f $@
%-slave.o: %-slave.S
	sw9gcc -mslave $(CFLAGS) $(CFLAGS_SLAVE) $< -c -o $@
%.dep: %.c
	sw9gcc -MM -D__sw_host__ -D__sw_slave__ $(CPPFLAGS) -MT $<.o $< -o $@
-include $(DEPS)
.PHONY: clean
clean:
	rm *.o *.dep *.a
.SUFFIXES:
