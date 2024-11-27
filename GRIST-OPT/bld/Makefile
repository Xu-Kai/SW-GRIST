ARGS  :=$(subst -, ,$(MAKECMDGOALS))
MODEL ?=$(word 1, $(ARGS))
MACH  ?=$(word 2, $(ARGS))
ifeq ($(MODEL), clean)
CLEAN :=Y
MODEL :=$(word 2, $(ARGS))
MACH  :=$(word 3, $(ARGS))
endif
EXE_ATM := ParGRIST-$(MODEL)-$(MACH).exe
EXE_PAR := partition-$(MODEL)-$(MACH).exe
BLD     := bld-$(MODEL)-$(MACH)

.PHONY:.DEFAULT
.DEFAULT:
	+if [ "$(CLEAN)" != "Y" ]; then \
		mkdir -p $(BLD); \
		cp Makefile.sub $(BLD)/Makefile && \
		cp path/Filepath.$(MODEL) $(BLD)/Filepath && \
		cp Macros.$(MACH) $(BLD)/Macros && \
		cd $(BLD) && ../mkSrcfiles && ../mkDepends Filepath Srcfiles > Depends && make MODEL=$(MODEL) MACHINE=$(MACH) $(MFLAGS) -r;\
	else \
		rm -rf bld-* ParGRIST-* partition-*; \
	fi
