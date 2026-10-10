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
		cp Filepath.$(MODEL) $(BLD)/Filepath && \
		cd $(BLD) && ../mkSrcfiles && ../mkDepends Filepath Srcfiles > Depends && make -j 48 MODEL=$(MODEL) MACHINE=$(MACH) $(MFLAGS) -r;\
	else \
		rm -rf bld-* ParGRIST-* partition-*; \
	fi
