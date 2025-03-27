FFLAGS =
FLIBS  =
FC     = f95

OBJS   = cross_tpcf.o external.o external90.o recipes.o rand_sky.o int2str.o

crossdr72: $(OBJS) crossdr72.o
	$(FC) -o crossdr72 $(FFLAGS) $(OBJS) crossdr72.o $(FLIBS)

cross_tpcf.o: cross_tpcf.f90
	$(FC) $(FFLAGS) -c cross_tpcf.f90

external.o: external.f
	$(FC) $(FFLAGS) -c external.f

external90.o: external90.f90
	$(FC) $(FFLAGS) -c external90.f90

clean:
	rm -f *.o

cleanest:
	rm -f *.o crossdr72

