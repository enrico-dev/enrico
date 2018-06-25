CXX = g++
FC =gfortran
CXXFLAGS = -g -I.
FCFLAGS = -g -I.

driver: driver.o fobj.o
	$(CXX) -o $@ $^

driver.o: driver.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^

fobj.o: fobj.f90
	$(FC) $(FCFLAGS) -c -o $@ $^

clean:
	rm -f driver driver.o fobj.o fobj.mod

.PHONY: clean

