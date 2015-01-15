all:
	gfortran thermo.f90 sandiego.f90 euler2d_3d.f90 weno5.f90 weno5s.f90 hllc.f90 rusanov_gas.f90 cons2prim.f90 castro.f90 drag_coeff.f90  -o sd.exe
clean:
	rm thermo.mod sd.exe *~
