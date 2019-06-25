FC = mpiifort
CFLAGS +=-O3
EXEC = a.exe

object = Element.o Particle.o Simulation.o Material.o FFI.o TimeFunction.o MPIVariable.o DataOut.o DataIn.o Inverse.o FEMCalc-gnu.f90 Preprocessing.o main.o libmetis.a

$(EXEC):$(object)
	$(FC)  $(CFLAGS) -o $(EXEC) $(object) 
$(filter %.o,$(object)): %.o: %.f90
	$(FC) $(CFLAGS) -c $< -o $@



clean:
	rm -vf *.o *.mod *.exe *.rst

clean1:
	rm -vf *.out *.err
