all: interpreter.o

compare: lama interpreter.o byterun performance
	lamac -b ../performance/Sort.lama
	`which time` -f "mySort\t%U" ./interpreter Sort.bc

interpreter.o: byterun runtime interpreter.cpp
	g++ -std=c++17 -m32 interpreter.cpp ../runtime/runtime.a ../byterun/bytelib.o -o interpreter

lama:
	$(MAKE) -C ..

byterun: runtime
	$(MAKE) -C ../byterun

runtime:
	$(MAKE) -C ../runtime

performance:
	$(MAKE) -C ../performance

clean:
	$(RM) *.a *.o *.bc *~ interpreter
