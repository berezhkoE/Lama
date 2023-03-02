all: byterun.o bytelib.o
	$(CC) -m32 -g -o byterun byterun.o bytelib.o ../runtime/runtime.a

byterun.o: runtime.o byterun.c bytelib.c bytelib.h
	$(CC) -g -fstack-protector-all -m32 -c byterun.c bytelib.c

runtime.o:
	$(MAKE) -C ../runtime

clean:
	$(RM) *.a *.o *~ byterun
