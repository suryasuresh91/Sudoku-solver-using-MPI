CFLAGS= -fopenmp

sudoku-mpi:
	mpicc -fopenmp -o sudoku-mpi list.c sudoku-mpi.c
	mpirun -np 4 sudoku-mpi input04.txt

sudoku-serial:
	gcc -o sudoku-serial sudoku-serial.c list.c
	./sudoku-serial input04.txt
clean:
	rm -f *.o *.~ sudoku *.gch
