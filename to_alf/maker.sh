gcc -I../util -Wall -O3   -c -o test_alf.o test_alf.c
gcc -I../util -Wall -O3   -c -o t_o.o t_o.c
gcc -I../util -Wall -O3   -c -o doubly_linked_list.o ../util/doubly_linked_list.c
gcc test_alf.o t_o.o doubly_linked_list.o ../util/epsilon.o ../util/memfunc.o   -o test_alf -lm -lX11
