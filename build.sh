seq="-DNOPARALLEL"

cd bitrank
gcc -c basic.c bitrankw32int.c
cd ..

gcc -O2 -c bit_array.c

echo "Compiling parallel algorithm"
gcc -std=gnu99 -o sg_par main.c defs.c planar_graph.c spanning_tree.c util.c\
    succinct_tree.c lookup_tables.c bit_array.o bitrank/basic.o\
    bitrank/bitrankw32int.o -fcilkplus -lcilkrts -lrt -lm 

echo "Compiling sequential algorithm"
gcc -std=gnu99 $seq -o sg_seq main.c defs.c planar_graph.c spanning_tree.c util.c\
    succinct_tree.c lookup_tables.c bit_array.o bitrank/basic.o\
    bitrank/bitrankw32int.o -fcilkplus -lcilkrts -lrt -lm 

echo "Compiling malloc_count version"
gcc -std=gnu99 -DMALLOC_COUNT -DNOPARALLEL -o sg_mem main.c defs.c planar_graph.c\
    spanning_tree.c util.c succinct_tree.c lookup_tables.c malloc_count.c\
    bit_array.o bitrank/basic.o bitrank/bitrankw32int.o -lrt -lm -ldl
