C++ compile: g++ -o statusToMatrix_simple -I/data/yyin/data/software/anaconda3/include/ statusToMatrix_simple.cpp

Update 12/21/2020:

Install: 
/data/yyin/data/software/anaconda3/bin/x86_64-conda_cos6-linux-gnu-g++ -o statusToMatrix_p2p -I/data/yyin/data/software/anaconda3/include/ statusToMatrix_p2p.cpp

Parameters:

-i: input file, contains all absolute path of cpg files.

-o: output file, contains a distance matrix.

-d: optional. If open(-d 1), program will print cell 2 cell information as this format:
Cellpair: region_id dist haplo_a->id haplo_a->mcode haplo_b->id haplo_b->mcode chr start_index ;

Test Case:

Only print the distance matrix:
./statusToMatrix_p2p -i cpptest/cppnames -o testtest

Print cell2cell information to stdout (cell2cell information will be very large):
./statusToMatrix_p2p -i cpptest/cppnames -o testtest -d 1 >p2plog


This command can retrieve the cell2cell information:
cat p2plog | grep Cellpair
