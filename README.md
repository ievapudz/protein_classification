
# Dataset Generator
### For protein classification into SCOP classes using recurrent neural network

1. Build with `make`;
2. Run with ` ./dataset_generator [training | training | testing].txt`;

Output will be a .csv file with m * n table, where m is equal to number of files and n is equal to 
2001 (phi and psi protein angle pairs, in which each angle is represented as a pair of its sin and
cos values, and the label of the SCOP class). 

Note: file that is given as a command line argument should contain a list of labelled protein 
chains.
All protein chains should be found in ./mmCIF_files directory.

Improvements that need to be done:
1. Command line argument that determines how many elements of each class our dataset will 
contain.
