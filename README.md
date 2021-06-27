
# Dataset Generator
### For protein classification into SCOP classes using recurrent neural network

1. Build with `make`;
2. Run with ` ./dataset_generator [training | training | testing].txt`;

Output will be a .csv file with m x n table, where  
m = number of files  
n = 2001 (phi and psi protein angle pairs, in which each angle is represented as a pair of its sin and
cos values plus the label of the SCOP class). 

Note: file that is given as a command line argument should contain a list of labelled protein 
chains.
All protein chains should be found in ./mmCIF_files directory.

Improvements that need to be done:
1. Command line argument that determines how many elements of each class our dataset will 
contain.
2. Download one SCOP class at a time and calculate phi and psi anglesand save to separate database files
    * all_alpha
    * all_beta
    * alpha_plus_beta
    * alpha_slash_beta
    * small_protein
3. Parse separate database files for each class (60% for training, 20% for validation, 20% for testing)
