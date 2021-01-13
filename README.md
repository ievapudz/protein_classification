# Structural Alignment Program
Program written to align pairs of protein chains using the values of dihedral angles (phi and psi). Currently it is being developed to perform protein classification to SCOP classes based on identity scores.

By default during each alignment process there are seven iterations made due to usage of different block lengths ( 1-7 ) for distance calculations between dihedral angles.

* Compile with `make all`;
* Run `./structural_alignment protein_chains.txt`

## Requirements for protein chains file

### For pairwise global protein chain alignment

List down all proteins and their chains of interest (chain can be labelled with a symbol of underscore: "_" ) to align with each other.
` 1DLW_A \n
  1UX8_A `

### For classification

List down representatives of each SCOP class and add protein chain of interest at the end of the list:
` all alpha \n 
 all beta \n
 alph / beta \n
 alpha + beta \n
 small protein \n
 POI `

## Results

### Alignment

The program will align all possible protein pairs in `protein_chains.txt` file and will place produced results in `./alignment_results_numeral` directory.

### Classification

Identity scores tables can be found in `./identity_scores` directory. The predicted SCOP class for each alignment will be printed to terminal window.



