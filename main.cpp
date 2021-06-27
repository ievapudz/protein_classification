#include <fstream>
#include <algorithm>
#include <vector>
#include "./src/scoring_matrix.h"
#include "./src/distance_score_matrix.h"
#include "./src/distance_matrix.h"
#include "./src/biological_structures.h"
#include "./src/cif_parser.h"
#include "./src/file_workflow.h"
#include "./src/block_distance.h"
#include "./src/constants.hpp"
#include "./src/identity_score_table.h"
#include "./src/dataset_generator.h"
#include "preparatory_phase.h"
#include "calculation_phase.h"
#include "representation_phase.h"

/*
    Program that extracts the dihedral angle chains from proteins.
    These chains are written into
    the csv file along with the SCOP class that a certain protein belongs to.
 
    Input labelled (with SCOP classes) protein chains file.
    Output: csv file with dihedral angles
*/

int main(int argc, const char * argv[]){

    try{
        DatasetGenerator dataset_gen(argv[1]);
        dataset_gen.readFile();
        dataset_gen.extractChains();
        dataset_gen.extractLabels();
        dataset_gen.extractClassClusterIndeces();
        
        int number_of_indeces = atoi(argv[2]);
        dataset_gen.extractUsedIndeces(number_of_indeces);
        dataset_gen.extractAuthSeqIds();
        dataset_gen.determineChainsFiles();
        std::vector<int> used_indeces = dataset_gen.getUsedIndeces();
        std::cout << used_indeces.size() << std::endl;
        for(int j = 0; j < used_indeces.size(); j++){
            std::cout << j << std::endl;
            std::ifstream infile(dataset_gen.getChains()[used_indeces[j]]);
            if(infile.good()){
                CIFParser parser;
            
                parser.setFilePath(dataset_gen.getChains()[used_indeces[j]]);
                Protein protein;
                parser.parseAtomSiteColumns();
                protein.setAllAtoms(parser.parseAtoms());
                
                protein.filterAtoms("N", dataset_gen.getAuthSeqIds()[used_indeces[j]]);
                protein.filterAtoms("CA", dataset_gen.getAuthSeqIds()[used_indeces[j]]);
                protein.filterAtoms("C", dataset_gen.getAuthSeqIds()[used_indeces[j]]);
                
                std::vector<double> phi_angles = protein.getPhiAngles();
                std::vector<double> psi_angles = protein.getPsiAngles();
                
                int iterations = phi_angles.size();
                if(psi_angles.size() < iterations){
                    iterations = psi_angles.size();
                }
                
                if((phi_angles.size() < number_of_indeces) || (psi_angles.size() < number_of_indeces)){
                    int additional_phi = number_of_indeces - phi_angles.size();
                    int additional_psi = number_of_indeces - psi_angles.size();
                    for(int k = 0; k < additional_phi; k++){
                        phi_angles.push_back(0.00);
                    }
                    for(int k = 0; k < additional_psi; k++){
                        psi_angles.push_back(0.00);
                    }
                }
                
                dataset_gen.createSample(number_of_indeces, phi_angles, psi_angles);
                dataset_gen.addSample(used_indeces[j]);
            }else{
                std::cout << "No file: " << dataset_gen.getChains()[used_indeces[j]] << std::endl;
            }
        }
        dataset_gen.writeSamples();

    }catch(const std::length_error& le){
        std::cerr << le.what() << std::endl;
    }
}
