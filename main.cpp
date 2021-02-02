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
#include "preparatory_phase.h"
#include "calculation_phase.h"
#include "representation_phase.h"

int main(int argc, const char * argv[]){

    try{
        std::ifstream input;
        input.open(argv[1]);
        std::string read;
        std::vector<std::string> file_content;
        std::vector<std::string> chains;
        std::vector<std::string> auth_seq_ids;
        std::vector<std::string> labels;
        std::vector<std::string> samples;
        
        while(!input.eof()){
            input >> read;
            file_content.push_back(read);
        }
        file_content.pop_back();
        
        for(int i = 0; i < file_content.size(); i++){
            if(i % 2 == 0){
                chains.push_back(file_content[i]);
            }else{
                labels.push_back(file_content[i]);
            }
        }
        const int number_of_classes = 5;
        
        int class_cluster_indeces [number_of_classes];
        const std::string classes [number_of_classes] = { "all_alpha", "all_beta", "alpha_slash_beta", "alpha_plus_beta", "small_protein" };
        bool class_index_found [number_of_classes];
        for(int i = 0; i < number_of_classes; i++){
            class_index_found[i] = false;
        }
        
        for(int i = 0; i < labels.size(); i++){
            for(int j = 0; j < number_of_classes; j++){
                if((labels[i] == classes[j])&&(!class_index_found[j])){
                    class_cluster_indeces[j] = i;
                    class_index_found[j] = true;
                }
            }
        }
        
        int number_of_used_indeces = 1000;
        std::vector<int> used_indeces;
        
        for(int i = 0; i < number_of_classes; i++){
            for(int j = 0; j < number_of_used_indeces; j++){
                used_indeces.push_back(class_cluster_indeces[i] + j);
            }
        }
        
        for(int i = 0; i < chains.size(); i++){
            std::string chain = chains[i].substr(0, chains[i].length() - 2);
            std::string auth_seq_id = chains[i].substr(chains[i].length() - 1, chains[i].length());
            auth_seq_ids.push_back(auth_seq_id);
            chains[i] = "./mmCIF_files/" + chain + ".cif";
        }
        
        for(int j = 0; j < used_indeces.size(); j++){
            std::cout << j << std::endl;
            std::ifstream infile(chains[used_indeces[j]]);
            if(infile.good()){
                CIFParser parser;
                //parser.setFilePath(chains[j]);
                parser.setFilePath(chains[used_indeces[j]]);
                Protein protein;
                parser.parseAtomSiteColumns();
                protein.setAllAtoms(parser.parseAtoms());
                /*
                protein.filterAtoms("N", auth_seq_ids[j]);
                protein.filterAtoms("CA", auth_seq_ids[j]);
                protein.filterAtoms("C", auth_seq_ids[j]);
                 */
                
                protein.filterAtoms("N", auth_seq_ids[used_indeces[j]]);
                protein.filterAtoms("CA", auth_seq_ids[used_indeces[j]]);
                protein.filterAtoms("C", auth_seq_ids[used_indeces[j]]);
                
                std::vector<double> phi_angles = protein.getPhiAngles();
                std::vector<double> psi_angles = protein.getPsiAngles();
                
                int iterations = phi_angles.size();
                if(psi_angles.size() < iterations){
                    iterations = psi_angles.size();
                }
                
                if((phi_angles.size() < 100) || (psi_angles.size() < 100)){
                    int additional_phi = 100 - phi_angles.size();
                    int additional_psi = 100 - psi_angles.size();
                    for(int k = 0; k < additional_phi; k++){
                        phi_angles.push_back(0.00);
                    }
                    for(int k = 0; k < additional_psi; k++){
                        psi_angles.push_back(0.00);
                    }
                }
                
                std::string sample = "";
                for(int i = 0; i < 100; i++){
                    sample = sample.append(std::to_string(phi_angles[i]) + ", " + std::to_string(psi_angles[i]) + ", ");
                }
                //sample = sample.append(labels[j]);
                sample = sample.append(labels[used_indeces[j]]);
                samples.push_back(sample);
            }else{
                std::cout << "No file: " << chains[used_indeces[j]] << std::endl;
            }
            
        }
        
        std::string csv_file_name = argv[1];
        csv_file_name = csv_file_name.substr(0, csv_file_name.length() - 4);
        csv_file_name = csv_file_name.append(".csv");
        CSVFile csv_file(csv_file_name);
        csv_file.writeData("./", samples);
        
    }catch(const std::length_error& le){
        std::cerr << le.what() << std::endl;
    }
}
