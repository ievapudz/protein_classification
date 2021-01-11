#include <algorithm>
#include "./src/scoring_matrix.h"
#include "./src/distance_score_matrix.h"
#include "./src/distance_matrix.h"
#include "./src/biological_structures.h"
#include "./src/cif_parser.h"
#include "./src/file_workflow.h"
#include "./src/block_distance.h"
#include "./src/constants.hpp"
#include "preparatory_phase.h"
#include "calculation_phase.h"
#include "representation_phase.h"

int main(int argc, const char * argv[]){

    try{
        PreparatoryPhase prep_phase(argv[1]);
        std::vector<std::string> protein_chains;
        
        double identity_scores[prep_phase.protein_chains_.size()][prep_phase.protein_chains_.size()];
        for(int i = 0; i < prep_phase.protein_chains_.size(); i++){
            for(int j = 0; j < prep_phase.protein_chains_.size(); j++){
                identity_scores[i][j] = 0.0;
            }
        }
        
        for(int j = 0; j < prep_phase.protein_chains_.size(); j++){
            protein_chains.push_back(prep_phase.protein_chains_[j]);
            for(int k = 0; k < prep_phase.protein_chains_.size(); k++){
                prep_phase.run(j, k);
                for(int i = prep_phase.constants_.minSubstructureLength(); i <= prep_phase.constants_.maxSubstructureLength(); i++){
                    CalculationPhase calc_phase(&prep_phase, i);
                    calc_phase.run();
                    RepresentationPhase repr_phase(&calc_phase);
                    repr_phase.representNumeralAlignment();
                    
                    std::string str = std::to_string(i);
                    if(str == argv[2]){
                        identity_scores[j][k] = repr_phase.representIdentityScore().first;
                    }
                }
            }
        }
        
        std::cout << "\nAligned proteins: " << std::endl;
        for(int i = 0; i < protein_chains.size(); i++){
            std::cout << protein_chains[i] << std::endl;
        }
        
        for(int i = 0; i < prep_phase.protein_chains_.size(); i++){
            for(int j = 0; j < prep_phase.protein_chains_.size(); j++){
                std::cout << identity_scores[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }catch(const std::length_error& le){
        std::cerr << le.what() << std::endl;
    }
    catch(const std::invalid_argument& ia){
        std::cerr << ia.what() << std::endl;
    }
}
