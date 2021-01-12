#include <algorithm>
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
        PreparatoryPhase prep_phase(argv[1]);
        
        std::vector<std::string> protein_chains;
        
        std::vector< IdentityScoreTable > tables(prep_phase.constants_.maxSubstructureLength() - prep_phase.constants_.minSubstructureLength() + 1, IdentityScoreTable(prep_phase.protein_chains_.size()));
        
        for(int j = 0; j < prep_phase.protein_chains_.size(); j++){
            protein_chains.push_back(prep_phase.getProteinChain(j));
            for(int k = 0; k < prep_phase.protein_chains_.size(); k++){
                prep_phase.run(j, k);
                int tables_index = 0;
                for(int i = prep_phase.constants_.minSubstructureLength(); i <= prep_phase.constants_.maxSubstructureLength(); i++){
                    CalculationPhase calc_phase(&prep_phase, i);
                    calc_phase.run();
                    RepresentationPhase repr_phase(&calc_phase);
                    repr_phase.representNumeralAlignment();
                    
                    tables[tables_index].setSubstructureLength(i);
                    tables[tables_index].setAt(j, k, repr_phase.representIdentityScore().first);
                    tables_index++;
                }
            }
        }
        std::string protein_chains_list_name(argv[1]);
        for(IdentityScoreTable table : tables){
            table.setProteins(protein_chains);
            table.printTableToFile("identity_scores_" + protein_chains_list_name.substr(0, protein_chains_list_name.size() - 4) + "_" + std::to_string(table.getSubstructureLength()) + ".txt");
        }
        
    }catch(const std::length_error& le){
        std::cerr << le.what() << std::endl;
    }
    catch(const std::invalid_argument& ia){
        std::cerr << ia.what() << std::endl;
    }
}
