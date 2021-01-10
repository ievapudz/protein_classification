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
        for(int j = 0; j < prep_phase.protein_chains_.size(); j++){
            for(int k = 0; k < prep_phase.protein_chains_.size(); k++){
                prep_phase.run(j, k);
                for(int i = prep_phase.constants_.minSubstructureLength(); i <= prep_phase.constants_.maxSubstructureLength(); i++){
                    CalculationPhase calc_phase(&prep_phase, i);
                    calc_phase.run();
                    RepresentationPhase repr_phase(&calc_phase);
                    repr_phase.representNumeralAlignment();
                    //repr_phase.representIdentityScore();
                }
            }
        }
    }catch(const std::length_error& le){
        std::cerr << le.what() << std::endl;
    }
    catch(const std::invalid_argument& ia){
        std::cerr << ia.what() << std::endl;
    }
}
