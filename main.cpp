#include "./src/scoring_matrix.h"
#include "./src/distance_score_matrix.h"
#include "./src/distance_matrix.h"
#include "./src/biological_structures.h"
#include "./src/cif_parser.h"
#include "./src/file_workflow.h"
#include "./src/block_distance.h"
//#define SUBSTRUCTURE_LENGTH 7

int main(int argc, const char * argv[]){
    
    try{
        for(int SUBSTRUCTURE_LENGTH = 1; SUBSTRUCTURE_LENGTH <= 7; SUBSTRUCTURE_LENGTH++){
            std::string distance_file_name = "./csv_distance_files/";
            distance_file_name.append(std::to_string(SUBSTRUCTURE_LENGTH));
            distance_file_name.append("_length_distances.csv");
            
            CSVFile distance_file(distance_file_name);
            distance_file.parseOneDataSet(1);
            
            StatisticCalculator sc;
            sc.setData(distance_file.getData(1));
            
            double mean = sc.getMean();
            double standard_deviation = sc.getStandardDeviation();
            
            mmCIFFile file_1(argv[1]);
            mmCIFFile file_2(argv[2]);
            
            CIFParser cp_1;
            cp_1.setFilePath(file_1.getFileName());
            cp_1.parseAtomSiteColumns();
            
            CIFParser cp_2;
            cp_2.setFilePath(file_2.getFileName());
            cp_2.parseAtomSiteColumns();
            
            Protein p_protein(file_1.getProteinName(), cp_1.parseAtoms());
            p_protein.filterAtoms("N");
            p_protein.filterAtoms("CA");
            p_protein.filterAtoms("C");
            
            Protein q_protein(file_2.getProteinName(), cp_2.parseAtoms());
            q_protein.filterAtoms("N");
            q_protein.filterAtoms("CA");
            q_protein.filterAtoms("C");
            
            int z_score_min = -4;
            int z_score_max = 4;
            
            BlockDistanceCalculator bdc(p_protein, q_protein);
            bdc.setSubstructureLength(SUBSTRUCTURE_LENGTH);
            
            DistanceMatrix dm(bdc.getSuccessiveDistances(1));
            
            DistanceScoreMatrix dsm(dm.calculateScoreMatrix(mean, standard_deviation));
            dsm.limitScores(z_score_min, z_score_max);
            dsm.reverseScoreSigns();
            
            std::vector<int> chain_P = bdc.getNumberedSubunitChain(1);
            std::vector<int> chain_Q = bdc.getNumberedSubunitChain(2);
            
            double gap_ext_penalty = 0.0;
            double gap_open_penalty = z_score_max;
            
            ScoringMatrix sm(chain_P.size() + 1, chain_Q.size() + 1, gap_open_penalty, gap_ext_penalty);
            sm.algorithmNeedlemanWunsch(chain_P, chain_Q, dsm, SUBSTRUCTURE_LENGTH);
        }
        
    }catch(const std::length_error& le){
        std::cerr << le.what() << std::endl;
    }
    catch(const std::invalid_argument& ia){
        std::cerr << ia.what() << std::endl;
    }
    
}
