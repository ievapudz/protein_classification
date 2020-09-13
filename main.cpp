#include "./src/scoring_matrix.h"
#include "./src/distance_score_matrix.h"
#include "./src/distance_matrix.h"
#include "./src/biological_structures.h"
#include "./src/cif_parser.h"
#include "./src/file_workflow.h"
#include "./src/block_distance.h"

int main(int argc, const char * argv[]){
    
    try{
        TXTFile aminoacid_code_file("aminoacid_codes.txt");
        std::vector<std::pair<std::string, std::string> > aminoacid_codes = aminoacid_code_file.parsePairedData();
        
        mmCIFFile file_1(argv[1]);
        mmCIFFile file_2(argv[3]);
        
        CIFParser cp_1;
        cp_1.setFilePath(file_1.getFileName());
        cp_1.parseAtomSiteColumns();
        
        CIFParser cp_2;
        cp_2.setFilePath(file_2.getFileName());
        cp_2.parseAtomSiteColumns();
        
        Protein p_protein(file_1.getProteinName(), cp_1.parseAtoms());
        
        p_protein.filterAtoms("N", argv[2]);
        p_protein.filterAtoms("CA", argv[2]);
        p_protein.filterAtoms("C", argv[2]);
        
        std::vector<std::string> p_aminoacid_sequence = p_protein.getAminoacidSequence(aminoacid_codes);
        
        Protein q_protein(file_2.getProteinName(), cp_2.parseAtoms());
        
        q_protein.filterAtoms("N", argv[4]);
        q_protein.filterAtoms("CA", argv[4]);
        q_protein.filterAtoms("C", argv[4]);
        
        std::vector<std::string> q_aminoacid_sequence = q_protein.getAminoacidSequence(aminoacid_codes);
        
        int z_score_min = -4;
        int z_score_max = 4;
        
        BlockDistanceCalculator bdc(p_protein, q_protein);
        
        for(int substructure_length = 1; substructure_length <= 7; substructure_length++){
            std::string distance_file_name = "./csv_distance_files/";
            distance_file_name.append(std::to_string(substructure_length));
            distance_file_name.append("_length_distances.csv");
            
            CSVFile distance_file(distance_file_name);
            distance_file.parseOneDataSet(1);
            
            StatisticCalculator sc;
            sc.setData(distance_file.getData(1));
            
            double mean = sc.getMean();
            double standard_deviation = sc.getStandardDeviation();
            
            bdc.setSubstructureLength(substructure_length);
            
            DistanceMatrix dm(bdc.getSuccessiveDistances(1));
            
            DistanceScoreMatrix dsm(dm.calculateScoreMatrix(mean, standard_deviation));
            
            dsm.limitScores(z_score_min, z_score_max);
            dsm.reverseScoreSigns();
            
            std::vector<int> chain_P = bdc.getNumberedSubunitChain(1);
            std::vector<int> chain_Q = bdc.getNumberedSubunitChain(2);
            
            double gap_open_penalty = z_score_max;
            double gap_ext_penalty = 1;
            
            ScoringMatrix sm(chain_P.size() + 1, chain_Q.size() + 1, gap_open_penalty, gap_ext_penalty);
            sm.algorithmNeedlemanWunsch(chain_P, chain_Q, dsm, substructure_length, p_aminoacid_sequence, q_aminoacid_sequence, 2);
        }
        
    }catch(const std::length_error& le){
        std::cerr << le.what() << std::endl;
    }
    catch(const std::invalid_argument& ia){
        std::cerr << ia.what() << std::endl;
    }
    
}
