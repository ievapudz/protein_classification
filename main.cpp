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

int main(int argc, const char * argv[]){
    
    try{
        /*
        PreparatoryPhase prep_phase(argv[1]);
        prep_phase.run(0, 1);
        CalculationPhase calc_phase(&prep_phase, 6);
        calc_phase.run();*/
        
        Constants constants("./mmCIF_files/", -4, 4, 4.0, 1.0, 1, 7);
        
        TXTFile aminoacid_code_file("aminoacid_codes.txt");
        std::vector<std::pair<std::string, std::string> > aminoacid_codes = aminoacid_code_file.parsePairedData();
        
        TXTFile protein_chain_file(argv[1]);
        std::vector< std::string > protein_chains = protein_chain_file.parseData();
        
        std::vector< std::string > auth_asym_ids;
        
        for(int i = 0; i < protein_chains.size(); i++){
            auth_asym_ids.push_back(std::string(1, protein_chains[i][5]));
            transform(protein_chains[i].begin(), protein_chains[i].end(), protein_chains[i].begin(), ::tolower);
            protein_chains[i] = protein_chains[i].substr(0, 4);
            std::string path = constants.cifFilePath();
            path.append(protein_chains[i]+".cif");
            protein_chains[i] = path;
            std::cout << protein_chains[i] << std::endl;
        }
        
        for(int i = 0; i < protein_chains.size() - 1; i++){
            long int start = time(NULL);
            
            mmCIFFile file_1(protein_chains[i]);
            mmCIFFile file_2(protein_chains[i+1]);
            
            CIFParser cp_1;
            cp_1.setFilePath(file_1.getFileName());
            cp_1.parseAtomSiteColumns();
            
            CIFParser cp_2;
            cp_2.setFilePath(file_2.getFileName());
            cp_2.parseAtomSiteColumns();
            
            Protein p_protein(file_1.getProteinName(), cp_1.parseAtoms());
            
            p_protein.filterAtoms("N", auth_asym_ids[i]);
            p_protein.filterAtoms("CA", auth_asym_ids[i]);
            p_protein.filterAtoms("C", auth_asym_ids[i]);
            p_protein.setAminoacidSequence(aminoacid_codes);
            
            Protein q_protein(file_2.getProteinName(), cp_2.parseAtoms());
            
            q_protein.filterAtoms("N", auth_asym_ids[i+1]);
            q_protein.filterAtoms("CA", auth_asym_ids[i+1]);
            q_protein.filterAtoms("C", auth_asym_ids[i+1]);
            q_protein.setAminoacidSequence(aminoacid_codes);
            
            BlockDistanceCalculator bdc(p_protein, q_protein);
            
            for(int substructure_length = constants.minSubstructureLength(); substructure_length <= constants.maxSubstructureLength(); substructure_length++){
                std::string distance_file_name = "./csv_distance_files/";
                distance_file_name.append(std::to_string(substructure_length));
                distance_file_name.append("_length_distances.csv");
                
                CSVFile distance_file(distance_file_name);
                distance_file.parseOneDataSet(1);
                
                // Calculations
                StatisticCalculator sc;
                sc.setData(distance_file.getData(1));
                
                double mean = sc.getMean();
                double standard_deviation = sc.getStandardDeviation();
                
                bdc.setSubstructureLength(substructure_length);
                
                DistanceMatrix dm(bdc.getSuccessiveDistances(1));
                
                DistanceScoreMatrix dsm(dm.calculateScoreMatrix(mean, standard_deviation));
                
                dsm.limitScores(constants.zScoreMin(), constants.zScoreMax());
                dsm.reverseScoreSigns();
                
                p_protein.setSubunitChain(bdc.getNumberedSubunitChain(1));
                q_protein.setSubunitChain(bdc.getNumberedSubunitChain(2));
                
                ScoringMatrix sm(p_protein, q_protein, constants.gapOpenPenalty(), constants.gapExtPenalty());
                sm.algorithmNeedlemanWunsch(dsm, substructure_length, p_protein, q_protein, 2);
                
            }
            std::cout << "Duration: " << time(NULL) - start << " seconds" << std::endl;
        }
        
    }catch(const std::length_error& le){
        std::cerr << le.what() << std::endl;
    }
    catch(const std::invalid_argument& ia){
        std::cerr << ia.what() << std::endl;
    }
    
}
