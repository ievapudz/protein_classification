#ifndef _PREPARATORY_PHASE_
#define _PREPARATORY_PHASE_
#include "./src/constants.hpp"
#include "./src/file_workflow.h"
#include "./src/cif_parser.h"
#include "./src/biological_structures.h"
#include <vector>
#include <iostream>
#include <algorithm>


class PreparatoryPhase{
public:
    Constants constants_;
    TXTFile aminoacid_codes_file_;
    std::vector<std::pair<std::string, std::string> > aminoacid_codes_;
    TXTFile protein_chains_list_;
    std::vector< std::string > protein_chains_;
    std::vector<std::string> auth_asym_ids_;
    mmCIFFile p_mmCIF_file_;
    mmCIFFile q_mmCIF_file_;
    CIFParser p_parser_;
    CIFParser q_parser_;
    Protein p_protein_;
    Protein q_protein_;
    CSVFile distance_file_;
    
    PreparatoryPhase(std::string protein_chains_list_file_name);
    void setProtein(char protein, int index);
    void setDistanceFile(std::string distance_file_name, int substructure_length);
    void run(int index_p, int index_q);
};

#endif
