#include "preparatory_phase.h"

PreparatoryPhase::PreparatoryPhase(std::string protein_chains_list_file_name):
    constants_("./mmCIF_files/", -4, 4, 4.0, 1.0, 1, 7),
    aminoacid_codes_file_("aminoacid_codes.txt"),
    protein_chains_list_(protein_chains_list_file_name){

        aminoacid_codes_ = aminoacid_codes_file_.parsePairedData();
        protein_chains_ = protein_chains_list_.parseData();

        for(int i = 0; i < protein_chains_.size(); i++){
            auth_asym_ids_.push_back(std::string(1, protein_chains_[i][protein_chains_[i].size()-1]));
            //auth_asym_ids_.push_back(std::string(1, protein_chains_[i][5]));
            transform(protein_chains_[i].begin(), protein_chains_[i].end(), protein_chains_[i].begin(), ::tolower);
            protein_chains_[i] = protein_chains_[i].substr(0, protein_chains_[i].size() - 2);
            std::string path = constants_.cifFilePath();
            path.append(protein_chains_[i]+".cif");
            protein_chains_[i] = path;
        }
}

void PreparatoryPhase::setProtein(char protein, int index){
    if( protein == 'P'){
        p_mmCIF_file_.setFileName(protein_chains_[index]);
        p_parser_.setFilePath(protein_chains_[index]);
        p_parser_.parseAtomSiteColumns();
        p_protein_.setName(p_mmCIF_file_.getProteinName());
        p_protein_.setAllAtoms(p_parser_.parseAtoms());
        p_protein_.filterAtoms("N", auth_asym_ids_[index]);
        p_protein_.filterAtoms("CA", auth_asym_ids_[index]);
        p_protein_.filterAtoms("C", auth_asym_ids_[index]);
        p_protein_.setAminoacidSequence(aminoacid_codes_);
    }else{
        q_mmCIF_file_.setFileName(protein_chains_[index]);
        q_parser_.setFilePath(protein_chains_[index]);
        q_parser_.parseAtomSiteColumns();
        q_protein_.setName(q_mmCIF_file_.getProteinName());
        q_protein_.setAllAtoms(q_parser_.parseAtoms());
        q_protein_.filterAtoms("N", auth_asym_ids_[index]);
        q_protein_.filterAtoms("CA", auth_asym_ids_[index]);
        q_protein_.filterAtoms("C", auth_asym_ids_[index]);
        q_protein_.setAminoacidSequence(aminoacid_codes_);
    }
}

void PreparatoryPhase::setDistanceFile(std::string distance_file_name, int substructure_length){
    distance_file_name.append(std::to_string(substructure_length));
    distance_file_name.append("_length_distances.csv");

    distance_file_.setFileName(distance_file_name);
    distance_file_.parseOneDataSet(1);
}

std::string PreparatoryPhase::getProteinChain(int protein_chains_list_index){
    std::string protein_name = protein_chains_[protein_chains_list_index];
    std::size_t found = protein_name.find_last_of("/\\");
    std::string buffer = protein_name.substr(found+1);
    buffer = buffer.erase(buffer.size()-4, buffer.size());
    buffer.append("_" + auth_asym_ids_[protein_chains_list_index]);
    protein_name = buffer;
    return protein_name;
}

void PreparatoryPhase::run(int index_p, int index_q){
    p_protein_.reset();
    q_protein_.reset();
    this->setProtein('P', index_p);
    this->setProtein('Q', index_q);
}
