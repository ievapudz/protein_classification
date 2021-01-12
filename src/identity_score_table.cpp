#include "identity_score_table.h"

IdentityScoreTable::IdentityScoreTable(int number_of_proteins):
    table_(number_of_proteins, std::vector<double>(number_of_proteins, 0.0)){
    
}

void IdentityScoreTable::setProteins(std::vector< std::string > proteins){
    proteins_ = proteins;
}

void IdentityScoreTable::setTable(std::vector< std::vector<double> > table){
    table_ = table;
}

void IdentityScoreTable::setAt(int index_row, int index_column, double element){
    table_[index_row][index_column] = element;
}

std::vector< std::vector<double> > IdentityScoreTable::getTable() const{
    return table_;
}

void IdentityScoreTable::printProteins() const{
    for(int i = 0; i < proteins_.size(); i++){
        std::cout << proteins_[i] << std::endl;
    }
}

void IdentityScoreTable::printTable() const{
    for(int i = 0; i < table_.size(); i++){
        for(int j = 0; j < table_[0].size(); j++){
            std::cout << table_[i][j] << " ";
        }
        std::cout << std::endl;
    }
}
