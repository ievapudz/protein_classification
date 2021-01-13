#include "identity_score_table.h"

IdentityScoreTable::IdentityScoreTable(int number_of_proteins):
    table_(number_of_proteins, std::vector<double>(number_of_proteins, 0.0)){
    
}

void IdentityScoreTable::setSubstructureLength(int substructure_length){
    substructure_length_ = substructure_length;
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

int IdentityScoreTable::getSubstructureLength() const{
    return substructure_length_;
}

std::vector< std::string > IdentityScoreTable::getProteins() const{
    return proteins_;
}

std::vector< std::vector<double> > IdentityScoreTable::getTable() const{
    return table_;
}

double IdentityScoreTable::getMaxScoreInRow(int index_row){
    double max_score = table_[index_row][0];
    std::string label = SCOP_class_labels[0];
    for(int i = 0; i < table_.size() - 1; i++){
        if(table_[index_row][i] > max_score){
            max_score = table_[index_row][i];
            label = SCOP_class_labels[i];
        }
    }
    std::cout << label << std::endl;
    return max_score;
}

double IdentityScoreTable::getMaxScoreInColumn(int index_column){
    double max_score = table_[0][index_column];
    std::string label = SCOP_class_labels[0];
    for(int i = 0; i < table_.size() - 1; i++){
        if(table_[i][index_column] > max_score){
            max_score = table_[i][index_column];
            label = SCOP_class_labels[i];
        }
    }
    std::cout << label << std::endl;
    return max_score;
}

std::string IdentityScoreTable::getLabelByRow(int index_row){
    double max_score = table_[index_row][0];
    std::string label = SCOP_class_labels[0];
    for(int i = 0; i < table_.size() - 1; i++){
        if(table_[index_row][i] > max_score){
            max_score = table_[index_row][i];
            label = SCOP_class_labels[i];
        }
    }
    return label;
}

std::string IdentityScoreTable::getLabelByColumn(int index_column){
    double max_score = table_[0][index_column];
    std::string label = SCOP_class_labels[0];
    for(int i = 0; i < table_.size() - 1; i++){
        if(table_[i][index_column] > max_score){
            max_score = table_[i][index_column];
            label = SCOP_class_labels[i];
        }
    }
    return label;
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

void IdentityScoreTable::addLabels(){
    for(int i = 0; i < proteins_.size(); i++){
        if(i >= 5){
            break;
        }
        proteins_[i].append(" " + SCOP_class_labels[i]);
    }
}

void IdentityScoreTable::printTableToFile(std::string file_name) const{
    TXTFile identity_score_table_file(file_name);
    identity_score_table_file.writeData("./identity_scores/", proteins_, table_);
}
