#include "dataset_generator.h"

DatasetGenerator::DatasetGenerator(std::string input_file) :
    input_file_name_(input_file){

}

void DatasetGenerator::readFile(){
    std::ifstream input;
    input.open(input_file_name_);
    std::string read;
    
    
    while(!input.eof()){
        input >> read;
        input_file_content_.push_back(read);
    }
    input_file_content_.pop_back();
}

void DatasetGenerator::extractChains(){
    for(int i = 0; i < input_file_content_.size(); i++){
        if(i % 2 == 0){
            chains_.push_back(input_file_content_[i]);
        }
    }
}

void DatasetGenerator::extractLabels(){
    for(int i = 0; i < input_file_content_.size(); i++){
        if(i % 2 == 1){
            labels_.push_back(input_file_content_[i]);
        }
    }
}
    
void DatasetGenerator::extractClassClusterIndeces(){
    bool class_index_found [NUMBER_OF_CLASSES];
    for(int i = 0; i < NUMBER_OF_CLASSES; i++){
        class_index_found[i] = false;
    }
    for(int i = 0; i < labels_.size(); i++){
        for(int j = 0; j < NUMBER_OF_CLASSES; j++){
            if((labels_[i] == classes_[j])&&(!class_index_found[j])){
                class_cluster_indeces_[j] = i;
                class_index_found[j] = true;
            }
        }
    }
}

void DatasetGenerator::extractUsedIndeces(int number_of_used_indeces){
    for(int i = 0; i < NUMBER_OF_CLASSES; i++){
        for(int j = 0; j < number_of_used_indeces; j++){
            used_indeces_.push_back(class_cluster_indeces_[i] + j);
        }
    }
}

void DatasetGenerator::extractAuthSeqIds(){
    for(int i = 0; i < chains_.size(); i++){
        std::string auth_seq_id = chains_[i].substr(chains_[i].length() - 1, chains_[i].length());
        auth_seq_ids_.push_back(auth_seq_id);
    }
}

void DatasetGenerator::determineChainsFiles(){
    for(int i = 0; i < chains_.size(); i++){
        std::string chain = chains_[i].substr(0, chains_[i].length() - 2);
        chains_[i] = "./mmCIF_files/" + chain + ".cif";
    }
}

void DatasetGenerator::createSample(int number_of_angles, std::vector<double>& phi, std::vector<double>& psi){
    
    current_sample_ = "";
    for(int i = 0; i < number_of_angles; i++){
//        current_sample_ = current_sample_.append(std::to_string(phi[i]) + ", " + std::to_string(psi[i]) + ", ");
        current_sample_ = current_sample_.append(std::to_string(sin(phi[i])) + ", " + std::to_string(cos(phi[i])) + ", " + std::to_string(sin(psi[i])) + ", " + std::to_string(cos(psi[i])) + ", ");
    }
}

void DatasetGenerator::addSample(int label_index){
    current_sample_ = current_sample_.append(labels_[label_index]);
    samples_.push_back(current_sample_);
}

void DatasetGenerator::writeSamples(){
    std::string csv_file_name = input_file_name_;
    csv_file_name = csv_file_name.substr(0, csv_file_name.length() - 4);
    csv_file_name = csv_file_name.append(".csv");
    CSVFile csv_file(csv_file_name);
    csv_file.writeData("./", samples_);
}

std::vector<int>& DatasetGenerator::getUsedIndeces(){
    return used_indeces_;
}

std::vector<std::string>& DatasetGenerator::getChains(){
    return chains_;
}

std::vector<std::string>& DatasetGenerator::getAuthSeqIds(){
    return auth_seq_ids_;
}
