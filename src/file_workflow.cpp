#include"file_workflow.h"

//---class File
File::File() : file_name_("-"){

}

File::File(std::string file_name) : file_name_(file_name){

}

void File::setFileName(std::string file_name){
    file_name_ = file_name;
}

std::string File::getFileName() const{
    return file_name_;
}

//---class mmCIFFile

mmCIFFile::mmCIFFile() : File(){

}

mmCIFFile::mmCIFFile(std::string file_name) : File(file_name){

}

std::string mmCIFFile::getProteinName(){
    std::string protein_name = this->getFileName();
    std::size_t found = protein_name.find_last_of("/\\");
    std::string buffer = protein_name.substr(found+1);
    buffer = buffer.erase(buffer.size()-4, buffer.size());
    protein_name = buffer;
    return protein_name;
}

//---class CSVFile

CSVFile::CSVFile() : File(){

}

CSVFile::CSVFile(std::string file_name) : File(file_name){

}

CSVFile::CSVFile(std::string file_name, std::vector<double> data_1, std::vector<double> data_2) : File(file_name), data_1_(data_1), data_2_(data_2){

}

CSVFile::CSVFile(std::string file_name, std::vector<double> data) : File(file_name), data_1_(data){

}

void CSVFile::parseOneDataSet(int data_set_choice){
    std::ifstream input;
    input.open(this->getFileName());

    switch (data_set_choice) {
        case 1:
            if(data_1_.size()>0){
                data_1_.clear();
            }
            break;
        case 2:
            if(data_2_.size()>0){
                data_2_.clear();
            }
            break;
    }

    std::string reading_label;
    double reading_data;

    input >> reading_label;

    do{
        input >> reading_data;
        if(data_set_choice == 1)
            data_1_.push_back(reading_data);
        if(data_set_choice == 2)
            data_2_.push_back(reading_data);
    }while(!input.eof());

    if(data_set_choice == 1)
        data_1_.pop_back();
    if(data_set_choice == 2)
        data_2_.pop_back();

}

void CSVFile::setData(std::vector<double> data_set, int data_set_choice){
    switch (data_set_choice) {
        case 1:{
            data_1_ = data_set;
            break;
        }
        case 2:{
            data_2_ = data_set;
            break;
        }
        default:{
            std::invalid_argument ia("Invalid argument for data set choice in CSVFile::writeData.");
            throw ia;
            break;
        }
    }
}

void CSVFile::writeData(std::string output_location, std::vector<std::string> data_set) const{
    std::ofstream output;
    output_location.append(this->getFileName());
    output.open(output_location);
    
    for(int i = 0; i < data_set.size(); i++){
        output << data_set[i] << std::endl;
    }
}

void CSVFile::writeData(std::string output_location, std::string data_label, int data_set_choice) const{
    std::ofstream output;
    output_location.append(this->getFileName());
    output.open(output_location);

    output << data_label << std::endl;

    switch (data_set_choice) {
        case 1:{
            for(int i = 0; i < data_1_.size(); i++){
                output << data_1_[i] << std::endl;
            }
            break;
        }
        case 2:{
            for(int i = 0; i < data_2_.size(); i++){
                output << data_2_[i] << std::endl;
            }
            break;
        }
        default:{
            std::invalid_argument ia("Invalid argument for data set choice in CSVFile::writeData.");
            throw ia;
            break;
        }
    }
}

void CSVFile::writeData(std::string output_location, std::string data_1_label, std::string data_2_label) const{
    std::ofstream output;
    output_location.append(this->getFileName());
    output.open(output_location);

    int min = data_1_.size();
    if(data_2_.size() < data_1_.size())
        min = data_2_.size();

    output << data_1_label << ", " << data_2_label << std::endl;

    for(int i = 0; i < min; i++){
        output << data_1_[i] << ", " << data_2_[i] << std::endl;
    }
}

void CSVFile::appendData(std::string output_location, int data_set_choice) const{
    std::ofstream output;
    output_location.append(this->getFileName());
    output.open(output_location, std::ios_base::app);

    switch (data_set_choice) {
        case 1:{
            for(int i = 0; i < data_1_.size(); i++){
                output << data_1_[i] << std::endl;
            }
            break;
        }
        case 2:{
            for(int i = 0; i < data_2_.size(); i++){
                output << data_2_[i] << std::endl;
            }
            break;
        }
        default:{
            std::invalid_argument ia("Invalid argument for data set choice in CSVFile::writeData.");
            throw ia;
            break;
        }
    }
}

std::vector<double> CSVFile::getData(int data_set_choice){
    switch (data_set_choice) {
        case 1:{
            return data_1_;
            break;
        }
        case 2:{
            return data_2_;
            break;
        }
        default:{
            std::invalid_argument ia("Invalid argument for data set choice in CSVFile::getData.");
            throw ia;
            break;
        }
    }
}

//---class TXTFile

TXTFile::TXTFile(std::string file_name) : File(file_name){

}

TXTFile::TXTFile(std::string file_name, std::vector<std::string> data) : File(file_name), data_(data){

}

std::vector< std::string > TXTFile::parseData(){
    std::vector< std::string > data;
    std::ifstream input;
    input.open(this->getFileName());
    std::string reading_1;
    while(!input.eof()){
        input >> reading_1;
        data.push_back(reading_1);
    }
    data.pop_back();
    return data;
}

std::vector< std::pair<std::string, std::string> > TXTFile::parsePairedData(){
    std::vector< std::pair<std::string, std::string> > paired_data;

    std::ifstream input;
    input.open(this->getFileName());
    std::string reading_1, reading_2;
    while(!input.eof()){
        input >> reading_1 >> reading_2;
        std::pair<std::string, std::string> pair = std::make_pair(reading_1, reading_2);
        paired_data.push_back(pair);
    }

    return paired_data;
}

void TXTFile::setData(std::vector<std::string> data){
    data_ = data;
}

void TXTFile::writeData(std::string output_location){
    output_location.append(this->getFileName());
    std::ofstream output;
    output.open(output_location);

    for(int i = 0; i < data_.size(); i++){
        output << data_[i] << std::endl;
    }
}

void TXTFile::writeData(std::string output_location, std::vector< std::string > labels, std::vector< std::vector<double> > data){

    output_location.append(this->getFileName());
    std::ofstream output;
    output.open(output_location);
    
    for(int i = 0; i < labels.size(); i++){
        output << labels[i] << "\n";
    }

    for(int i = 0; i < data.size(); i++){
        for(int j = 0; j < data[0].size(); j++){
            output << data[i][j] << " ";
        }
        output << "\n";
    }
}

void TXTFile::writeData(std::string output_location, std::vector<std::string> data_1, std::vector<std::string> data_2){

    output_location.append(this->getFileName());
    std::ofstream output;
    output.open(output_location);

    for(int i = 0; i < data_1.size(); i++){
        output << data_1[i];
    }
    output << std::endl;
    for(int i = 0; i < data_1.size(); i++){
        output << data_2[i];
    }
}

void TXTFile::writeDataVertically(std::string output_location, std::vector<std::string> data_1, std::vector<std::string> data_2){
    output_location.append(this->getFileName());
    std::ofstream output;
    output.open(output_location);

    int min = data_1.size();
    if(min > data_2.size()){
      min = data_2.size();
    }

    for(int i = 0; i < min; i++){
        output << data_1[i] << " " << data_2[i] << std::endl;
    }
}
