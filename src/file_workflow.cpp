#include"file_workflow.h"

//---class File

File::File(std::string file_name) : file_name_(file_name){
    
}

void File::setFileName(std::string file_name){
    file_name_ = file_name;
}

std::string File::getFileName() const{
    return file_name_;
}

//---class mmCIFFile

mmCIFFile::mmCIFFile(std::string file_name) : File(file_name){
    
}

std::string mmCIFFile::getProteinName(){
    std::string protein_name = this->getFileName();
    protein_name.erase(0, this->getFileName().size() - 8);
    protein_name.erase(4, 4);
    return protein_name;
}

//---class CSVFile

CSVFile::CSVFile(std::string file_name) : File(file_name){
    
}

CSVFile::CSVFile(std::string file_name, std::vector<double> data_1, std::vector<double> data_2) : File(file_name), data_1_(data_1), data_2_(data_2){
    
}

CSVFile::CSVFile(std::string file_name, std::vector<double> data) : File(file_name), data_1_(data){
    
}

void CSVFile::parseOneDataSet(int data_set_choice){
    std::ifstream input;
    input.open(this->getFileName());
    
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

TXTFile::TXTFile(std::string file_name, std::vector<std::string> data) : File(file_name), data_(data){
    
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
