#ifndef _FILE_WORKFLOW_
#define _FILE_WORKFLOW_
#include <fstream>
#include <vector>
#include <utility>
#include <stdexcept>
#include "biological_structures.h"

class File{
        std::string file_name_;
    public:
        File();
        File(std::string file_name);
        void setFileName(std::string file_name);
        std::string getFileName() const;
};

class mmCIFFile : public File{
    public:
        mmCIFFile();
        mmCIFFile(std::string file_name);
        std::string getProteinName();
};

class CSVFile : public File{
        std::vector<double> data_1_;
        std::vector<double> data_2_;
    public:
        CSVFile();
        CSVFile(std::string file_name);
        CSVFile(std::string file_name, std::vector<double> data);
        CSVFile(std::string file_name, std::vector<double> data_1_, std::vector<double> data_2_);

        void parseOneDataSet(int data_set_choice);
        void setData(std::vector<double> data_set, int data_set_choice);
        void writeData(std::string output_location, std::string data_label, int data_set_choice) const;
        void writeData(std::string output_location, std::string data_1_label, std::string data_2_label) const;
        void appendData(std::string output_location, int data_set_choice) const;

        std::vector<double> getData(int data_set_choice);
};

class TXTFile : public File{
        std::vector<std::string> data_;
    public:
        TXTFile(std::string file_name);
        TXTFile(std::string file_name, std::vector<std::string> data);
        std::vector< std::string > parseData();
        std::vector< std::pair<std::string, std::string> > parsePairedData();
        void setData(std::vector<std::string> data);
        void writeData(std::string output_location);
        void writeData(std::string output_location, std::vector<std::string> data_1, std::vector<std::string> data_2);
        void writeDataVertically(std::string output_location, std::vector<std::string> data_1, std::vector<std::string> data_2);
};

#endif
