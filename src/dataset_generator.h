#ifndef _DATASET_GENERATOR_
#define _DATASET_GENERATOR_
#define NUMBER_OF_CLASSES 5
#include <cmath>
#include <vector>
#include <fstream>
#include "file_workflow.h"

class DatasetGenerator{
    private:
        const std::string classes_ [NUMBER_OF_CLASSES] = { "all_alpha", "all_beta", "alpha_slash_beta", "alpha_plus_beta", "small_protein" };
    
        std::string input_file_name_;
        std::vector<std::string> input_file_content_;
        std::vector<std::string> chains_;
        std::vector<std::string> auth_seq_ids_;
        std::vector<std::string> labels_;
    
        std::string current_sample_;
        std::vector<std::string> samples_;
    
        int class_cluster_indeces_ [NUMBER_OF_CLASSES];
        std::vector<int> used_indeces_;
    public:
        DatasetGenerator(std::string input_file);
        void readFile();
        void extractChains();
        void extractLabels();
        void extractClassClusterIndeces();
        void extractUsedIndeces(int number_of_used_indeces);
        void extractAuthSeqIds();
        void determineChainsFiles();
        void createSample(int number_of_angles, std::vector<double>& phi, std::vector<double>& psi);
        void addSample(int label_index);
        void writeSamples();
        
        std::vector<int>& getUsedIndeces();
        std::vector<std::string>& getChains();
        std::vector<std::string>& getAuthSeqIds();
    
    
};

#endif
