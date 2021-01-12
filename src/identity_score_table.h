#ifndef _IDENTITY_SCORE_TABLE_
#define _IDENTITY_SCORE_TABLE_
#include "file_workflow.h"
#include <iostream>
#include <vector>

class IdentityScoreTable{
    private:
        const std::string SCOP_class_labels[5] = { "all_alpha", "all_beta", "alpha/beta", "alpha+beta", "small_protein" };
        int substructure_length_;
        std::vector< std::string > proteins_;
        std::vector< std::vector<double> > table_;
    public:
        IdentityScoreTable(int number_of_proteins);
    
        void setSubstructureLength(int substructure_length);
        void setProteins(std::vector< std::string > proteins);
        void setTable(std::vector< std::vector<double> > table);
        void setAt(int index_row, int index_column, double element);
        
        int getSubstructureLength() const;
        std::vector< std::string > getProteins() const;
        std::vector< std::vector<double> > getTable() const;
        double getMaxScoreInRow(int index_row);
        double getMaxScoreInColumn(int index_column);
        std::string getLabelByRow(int index_row);
        std::string getLabelByColumn(int index_column);
    
        void printProteins() const;
        void printTable() const;
        void printTableToFile(std::string file_name) const;
};

#endif
