#ifndef _IDENTITY_SCORE_TABLE_
#define _IDENTITY_SCORE_TABLE_
#include <iostream>
#include <vector>

class IdentityScoreTable{
    private:
    std::vector< std::string > proteins_;
        std::vector< std::vector<double> > table_;
    public:
        IdentityScoreTable(int number_of_proteins);
    
        void setProteins(std::vector< std::string > proteins);
        void setTable(std::vector< std::vector<double> > table);
        void setAt(int index_row, int index_column, double element);
        std::vector< std::vector<double> > getTable() const;
    
        void printProteins() const;
        void printTable() const;
        void printTableToFile(std::string file_name) const;
};

#endif
