#ifndef _PROTEIN_ALIGNMENT_MATRICES_
#define _PROTEIN_ALIGNMENT_MATRICES_
#include <iostream>
#include <string>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <utility>
#include <stdexcept>
#include "statistic_calculator.h"
#include "traceback.h"
#include "sequence_aligner.h"
#define NORTH 360
#define NORTH_WEST 315
#define WEST 270

class BLOSUMMatrix{
    private:
        std::map<std::string, int> matrix_;
        std::map<std::string, int>::iterator it_;
    public:
        void fillMatrix();
        void displayMatrix();
        int returnScore(const std::string& pair) const;
};

class DistanceMatrix{
    private:
        int rows_;
        int columns_;
        std::vector< std::vector<double> > matrix_;
    public:
        DistanceMatrix(std::vector< std::vector<double> > matrix);
    
        void displayMatrix(int neighbour_number_left) const;
    
        void setMatrix(std::vector< std::vector<double> > matrix);
    
        std::vector< std::vector<double> > getMatrix() const;
        double getElement(int index_row, int index_column) const;
        double getMaxValue();
        double getMinValue();
    
        double getMean(StatisticCalculator& statistic_calculator);
        double getStandardDeviation(StatisticCalculator& statistic_calculator);
    
        std::vector< std::vector<double> > calculateScoreMatrix(double mean, double standard_deviation);
};

class DistanceScoreMatrix{
    private:
        int rows_;
        int columns_;
        std::vector< std::vector<double> > matrix_;
    public:
        DistanceScoreMatrix(int rows, int columns);
        DistanceScoreMatrix(std::vector< std::vector<double> > matrix);
    
        void displayMatrix(int neighbour_number_left) const;
    
        void setMatrix(std::vector< std::vector<double> >& score_matrix);
        
        double getMaxValue();
        double getMinValue();
        std::vector< std::vector<double> > getMatrix() const;
        double getScore(std::pair<int, int> pair) const;
    
        void limitScores(int z_score_min, int z_score_max);
        void reverseScoreSigns();
};

class DirectionMatrix{
    private:
        int rows_;
        int columns_;
        std::vector< std::vector<int> > direction_matrix_;
    public:
        DirectionMatrix(int rows, int columns);
        void setDirection(int rows_coordinate, int column_coordinate, int direction);
    
        std::vector< std::vector<int> > returnMatrix() const;
        void displayMatrix() const;
        std::vector< std::pair<int, int> > returnNonZeroCoords() const;
        std::vector<int> returnDirections() const;
};

class ScoringMatrix{
    private:
        int rows_;
        int columns_;
        std::vector< std::vector<double> > scoring_matrix_;
    
        double gap_open_penalty_;
        double gap_ext_penalty_;
    public:
        ScoringMatrix(int rows, int columns, double gap_open_penalty, double gap_ext_penalty);
    
        void setGapOpenPenalty(double gap_open_penalty);
        void setGapExtPenalty(double gap_ext_penalty);
    
        int getRows() const;
        int getColumns() const;
        std::vector< std::vector<double> > getScoringMatrix() const;
        double getGapOpenPenalty() const;
        double getGapExtPenalty() const;

        void printScoringMatrix() const;
        void fillWithGapPenalties();
    
        std::pair<int, int> getPair(int p_subunit, int q_subunit);
        double getNorthResult(double north_arg, double penalty);
        double getNorthWestResult(double north_west_arg, double score);
        double getWestResult(double west_arg, double score);
        double getInitialMaxOper(double north_arg, double north_west_arg, double west_arg);
        double getMaxOperationValue(double max_oper, double north, double north_west, double west);
        
        void fillWithScores(char algorithm_choice, std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q, const DistanceScoreMatrix& matrix);
    
        void setTracebackStart(char algorithm_choice, double& start_value, int& start_row_coord, int& start_col_coord);
        void setTracebackInitialMaxOper(char algorithm_choice, double& max_oper);
        void tracebackNorthOperation(int row_coord, int col_coord, double& north_value, int& north_row_coord, int& north_col_coord, char& zero_coordinate);
        void tracebackNorthWestOperation(int row_coord, int col_coord, double& north_west_value, int& north_west_row_coord, int& north_west_col_coord);
        void tracebackWestOperation(int row_coord, int col_coord, double& west_value, int& west_row_coord, int& west_col_coord, char& zero_coordinate);
        void tracebackMaxOper(double& max_oper, int& max_value_row_coord, int& max_value_col_coord, double value, int row_coord, int col_coord);
        char getIndicator(double& max_oper, int& max_value_row_coord, int& max_value_col_coord, double west_value, int west_row_coord, int west_col_coord, double north_value, int north_row_coord, int north_col_coord, double north_west_value, int north_west_row_coord, int north_west_col_coord);
    
        void setTracebackDirections(char zero_coordinate, char indicator, DirectionMatrix& direction_matrix, int west_row_coord, int west_col_coord, int north_row_coord, int north_col_coord, int north_west_row_coord, int north_west_col_coord);
    
        DirectionMatrix getDirectionMatrix(char algorithm_choice, std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q);
    
        void algorithmSmithWaterman(std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q, DistanceScoreMatrix& matrix);
        void algorithmNeedlemanWunsch(std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q, DistanceScoreMatrix& matrix);
        
};
#endif



