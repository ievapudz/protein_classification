#ifndef _SCORING_MATRIX_
#define _SCORING_MATRIX_
#include <iostream>
#include <vector>
#include <utility>
#include <iomanip>
#include <stdexcept>
#include "biological_structures.h"
#include "traceback.h"
#include "distance_matrix.h"
#include "distance_score_matrix.h"
#include "direction_matrix.h"
#include "sequence_aligner.h"
#include "file_workflow.h"
#include "penalty_decision_matrix.h"
#define NORTH 360
#define NORTH_WEST 315
#define WEST 270

class ScoringMatrix{
    private:
        int rows_;
        int columns_;
        std::vector< std::vector<double> > scoring_matrix_;
    
        double gap_open_penalty_;
        double gap_ext_penalty_;
    public:
        ScoringMatrix(int rows, int columns, double gap_open_penalty, double gap_ext_penalty);
        ScoringMatrix(Protein p_protein, Protein q_protein, double gap_open_penalty, double gap_ext_penalty);
    
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
        double getMaxOperationValue(double max_oper, double north, double north_west, double west, int& penalty_decision_matrix_element);
    
        void fillWithScores(char algorithm_choice, std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q, const DistanceScoreMatrix& matrix);
    void fillWithScores(char algorithm_choice, std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q, const DistanceMatrix& matrix);
    
        DirectionMatrix getDirectionMatrix(char algorithm_choice);
    
        void algorithmSmithWaterman(std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q, DistanceScoreMatrix& matrix, int substructure_length, std::vector<std::string> p_aminoacid_sequence, std::vector<std::string> q_aminoacid_sequence, int alignment_representation_choice);
    
        void algorithmNeedlemanWunsch(DistanceScoreMatrix& matrix, int substructure_length, Protein p_protein, Protein q_protein, int alignment_representation_choice);
};

#endif
