#ifndef _SEQUENCE_ALIGNER_
#define _SEQUENCE_ALIGNER_
#include <iostream>
#include <string>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <vector>
#include <utility>
#include "distance_matrix.h"
#include "distance_score_matrix.h"
#include "constants.hpp"
#define NORTH 360
#define NORTH_WEST 315
#define WEST 270

class SequenceAligner{
    private:
        // coordinates: (rows; columns)
        std::vector<int> directions_;
        std::vector< std::pair<int, int> > coordinates_;
        std::vector<int> subunit_chain_P_;
        std::vector<int> subunit_chain_Q_;
        double identity_score_;
        DistanceMatrix* score_matrix_;
    
        double gap_open_penalty_;
        double gap_ext_penalty_;
    public:
    SequenceAligner();
    SequenceAligner(std::vector<int> directions, std::vector< std::pair<int, int> > coordinates, std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q, double gap_open_penalty, double gap_ext_penalty);
        //SequenceAligner(std::vector<int> directions, std::vector< std::pair<int, int> > coordinates, std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q, DistanceScoreMatrix& score_matrix, double gap_open_penalty, double gap_ext_penalty);
    
    void setDirections(std::vector<int> directions);
    void setCoordinates(std::vector< std::pair<int, int> > coordinates);
    void setSubunitChainP(std::vector<int> subunit_chain_P);
    void setSubunitChainQ(std::vector<int> subunit_chain_Q);
    void setScoreMatrix(DistanceMatrix* score_matrix);
    void setGapOpenPenalty(double gap_open_penalty);
    void setGapExtPenalty(double gap_ext_penalty);
        
        void displayDirections() const;
        void displayCoordinates() const;
        void alignSequences();
    
        std::vector<std::string> getAlignedSequenceP();
        std::vector<std::string> getAlignedSequenceQ();
        std::vector<std::string> getAlignedSequences();
    
        std::vector<std::string> getAlignedSequenceP(std::vector<std::string> p_aminoacid_sequence);
        std::vector<std::string> getAlignedSequenceQ(std::vector<std::string> q_aminoacid_sequence);
    std::vector<std::string> getAlignedSequenceP(std::vector<std::string> p_aminoacid_sequence, double& score_normed_by_P);
    std::vector<std::string> getAlignedSequenceQ(std::vector<std::string> q_aminoacid_sequence, double& score_normed_by_Q);
    
        double getIdentity();
        void increaseIdentityScore(std::pair<int, int> coordinates, bool& is_gap_start);
        void decreaseIdentityScore(std::pair<int, int> coordinates, bool& is_gap_start);
    void normalizeIdentityScore(int chain_length);
};

#endif
