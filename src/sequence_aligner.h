#ifndef _SEQUENCE_ALIGNER_
#define _SEQUENCE_ALIGNER_
#include <iostream>
#include <string>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <vector>
#include <utility>
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
    public:
        SequenceAligner(std::vector<int> directions, std::vector< std::pair<int, int> > coordinates, std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q);
        
        void displayDirections() const;
        void displayCoordinates() const;
        void alignSequences();
    
        std::vector<std::string> getAlignedSequenceP();
        std::vector<std::string> getAlignedSequenceQ();
        std::vector<std::string> getAlignedSequences();
    
        std::vector<std::string> getAlignedSequenceP(std::vector<std::string> p_aminoacid_sequence);
        std::vector<std::string> getAlignedSequenceQ(std::vector<std::string> q_aminoacid_sequence);
    
        double getIdentity();
};

#endif
