#include "sequence_aligner.h"

//---SequenceAligner class

SequenceAligner::SequenceAligner(std::vector<int> directions, std::vector< std::pair<int, int> > coordinates, std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q) :
    directions_(directions),
    coordinates_(coordinates),
    subunit_chain_P_(subunit_chain_P),
    subunit_chain_Q_(subunit_chain_Q){
        
}

void SequenceAligner::displayDirections() const{
    for(int k = 0; k < directions_.size(); k++){
        std::cout << directions_[k] << std::endl;
    }
}

void SequenceAligner::displayCoordinates() const{
    for(int k = 0; k < coordinates_.size(); k++){
        std::cout << "(" << coordinates_[k].first << "; " << coordinates_[k].second << ")\n";
    }
}

void SequenceAligner::alignSequences(){
    std::vector<std::string> aligned_P = this->getAlignedSequenceP();
    std::vector<std::string> aligned_Q = this->getAlignedSequenceQ();
    
    int alignment_size = aligned_P.size();
    if(aligned_Q.size() < alignment_size)
        alignment_size = aligned_Q.size();
    
    for(int k = 0; k < alignment_size; k++){
        std::cout << std::setw(4) << std::right << aligned_P[k] << " " << aligned_Q[k] << std::endl;
    }
}

std::vector<std::string> SequenceAligner::getAlignedSequenceP(){
    std::vector<std::string> aligned_p;
    
    for(int k = 0; k < coordinates_.size(); k++){
        if(directions_[k] == -1){
            if(subunit_chain_Q_[coordinates_[k].second - 1] == subunit_chain_P_[coordinates_[k].first - 1]){
                aligned_p.push_back(std::to_string(subunit_chain_P_[coordinates_[k].first - 1]));
            }
            else aligned_p.push_back("-");
        }
        if(directions_[k] == NORTH){
            aligned_p.push_back(std::to_string(subunit_chain_P_[coordinates_[k].first - 1]));
        }
        if(directions_[k] == NORTH_WEST){
            aligned_p.push_back(std::to_string(subunit_chain_P_[coordinates_[k].first - 1]));
        }
        if(directions_[k] == WEST){
            aligned_p.push_back("-");
        }
    }
    return aligned_p;
}

std::vector<std::string> SequenceAligner::getAlignedSequenceQ(){
    std::vector<std::string> aligned_q;
    
    for(int k = 0; k < coordinates_.size(); k++){
        if(directions_[k] == -1){
            aligned_q.push_back(std::to_string(subunit_chain_Q_[coordinates_[k].second - 1]));
        }
        if(directions_[k] == NORTH){
            aligned_q.push_back("-");
        }
        if(directions_[k] == NORTH_WEST){
            aligned_q.push_back(std::to_string(subunit_chain_Q_[coordinates_[k].second - 1]));
        }
        if(directions_[k] == WEST){
            aligned_q.push_back(std::to_string(subunit_chain_Q_[coordinates_[k].second - 1]));
        }
    }
    
    return aligned_q;
}

std::vector<std::string> SequenceAligner::getAlignedSequences(){
    std::vector<std::string> aligned_sequences;
    
    std::vector<std::string> aligned_P = this->getAlignedSequenceP();
    std::vector<std::string> aligned_Q = this->getAlignedSequenceQ();
    
    int alignment_size = aligned_P.size();
    if(aligned_Q.size() < alignment_size)
        alignment_size = aligned_Q.size();
    
    for(int k = 0; k < alignment_size; k++){
        std::string alignment_subunit = aligned_P[k];
        alignment_subunit.append(" ");
        alignment_subunit.append(aligned_Q[k]);
        aligned_sequences.push_back(alignment_subunit);
    }
    
    return aligned_sequences;
}

