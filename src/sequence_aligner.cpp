#include "sequence_aligner.h"

//---SequenceAligner class
SequenceAligner::SequenceAligner() : identity_score_(0.0){
    
}

SequenceAligner::SequenceAligner(std::vector<int> directions, std::vector< std::pair<int, int> > coordinates, std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q, double gap_open_penalty, double gap_ext_penalty) :
    directions_(directions),
    coordinates_(coordinates),
    subunit_chain_P_(subunit_chain_P),
    subunit_chain_Q_(subunit_chain_Q), identity_score_(0.0),
    gap_open_penalty_(gap_open_penalty),
    gap_ext_penalty_(gap_ext_penalty){
    
}
/*
SequenceAligner::SequenceAligner(std::vector<int> directions, std::vector< std::pair<int, int> > coordinates, std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q, DistanceScoreMatrix& score_matrix, double gap_open_penalty, double gap_ext_penalty) :
    directions_(directions),
    coordinates_(coordinates),
    subunit_chain_P_(subunit_chain_P),
    subunit_chain_Q_(subunit_chain_Q), identity_score_(0.0), score_matrix_(score_matrix),
    gap_open_penalty_(gap_open_penalty),
    gap_ext_penalty_(gap_ext_penalty){
        
}*/

void SequenceAligner::setDirections(std::vector<int> directions){
    directions_ = directions;
}

void SequenceAligner::setCoordinates(std::vector< std::pair<int, int> > coordinates){
    coordinates_ = coordinates;
}

void SequenceAligner::setSubunitChainP(std::vector<int> subunit_chain_P){
    subunit_chain_P_ = subunit_chain_P;
}

void SequenceAligner::setSubunitChainQ(std::vector<int> subunit_chain_Q){
    subunit_chain_Q_ = subunit_chain_Q;
}

void SequenceAligner::setScoreMatrix(DistanceMatrix* score_matrix){
    score_matrix_ = score_matrix;
}

void SequenceAligner::setGapOpenPenalty(double gap_open_penalty){
    gap_open_penalty_ = gap_open_penalty;
}

void SequenceAligner::setGapExtPenalty(double gap_ext_penalty){
    gap_ext_penalty_ = gap_ext_penalty;
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
        
        if(directions_[k] == 0){
            if(subunit_chain_Q_[coordinates_[k].second - 1] == subunit_chain_P_[coordinates_[k].first - 1]){
                aligned_p.push_back(std::to_string(subunit_chain_P_[coordinates_[k].first - 1]));
            }
            else aligned_p.push_back("-");
        }
        else if(directions_[k] == NORTH){
            aligned_p.push_back(std::to_string(subunit_chain_P_[coordinates_[k].first - 1]));
        }
        else if(directions_[k] == NORTH_WEST){
            aligned_p.push_back(std::to_string(subunit_chain_P_[coordinates_[k].first - 1]));
        }
        else if(directions_[k] == WEST){
            aligned_p.push_back("-");
        }
    }
    return aligned_p;
}

std::vector<std::string> SequenceAligner::getAlignedSequenceQ(){
    std::vector<std::string> aligned_q;
    
    for(int k = 0; k < coordinates_.size(); k++){
        
        if(directions_[k] == 0){
            aligned_q.push_back(std::to_string(subunit_chain_Q_[coordinates_[k].second - 1]));
        }
        else if(directions_[k] == NORTH){
            aligned_q.push_back("-");
        }
        else if(directions_[k] == NORTH_WEST){
            aligned_q.push_back(std::to_string(subunit_chain_Q_[coordinates_[k].second - 1]));
        }
        else if(directions_[k] == WEST){
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

std::vector<std::string> SequenceAligner::getAlignedSequenceP(std::vector<std::string> p_aminoacid_sequence){
    std::vector<std::string> aligned_p;
    
    bool is_gap_start = false;
    for(int k = 0; k < coordinates_.size(); k++){
        if(directions_[k] == 0){
            if(subunit_chain_Q_[coordinates_[k].second - 1] == subunit_chain_P_[coordinates_[k].first - 1]){
                aligned_p.push_back(p_aminoacid_sequence[ subunit_chain_P_[coordinates_[k].first - 1] - 1 ]);
            }
            else{
                aligned_p.push_back("-");
            }
        }
        else if(directions_[k] == NORTH){
            aligned_p.push_back(p_aminoacid_sequence[ subunit_chain_P_[coordinates_[k].first - 1] - 1 ]);
        }
        else if(directions_[k] == NORTH_WEST){
            aligned_p.push_back(p_aminoacid_sequence[ subunit_chain_P_[coordinates_[k].first - 1] - 1 ]);
        }
        else if(directions_[k] == WEST){
            aligned_p.push_back("-");
        }
    }
    return aligned_p;
}

std::vector<std::string> SequenceAligner::getAlignedSequenceQ(std::vector<std::string> q_aminoacid_sequence){
    
    std::vector<std::string> aligned_q;
    
    bool is_gap_start = false;
    for(int k = 0; k < coordinates_.size(); k++){
        if(directions_[k] == 0){
            aligned_q.push_back(q_aminoacid_sequence[ subunit_chain_Q_[coordinates_[k].second - 1] - 1 ]);
        }
        else if(directions_[k] == NORTH){
            aligned_q.push_back("-");
        }
        else if(directions_[k] == NORTH_WEST){
            aligned_q.push_back(q_aminoacid_sequence[ subunit_chain_Q_[coordinates_[k].second - 1] - 1 ]);
        }
        else if(directions_[k] == WEST){
            aligned_q.push_back(q_aminoacid_sequence[ subunit_chain_Q_[coordinates_[k].second - 1] - 1 ]);
        }
    }
    return aligned_q;
}

std::vector<std::string> SequenceAligner::getAlignedSequenceP(std::vector<std::string> p_aminoacid_sequence, double& score_normed_by_P){
   
    std::vector<std::string> aligned_p;
    identity_score_ = 0.0;
    
    bool is_gap_start = false;
    for(int k = 0; k < coordinates_.size(); k++){
        if(directions_[k] == 0){
            if(subunit_chain_Q_[coordinates_[k].second - 1] == subunit_chain_P_[coordinates_[k].first - 1]){
                aligned_p.push_back(p_aminoacid_sequence[ subunit_chain_P_[coordinates_[k].first - 1] - 1 ]);
                this->increaseIdentityScore(coordinates_[k], is_gap_start);
            }
            else{
                aligned_p.push_back("-");
                this->decreaseIdentityScore(coordinates_[k], is_gap_start);
            }
        }
        else if(directions_[k] == NORTH){
            aligned_p.push_back(p_aminoacid_sequence[ subunit_chain_P_[coordinates_[k].first - 1] - 1 ]);
            this->decreaseIdentityScore(coordinates_[k], is_gap_start);
        }
        else if(directions_[k] == NORTH_WEST){
            aligned_p.push_back(p_aminoacid_sequence[ subunit_chain_P_[coordinates_[k].first - 1] - 1 ]);
            this->increaseIdentityScore(coordinates_[k], is_gap_start);
        }
        else if(directions_[k] == WEST){
            aligned_p.push_back("-");
            this->decreaseIdentityScore(coordinates_[k], is_gap_start);
        }
    }
   
    this->normalizeIdentityScore(subunit_chain_P_.size());
    
    score_normed_by_P = identity_score_;
    
    return aligned_p;
}

std::vector<std::string> SequenceAligner::getAlignedSequenceQ(std::vector<std::string> q_aminoacid_sequence, double& score_normed_by_Q){
    
    std::vector<std::string> aligned_q;
    identity_score_ = 0.0;
    
    bool is_gap_start = false;
    for(int k = 0; k < coordinates_.size(); k++){
        if(directions_[k] == 0){
            aligned_q.push_back(q_aminoacid_sequence[ subunit_chain_Q_[coordinates_[k].second - 1] - 1 ]);
        }
        else if(directions_[k] == NORTH){
            aligned_q.push_back("-");
            this->decreaseIdentityScore(coordinates_[k], is_gap_start);
        }
        else if(directions_[k] == NORTH_WEST){
            aligned_q.push_back(q_aminoacid_sequence[ subunit_chain_Q_[coordinates_[k].second - 1] - 1 ]);
            this->increaseIdentityScore(coordinates_[k], is_gap_start);
        }
        else if(directions_[k] == WEST){
            aligned_q.push_back(q_aminoacid_sequence[ subunit_chain_Q_[coordinates_[k].second - 1] - 1 ]);
            this->decreaseIdentityScore(coordinates_[k], is_gap_start);
        }
    }
    this->normalizeIdentityScore(subunit_chain_Q_.size());
    
    score_normed_by_Q = identity_score_;
    
    return aligned_q;
}

std::vector<std::string> SequenceAligner::getAlignedSequencePNumeral(std::vector<std::string> p_aminoacid_sequence, double& score_normed_by_P){
    
    std::vector<std::string> aligned_p;
    identity_score_ = 0.0;
    
    bool is_gap_start = false;
    for(int k = 0; k < coordinates_.size(); k++){
        if(directions_[k] == 0){
            if(subunit_chain_Q_[coordinates_[k].second - 1] == subunit_chain_P_[coordinates_[k].first - 1]){
                aligned_p.push_back(std::to_string(subunit_chain_P_[coordinates_[k].first - 1]));
                this->increaseIdentityScore(coordinates_[k], is_gap_start);
            }
            else{
                aligned_p.push_back("-");
                this->decreaseIdentityScore(coordinates_[k], is_gap_start);
            }
        }
        else if(directions_[k] == NORTH){
            aligned_p.push_back(std::to_string(subunit_chain_P_[coordinates_[k].first - 1]));
            this->decreaseIdentityScore(coordinates_[k], is_gap_start);
        }
        else if(directions_[k] == NORTH_WEST){
            aligned_p.push_back(std::to_string(subunit_chain_P_[coordinates_[k].first - 1]));
            this->increaseIdentityScore(coordinates_[k], is_gap_start);
        }
        else if(directions_[k] == WEST){
            aligned_p.push_back("-");
            this->decreaseIdentityScore(coordinates_[k], is_gap_start);
        }
    }
   
    this->normalizeIdentityScore(subunit_chain_P_.size());
    
    score_normed_by_P = identity_score_;
    
    return aligned_p;
    
}

std::vector<std::string> SequenceAligner::getAlignedSequenceQNumeral(std::vector<std::string> q_aminoacid_sequence, double& score_normed_by_Q){
    
    std::vector<std::string> aligned_q;
    identity_score_ = 0.0;
    
    bool is_gap_start = false;
    for(int k = 0; k < coordinates_.size(); k++){
        if(directions_[k] == 0){
            aligned_q.push_back(std::to_string(subunit_chain_Q_[coordinates_[k].second - 1]));
        }
        else if(directions_[k] == NORTH){
            aligned_q.push_back("-");
            this->decreaseIdentityScore(coordinates_[k], is_gap_start);
        }
        else if(directions_[k] == NORTH_WEST){
            aligned_q.push_back(std::to_string(subunit_chain_Q_[coordinates_[k].second - 1]));
            this->increaseIdentityScore(coordinates_[k], is_gap_start);
        }
        else if(directions_[k] == WEST){
            aligned_q.push_back(std::to_string(subunit_chain_Q_[coordinates_[k].second - 1]));
            this->decreaseIdentityScore(coordinates_[k], is_gap_start);
        }
    }
    this->normalizeIdentityScore(subunit_chain_Q_.size());
    
    score_normed_by_Q = identity_score_;
    
    return aligned_q;
    
}

double SequenceAligner::getIdentity(){
    return identity_score_;
}

void SequenceAligner::increaseIdentityScore(std::pair<int, int> coordinates, bool& is_gap_start){
    is_gap_start = false;
    identity_score_ += this->score_matrix_->getScore(coordinates);
}

void SequenceAligner::decreaseIdentityScore(std::pair<int, int> coordinates, bool& is_gap_start){
    if(!is_gap_start){
        is_gap_start = true;
        identity_score_ -= gap_open_penalty_;
    }else{
        identity_score_ -= gap_ext_penalty_;
    }
}

void SequenceAligner::normalizeIdentityScore(int chain_length){
    identity_score_ /= chain_length;
}

