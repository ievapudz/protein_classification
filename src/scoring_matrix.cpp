#include "scoring_matrix.h"

ScoringMatrix::ScoringMatrix(int rows, int columns, double gap_open_penalty, double gap_ext_penalty) :
    rows_(rows),
    columns_(columns),
    gap_open_penalty_(gap_open_penalty),
    gap_ext_penalty_(gap_ext_penalty),
    scoring_matrix_(rows, std::vector<double>(columns, 0)){
    
}

void ScoringMatrix::setGapOpenPenalty(double gap_open_penalty){
    gap_open_penalty_ = gap_open_penalty;
}

void ScoringMatrix::setGapExtPenalty(double gap_ext_penalty){
    gap_ext_penalty_ = gap_ext_penalty;
}

int ScoringMatrix::getRows() const{
    return rows_;
}

int ScoringMatrix::getColumns() const{
    return columns_;
}

std::vector< std::vector<double> > ScoringMatrix::getScoringMatrix() const{
    return scoring_matrix_;
}

double ScoringMatrix::getGapOpenPenalty() const{
    return gap_open_penalty_;
}

double ScoringMatrix::getGapExtPenalty() const{
    return gap_ext_penalty_;
}

void ScoringMatrix::printScoringMatrix() const{
    int i, j;
    for (i = 0; i < rows_; i++){
        for (j = 0; j < columns_; j++){
            std::cout << std::setw(10) << std::right << scoring_matrix_[i][j];
        }
        std::cout << std::endl;
    }
}

void ScoringMatrix::fillWithGapPenalties(){
    // Method that fills "0th" row and "0th" column with gap extension penalty scores.
    int i, j;
    for(i = 0; i < rows_; i++){
        scoring_matrix_[i][0] = gap_ext_penalty_ * i * (-1);
    }
    for(j = 0; j < columns_; j++){
        scoring_matrix_[0][j] = gap_ext_penalty_ * j * (-1);
    }
}

std::pair<int, int> ScoringMatrix::getPair(int p_subunit, int q_subunit){
    // Method that returns a pair made of subunit numbers in P and Q chains.
    std::pair<int, int> pair = std::make_pair(p_subunit, q_subunit);
    return pair;
}

double ScoringMatrix::getNorthResult(double north_arg, double penalty){
    double north_result = north_arg - penalty;
    return north_result;
}

double ScoringMatrix::getNorthWestResult(double north_west_arg, double score){
    double north_west_result = north_west_arg + score;
    return north_west_result;
}

double ScoringMatrix::getWestResult(double west_arg, double penalty){
    double west_result = west_arg - penalty;
    return west_result;
}

double ScoringMatrix::getInitialMaxOper(double north_arg, double north_west_arg, double west_arg){
    double initial_max_oper = 0;
    
    if(north_arg < initial_max_oper)
       initial_max_oper = north_arg;
    if(north_west_arg < initial_max_oper)
        initial_max_oper = north_west_arg;
    if(west_arg < initial_max_oper)
        initial_max_oper = west_arg;
    
    return initial_max_oper;
}

double ScoringMatrix::getMaxOperationValue(double max_oper, double north, double north_west, double west){
    
    if(north > max_oper)
        max_oper = north;
    if(north_west > max_oper)
        max_oper = north_west;
    if(west > max_oper)
        max_oper = west;
    
    return max_oper;
}

double ScoringMatrix::getMaxOperationValue(double max_oper, double north, double north_west, double west, int& penalty_decision_matrix_element){
    
    if(north > max_oper){
        max_oper = north;
        penalty_decision_matrix_element = 1;
    }
    if(north_west > max_oper){
        max_oper = north_west;
    }
    if(west > max_oper){
        max_oper = west;
        penalty_decision_matrix_element = 2;
    }
    
    return max_oper;
}

void ScoringMatrix::fillWithScores(char algorithm_choice, std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q, const DistanceScoreMatrix& matrix, char initial){
    
    double north, north_west, west;
    int gap_ext_count = 0;
    double max_oper;
    
    for(int i = 1; i < rows_; i++){
        for(int j = 1; j < columns_; j++){
            // Pair holds numbers of subunits in aminoacid sequence.
            std::pair<int, int> pair = this->getPair(subunit_chain_P[i-1], subunit_chain_Q[j-1]);
            
            north_west = this->getNorthWestResult(scoring_matrix_[i-1][j-1], matrix.getScore(pair));
            
            if(pair.first != pair.second){
                gap_ext_count++;
                if(gap_ext_count == 1){
                    north = this->getNorthResult(scoring_matrix_[i-1][j], gap_open_penalty_);
                    west = this->getWestResult(scoring_matrix_[i][j-1], gap_open_penalty_);
                }
                if(gap_ext_count > 1){
                    north = this->getNorthResult(scoring_matrix_[i-1][j], gap_ext_penalty_);
                    west = this->getWestResult(scoring_matrix_[i][j-1], gap_ext_penalty_);
                }
            }
            else{
                gap_ext_count = 0;
                north = this->getNorthResult(scoring_matrix_[i-1][j], gap_open_penalty_);
                west = this->getWestResult(scoring_matrix_[i][j-1], gap_open_penalty_);
            }
            
            if(algorithm_choice == '1'){
                max_oper = 0;
            }
            else if(algorithm_choice == '2'){
                max_oper = this->getInitialMaxOper(scoring_matrix_[i][j-1],scoring_matrix_[i-1][j-1], scoring_matrix_[i-1][j]);
            }
            else{
                std::invalid_argument ia("Exception in ScoringMatrix::fillWithScores() invalid algorithm choice.");
                throw ia;
            }
            // Finding and setting the max value from operations.
            scoring_matrix_[i][j] = this->getMaxOperationValue(max_oper, north, north_west, west);
        }
    }
}


void ScoringMatrix::fillWithScores(char algorithm_choice, std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q, const DistanceScoreMatrix& matrix){
    
    double north, north_west, west;
    int gap_ext_count = 0;
    double max_oper;
    
    PenaltyDecisionMatrix penalty_decision_matrix(rows_, columns_);
    penalty_decision_matrix.fillMatrix();
    
    for(int i = 1; i < rows_; i++){
        for(int j = 1; j < columns_; j++){
            // Pair holds numbers of subunits in aminoacid sequence.
            std::pair<int, int> pair = this->getPair(subunit_chain_P[i-1], subunit_chain_Q[j-1]);
            
            north_west = this->getNorthWestResult(scoring_matrix_[i-1][j-1], matrix.getScore(pair));
            
            if(penalty_decision_matrix.getMatrix()[i-1][j] != 1){
                north = this->getNorthResult(scoring_matrix_[i-1][j], gap_open_penalty_);
            }else{
                north = this->getNorthResult(scoring_matrix_[i-1][j], gap_ext_penalty_);
            }
            
            if(penalty_decision_matrix.getMatrix()[i][j-1] != 2){
                west = this->getWestResult(scoring_matrix_[i][j-1], gap_open_penalty_);
            }else{
                west = this->getWestResult(scoring_matrix_[i][j-1], gap_ext_penalty_);
            }
            
            if(algorithm_choice == '1'){
                max_oper = 0;
            }
            else if(algorithm_choice == '2'){
                max_oper = this->getInitialMaxOper(scoring_matrix_[i][j-1], scoring_matrix_[i-1][j-1], scoring_matrix_[i-1][j]);
            }
            else{
                std::invalid_argument ia("Exception in ScoringMatrix::fillWithScores() invalid algorithm choice.");
                throw ia;
            }
            // Finding and setting the max value from operations.
            scoring_matrix_[i][j] = this->getMaxOperationValue(max_oper, north, north_west, west, penalty_decision_matrix.getMatrix()[i][j]);
            
        }
    }
}

DirectionMatrix ScoringMatrix::getDirectionMatrix(char algorithm_choice){
    
    DirectionMatrix direction_matrix(rows_, columns_);
    
    Traceback traceback(rows_, columns_, scoring_matrix_);
    traceback.setAlgorithmChoice(algorithm_choice);
    
    direction_matrix = traceback.getDirectionMatrix();
    
    return direction_matrix;
}

void ScoringMatrix::algorithmSmithWaterman(std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q, DistanceScoreMatrix& matrix, int substructure_length, std::vector<std::string> p_aminoacid_sequence, std::vector<std::string> q_aminoacid_sequence, int alignment_representation_choice){
    try{
        this->fillWithScores('1', subunit_chain_P, subunit_chain_Q, matrix);
        
        DirectionMatrix direction_matrix = this->getDirectionMatrix('1');
        
        SequenceAligner seq_al(direction_matrix.returnDirections(), direction_matrix.returnNonZeroCoords(), subunit_chain_P, subunit_chain_Q);

        switch(alignment_representation_choice){
            case 1:{
                std::string alignment_file_name = std::to_string(substructure_length);
                alignment_file_name.append("_alignment_file.txt");
                TXTFile alignment_file(alignment_file_name, seq_al.getAlignedSequences());
                alignment_file.writeData("./alignment_results/");
                break;
            }
            case 2:{
                std::string traditional_alignment_file_name = std::to_string(substructure_length);
                traditional_alignment_file_name.append("_traditional_alignment_file.txt");
                TXTFile traditional_alignment_file(traditional_alignment_file_name);
                traditional_alignment_file.writeData("./alignment_results_traditional/", seq_al.getAlignedSequenceP(p_aminoacid_sequence), seq_al.getAlignedSequenceQ(q_aminoacid_sequence));
                break;
            }
            default:{
                std::invalid_argument ia("Exception in ScoringMatrix::algorithmSmithWaterman() invalid algorithm choice.");
                throw ia;
            }
        }
    }catch(std::out_of_range& oor){
        std::cerr << "Out of range: " << oor.what() << std::endl;
    }
    catch(std::invalid_argument& ia){
        std::cerr << ia.what() << std::endl;
    }
}

void ScoringMatrix::algorithmNeedlemanWunsch(std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q, DistanceScoreMatrix& matrix, int substructure_length, std::vector<std::string> p_aminoacid_sequence, std::vector<std::string> q_aminoacid_sequence, int alignment_representation_choice){
    this->fillWithGapPenalties();
    try{
        this->fillWithScores('2', subunit_chain_P, subunit_chain_Q, matrix);
        
        DirectionMatrix direction_matrix = this->getDirectionMatrix('2');
        //direction_matrix.displayMatrix();
        
        SequenceAligner seq_al(direction_matrix.returnDirections(), direction_matrix.returnNonZeroCoords(), subunit_chain_P, subunit_chain_Q);
        
        switch(alignment_representation_choice){
            case 1:{
                std::string alignment_file_name = std::to_string(substructure_length);
                alignment_file_name.append("_alignment_file.txt");
                TXTFile alignment_file(alignment_file_name, seq_al.getAlignedSequences());
                alignment_file.writeData("./alignment_results/");
                break;
            }
            case 2:{
                std::string traditional_alignment_file_name = std::to_string(substructure_length);
                traditional_alignment_file_name.append("_traditional_alignment_file.txt");
                TXTFile traditional_alignment_file(traditional_alignment_file_name);
                traditional_alignment_file.writeData("./alignment_results_traditional/", seq_al.getAlignedSequenceP(p_aminoacid_sequence), seq_al.getAlignedSequenceQ(q_aminoacid_sequence));
                break;
            }
            default:{
                std::invalid_argument ia("Exception in ScoringMatrix::algorithmNeedlemanWunsch() invalid algorithm choice.");
                throw ia;
            }
        }
            
    }catch(std::out_of_range& oor){
        std::cerr << "Out of range: " << oor.what() << std::endl;
    }
    catch(std::invalid_argument& ia){
        std::cerr << ia.what() << std::endl;
    }
}






