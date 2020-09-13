#include "penalty_decision_matrix.h"

PenaltyDecisionMatrix::PenaltyDecisionMatrix(int rows, int columns) : rows_(rows), columns_(columns), matrix_(rows_, std::vector<int>(columns_, 0)){
    
}

void PenaltyDecisionMatrix::fillMatrix(){
    for(int i = 1; i < rows_; i++){
        matrix_[i][0] = 1;
    }
    for(int j = 1; j < columns_; j++){
        matrix_[0][j] = 2;
    }
}

void PenaltyDecisionMatrix::displayMatrix() const{
    int i, j;
    for (i = 0; i < rows_; i++){
        for (j = 0; j < columns_; j++){
            std::cout << std::setw(4) << std::right << matrix_[i][j];
        }
        std::cout << std::endl;
    }
}

std::vector< std::vector<int> > PenaltyDecisionMatrix::getMatrix(){
    return matrix_;
}


