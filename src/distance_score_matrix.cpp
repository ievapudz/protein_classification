#include "distance_score_matrix.h"

DistanceScoreMatrix::DistanceScoreMatrix(int rows, int columns) : rows_(rows), columns_(columns),
    matrix_(rows, std::vector<double>(columns, 0.0)){
        
    if(rows_ == 0){
        std::length_error le("Zero length of rows in distance score matrix.");
        throw le;
    }
    if(columns_ == 0){
        std::length_error le("Zero length of columns in distance score matrix.");
        throw le;
    }
    
}

DistanceScoreMatrix::DistanceScoreMatrix(std::vector< std::vector<double> > matrix) : matrix_(matrix){
    
    rows_ = matrix.size();
    columns_ = matrix[0].size();
    
    if(rows_ == 0){
        std::length_error le("Zero length of rows in distance score matrix.");
        throw le;
    }
    if(columns_ == 0){
        std::length_error le("Zero length of columns in distance score matrix.");
        throw le;
    }
}

void DistanceScoreMatrix::displayMatrix(int neighbour_number_left) const{
    //Method that displays the matrix of distances between chain subunits (numbers represent accurate subunit number in the aminoacid chain).
    
    for(int h = -1; h < columns_; h++){
        if(h == -1)
            std::cout << std::setw(10) << std::right << " ";
        else
            std::cout << std::setw(10) << std::right << h + neighbour_number_left + 1;
    }
    
    std::cout << std::endl;
    
    for(int i = 0; i < rows_; i++){
        std::cout << std::setw(10) << std::right << i + neighbour_number_left + 1;
        for(int j = 0; j < columns_; j++){
            std::cout << std::setw(10) << std::right << matrix_[i][j];
        }
        std::cout << std::endl;
    }
}

void DistanceScoreMatrix::setMatrix(std::vector< std::vector<double> >& score_matrix){
    matrix_ = score_matrix;
}

double DistanceScoreMatrix::getMaxValue(){
    
    double max_value = matrix_[0][0];
    
    for(int i = 0; i < rows_; i++){
        for(int j = 0; j < columns_; j++){
            if(matrix_[i][j] > max_value)
                max_value = matrix_[i][j];
        }
    }
    
    return max_value;
}

double DistanceScoreMatrix::getMinValue(){
    
    double min_value = matrix_[0][0];
    
    for(int i = 0; i < rows_; i++){
        for(int j = 0; j < columns_; j++){
            if(matrix_[i][j] < min_value)
                min_value = matrix_[i][j];
        }
    }
    
    return min_value;
}

const std::vector< std::vector<double> >& DistanceScoreMatrix::getMatrix() const{
    return matrix_;
}

double DistanceScoreMatrix::getScore(std::pair<int, int> pair) const{
    // Method that returns score of the corresponding pair of subunits.
    
    double score;
    
    for(int i = 0; i < rows_; i++){
        for(int j = 0; j < columns_; j++){
            if((i+1 == pair.first) && (j+1 == pair.second)){
                score = matrix_[i][j];
                break;
            }
        }
    }
    return score;
}

void DistanceScoreMatrix::limitScores(int z_score_min, int z_score_max){
    for(int i = 0; i < rows_; i++){
        for(int j = 0; j < columns_; j++){
            if(matrix_[i][j] < z_score_min)
                matrix_[i][j] = z_score_min;
            if(matrix_[i][j] > z_score_max)
                matrix_[i][j] = z_score_max;
        }
    }
}

void DistanceScoreMatrix::reverseScoreSigns(){
    for(int i = 0; i < rows_; i++){
        for(int j = 0; j < columns_; j++){
            matrix_[i][j] = matrix_[i][j] * (-1);
        }
    }
}

