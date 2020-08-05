#include "distance_matrix.h"

DistanceMatrix::DistanceMatrix(std::vector< std::vector<double> > matrix) : matrix_(matrix){
    
    rows_ = matrix.size();
    columns_ = matrix[0].size();
    
    if(rows_ == 0){
        std::length_error le("Zero length of rows in distance matrix.");
        throw le;
    }
    if(columns_ == 0){
        std::length_error le("Zero length of columns in distance matrix.");
        throw le;
    }
    
}

void DistanceMatrix::displayMatrix(int neighbour_number_left) const{
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


void DistanceMatrix::setMatrix(std::vector< std::vector<double> > matrix){
    matrix_ = matrix;
}

std::vector< std::vector<double> > DistanceMatrix::getMatrix() const{
    return matrix_;
}

double DistanceMatrix::getElement(int index_row, int index_column) const{
    return matrix_[index_row][index_column];
}

double DistanceMatrix::getMaxValue(){
    
    double max_value = matrix_[0][0];
    
    for(int i = 0; i < rows_; i++){
        for(int j = 0; j < columns_; j++){
            if(matrix_[i][j] > max_value)
                max_value = matrix_[i][j];
        }
    }
    
    return max_value;
}

double DistanceMatrix::getMinValue(){
    
    double min_value = matrix_[0][0];
    
    for(int i = 0; i < rows_; i++){
        for(int j = 0; j < columns_; j++){
            if(matrix_[i][j] < min_value)
                min_value = matrix_[i][j];
        }
    }
    
    return min_value;
}

double DistanceMatrix::getMean(StatisticCalculator& statistic_calculator){
    
    double mean = statistic_calculator.getMean();
    
    return mean;
}

double DistanceMatrix::getStandardDeviation(StatisticCalculator& statistic_calculator){
    
    double standard_deviation = statistic_calculator.getStandardDeviation();
    
    return standard_deviation;
}

std::vector< std::vector<double> > DistanceMatrix::calculateScoreMatrix(double mean, double standard_deviation){
    /* Method that calculates z-scores for each distance in the distance matrix. Formula: https://en.wikipedia.org/wiki/Standard_score */
    
    std::vector< std::vector<double> > score_matrix(rows_, std::vector<double>(columns_, 0.0));
    
    for(int i = 0; i < rows_; i++){
        for(int j = 0; j < columns_; j++){
            double z_score = ( matrix_[i][j] - mean ) / standard_deviation;
            score_matrix[i][j] = z_score;
        }
    }
    
    return score_matrix;
}
