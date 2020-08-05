#include "direction_matrix.h"

DirectionMatrix::DirectionMatrix(int rows, int columns) :
    rows_(rows),
    columns_(columns),
    direction_matrix_(rows, std::vector<int>(columns, 0)){
    
}

void DirectionMatrix::setDirection(int row_coordinate, int column_coordinate, int direction){
    direction_matrix_[row_coordinate][column_coordinate] = direction;
}

std::vector< std::vector<int> > DirectionMatrix::returnMatrix() const{
    return this->direction_matrix_;
}

void DirectionMatrix::displayMatrix() const{
    int i, j;
    for (i = 0; i < rows_; i++){
        for (j = 0; j < columns_; j++){
            std::cout << std::setw(4) << std::right << direction_matrix_[i][j];
        }
        std::cout << std::endl;
    }
}

std::vector< std::pair<int, int> > DirectionMatrix::returnNonZeroCoords() const{
    int i, j;
    std::vector< std::pair<int, int> > non_zero_coords;
    for (i = 0; i < rows_; i++){
        for (j = 0; j < columns_; j++){
            if(direction_matrix_[i][j] != 0){
                std::pair<int, int> coordinates;
                coordinates = std::make_pair (i, j);
                non_zero_coords.push_back(coordinates);
            }
        }
    }
    return non_zero_coords;
}

std::vector<int> DirectionMatrix::returnDirections() const{
    int i, j;
    std::vector<int> directions;
    for (i = 0; i < rows_; i++){
        for (j = 0; j < columns_; j++){
            if(direction_matrix_[i][j] != 0){
                directions.push_back(direction_matrix_[i][j]);
            }
        }
    }
    return directions;
}

