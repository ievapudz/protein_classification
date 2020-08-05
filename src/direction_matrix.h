#ifndef _DIRECTION_MATRIX_
#define _DIRECTION_MATRIX_
#include <iostream>
#include <iomanip>
#include <vector>
#include <utility>

class DirectionMatrix{
    private:
        int rows_;
        int columns_;
        std::vector< std::vector<int> > direction_matrix_;
    public:
        DirectionMatrix(int rows, int columns);
        void setDirection(int rows_coordinate, int column_coordinate, int direction);
    
        std::vector< std::vector<int> > returnMatrix() const;
        void displayMatrix() const;
        std::vector< std::pair<int, int> > returnNonZeroCoords() const;
        std::vector<int> returnDirections() const;
};

#endif
