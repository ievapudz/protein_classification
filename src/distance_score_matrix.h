#ifndef _DISTANCE_SCORE_MATRIX_
#define _DISTANCE_SCORE_MATRIX_
#include <iostream>
#include <vector>
#include <utility>
#include <stdexcept>
#include <iomanip>

class DistanceScoreMatrix{
    private:
        int rows_;
        int columns_;
        std::vector< std::vector<double> > matrix_;
    public:
        DistanceScoreMatrix(int rows, int columns);
        DistanceScoreMatrix(std::vector< std::vector<double> > matrix);
    
        void displayMatrix(int neighbour_number_left) const;
    
        void setMatrix(std::vector< std::vector<double> >& score_matrix);
        
        double getMaxValue();
        double getMinValue();
        std::vector< std::vector<double> > getMatrix() const;
        double getScore(std::pair<int, int> pair) const;
    
        void limitScores(int z_score_min, int z_score_max);
        void reverseScoreSigns();
};

#endif
