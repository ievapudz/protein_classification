#ifndef _DISTANCE_MATRIX_
#define _DISTANCE_MATRIX_
#include <iostream>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <utility>
#include "statistic_calculator.h"

class DistanceMatrix{
    private:
        int rows_;
        int columns_;
        std::vector< std::vector<double> > matrix_;
    public:
        DistanceMatrix(std::vector< std::vector<double> > matrix);
    
        void displayMatrix(int neighbour_number_left) const;
    
        void setMatrix(std::vector< std::vector<double> > matrix);
    
        const std::vector< std::vector<double> >& getMatrix() const;
        double getElement(int index_row, int index_column) const;
        double getMaxValue();
        double getMinValue();
    
        double getMean(StatisticCalculator& statistic_calculator);
        double getStandardDeviation(StatisticCalculator& statistic_calculator);
    
        std::vector< std::vector<double> > calculateScoreMatrix(double mean, double standard_deviation);
    void modifyToDistanceScoreMatrix(double mean, double standard_deviation, double z_score_min, double z_score_max);
    double getScore(std::pair<int, int> pair) const;
};

#endif
