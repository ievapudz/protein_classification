#ifndef _STATISTIC_CALCULATOR_
#define _STATISTIC_CALCULATOR_
#include <vector>
#include <iostream>
#include <math.h>

class StatisticCalculator{
        std::vector<double> data_;
        std::vector< std::vector<double> > data_2D_;
    public:
        StatisticCalculator();
        StatisticCalculator(std::vector<double> data);
        StatisticCalculator(std::vector< std::vector<double> > data_2D);
        void setData(std::vector<double> data);
        void setData(std::vector< std::vector<double> > data_2D);
    
        std::vector<double> getData() const;
        double getVariance();
        double getMean();
        double getStandardDeviation();
};

#endif
