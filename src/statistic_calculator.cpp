#include "statistic_calculator.h"

StatisticCalculator::StatisticCalculator(){
    
}

StatisticCalculator::StatisticCalculator(std::vector<double> data) : data_(data){
    
}

StatisticCalculator::StatisticCalculator(std::vector< std::vector<double> > data_2D) : data_2D_(data_2D){
    
}

void StatisticCalculator::setData(std::vector<double> data){
    data_ = data;
}

void StatisticCalculator::setData(std::vector< std::vector<double> > data_2D){
    
    std::vector<double> data;
    
    for(int i = 0; i < data_2D.size(); i++){
        for(int j = 0; j < data_2D[0].size(); j++){
            data.push_back(data_2D[i][j]);
        }
    }
    
    data_ = data;
}

std::vector<double> StatisticCalculator::getData() const{
    return data_;
}

double StatisticCalculator::getVariance(){
    double variance;
    double mean = this->getMean();
    double sum = 0;
    for (int i = 0; i < data_.size(); i++) {
        sum = sum + (data_[i] - mean)*(data_[i] - mean);
    }
    variance = sum / data_.size();
    return variance;
}

double StatisticCalculator::getMean(){
    double mean;
    double sum = 0;
    for(int i = 0; i < data_.size(); i++){
        sum += data_[i];
    }
    return sum / data_.size();
}

double StatisticCalculator::getStandardDeviation(){
    double standard_deviation = sqrt(this->getVariance());
    return standard_deviation;
}

