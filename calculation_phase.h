#ifndef _CALCULATION_PHASE_
#define _CALCULATION_PHASE_
#include "preparatory_phase.h"
#include "./src/block_distance.h"
#include "./src/statistic_calculator.h"
#include "./src/distance_matrix.h"
#include "./src/scoring_matrix.h"

class CalculationPhase{
    PreparatoryPhase* preparatory_;
    int substructure_length_;
    
    BlockDistanceCalculator block_distance_calc_;
    StatisticCalculator statistic_calc_;
    
    double mean_;
    double standard_deviation_;
    
public:
    CalculationPhase(PreparatoryPhase* preparatory, int substructure_length);
    void setMean();
    void setStandardDeviation();
    DistanceMatrix calculateDistanceMatrix();
    void run();
};

#endif

