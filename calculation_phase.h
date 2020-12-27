#ifndef _CALCULATION_PHASE_
#define _CALCULATION_PHASE_
#include <utility>
#include "preparatory_phase.h"
#include "./src/block_distance.h"
#include "./src/statistic_calculator.h"
#include "./src/distance_matrix.h"
#include "./src/scoring_matrix.h"
#include "./src/direction_matrix.h"
#include "./src/sequence_aligner.h"

class CalculationPhase{
        PreparatoryPhase* preparatory_;
        int substructure_length_;
        
        BlockDistanceCalculator block_distance_calc_;
        StatisticCalculator statistic_calc_;
        SequenceAligner seq_al_;
        
        double mean_;
        double standard_deviation_;
    
        std::vector< std::string > alignment_2_;
        std::pair< std::vector<std::string>, std::vector<std::string> > alignment_;
        std::pair<double, double> identity_;
        
    public:
        CalculationPhase(PreparatoryPhase* preparatory, int substructure_length);
        void setMean();
        void setStandardDeviation();
    
        PreparatoryPhase* getPreparatory();
        int getSubstructureLength();
        std::pair< std::vector<std::string>, std::vector<std::string> >& getAlignment();
        std::pair<double, double> getIdentity();
    
        std::vector<std::string>& getAlignment2();
    
        DistanceMatrix calculateDistanceScoreMatrix();
        
        std::vector<std::string> algorithmNeedlemanWunsch(DistanceMatrix& match);
        void align(DirectionMatrix& direction_matrix, DistanceMatrix* distance_matrix);
    void alignNumerally(DirectionMatrix& direction_matrix, DistanceMatrix* distance_matrix);
        void run();
};

#endif

