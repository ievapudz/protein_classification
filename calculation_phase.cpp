#include "calculation_phase.h"

CalculationPhase::CalculationPhase(PreparatoryPhase* preparatory, int substructure_length) : preparatory_(preparatory), substructure_length_(substructure_length), block_distance_calc_(preparatory->p_protein_, preparatory->q_protein_){

    preparatory_->setDistanceFile("./csv_distance_files/", substructure_length_);
    statistic_calc_.setData(preparatory_->distance_file_.getData(1));
    block_distance_calc_.setSubstructureLength(substructure_length_);
}

void CalculationPhase::setMean(){
    mean_ = statistic_calc_.getMean();
}

void CalculationPhase::setStandardDeviation(){
    standard_deviation_ = statistic_calc_.getStandardDeviation();
}

DistanceMatrix CalculationPhase::calculateDistanceMatrix(){
    DistanceMatrix dm(block_distance_calc_.getSuccessiveDistances(1));
    return dm;
}

void CalculationPhase::run(){
    this->setMean();
    std::cout << mean_ << " " << standard_deviation_ << std::endl;
    this->setStandardDeviation();
    // 1. get distance matrix (by calculating successive distances)
    DistanceMatrix dm = this->calculateDistanceMatrix();
    
    // 2. modify distance matrix to distance score matrix (using mean and standard deviation, limiting to z_score_min and z_score_max, reversing signs)
    // 3. numbering chains according to substructure_length
    // 4. creating matrix for alignment (scoring_matrix)
    // 5. algorithm for alignment
}

