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

DistanceMatrix CalculationPhase::calculateDistanceScoreMatrix(){
    DistanceMatrix dm(block_distance_calc_.getSuccessiveDistances(1));
    // Calculating scores for distances in the matrix, limiting them to the interval and reversing the signs.
    dm.modifyToDistanceScoreMatrix(mean_, standard_deviation_, preparatory_->constants_.zScoreMin(), preparatory_->constants_.zScoreMax());
    return dm;
}

DirectionMatrix CalculationPhase::algorithmNeedlemanWunsch(DistanceMatrix& matrix){
    ScoringMatrix sm(preparatory_->p_protein_, preparatory_->q_protein_, preparatory_->constants_.gapOpenPenalty(), preparatory_->constants_.gapExtPenalty());
    sm.fillWithGapPenalties();
    sm.fillWithScores('2', preparatory_->p_protein_.getSubunitChain(), preparatory_->q_protein_.getSubunitChain(), matrix);
    return sm.getDirectionMatrix('2');
}

void CalculationPhase::align(DirectionMatrix& matrix){
    
    SequenceAligner seq_al_(matrix.returnDirections(), matrix.returnNonZeroCoords(), preparatory_->p_protein_.getSubunitChain(), preparatory_->q_protein_.getSubunitChain(), preparatory_->constants_.gapOpenPenalty(), preparatory_->constants_.gapExtPenalty());
   
    alignment_ = std::make_pair( seq_al_.getAlignedSequenceP(preparatory_->p_protein_.getSequence()), seq_al_.getAlignedSequenceQ(preparatory_->q_protein_.getSequence()));
    
    for(int i = 0; i < alignment_.first.size(); i++){
        std::cout << alignment_.first[i];
    }
    std::cout << std::endl;
    
    for(int i = 0; i < alignment_.second.size(); i++){
        std::cout << alignment_.second[i];
    }
}

void CalculationPhase::evaluateIdentity(DistanceMatrix& distance_matrix, DirectionMatrix& direction_matrix){
    
   
}

void CalculationPhase::run(){
    this->setMean();
    this->setStandardDeviation();
    
    DistanceMatrix distance_matrix = this->calculateDistanceScoreMatrix();

    preparatory_->p_protein_.setSubunitChain(block_distance_calc_.getNumberedSubunitChain(1));
    preparatory_->q_protein_.setSubunitChain(block_distance_calc_.getNumberedSubunitChain(2));
    
    DirectionMatrix direction_matrix = this->algorithmNeedlemanWunsch(distance_matrix);
    this->align(direction_matrix);
}

