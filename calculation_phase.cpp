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

PreparatoryPhase* CalculationPhase::getPreparatory(){
    return preparatory_;
}

int CalculationPhase::getSubstructureLength(){
    return substructure_length_;
}

std::pair< std::vector<std::string>, std::vector<std::string> >& CalculationPhase::getAlignment(){
    return alignment_;
}

std::pair<double, double> CalculationPhase::getIdentity(){
    return identity_;
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

void CalculationPhase::align(DirectionMatrix& direction_matrix, DistanceMatrix* distance_matrix){
    
    seq_al_.setDirections(direction_matrix.returnDirections());
    seq_al_.setCoordinates(direction_matrix.returnNonZeroCoords());
    seq_al_.setSubunitChainP(preparatory_->p_protein_.getSubunitChain());
    seq_al_.setSubunitChainQ(preparatory_->q_protein_.getSubunitChain());
    seq_al_.setScoreMatrix(distance_matrix);
    seq_al_.setGapOpenPenalty(preparatory_->constants_.gapOpenPenalty());
    seq_al_.setGapExtPenalty(preparatory_->constants_.gapExtPenalty());
    
    double identity_score_by_p = 0.0;
    double identity_score_by_q = 0.0;
   
    alignment_ = std::make_pair( seq_al_.getAlignedSequenceP(preparatory_->p_protein_.getSequence(), identity_score_by_p), seq_al_.getAlignedSequenceQ(preparatory_->q_protein_.getSequence(), identity_score_by_q));
    
    identity_ = std::make_pair(identity_score_by_p, identity_score_by_q);
}

void CalculationPhase::alignNumerally(DirectionMatrix& direction_matrix, DistanceMatrix* distance_matrix){
    seq_al_.setDirections(direction_matrix.returnDirections());
    seq_al_.setCoordinates(direction_matrix.returnNonZeroCoords());
    seq_al_.setSubunitChainP(preparatory_->p_protein_.getSubunitChain());
    seq_al_.setSubunitChainQ(preparatory_->q_protein_.getSubunitChain());
    seq_al_.setScoreMatrix(distance_matrix);
    seq_al_.setGapOpenPenalty(preparatory_->constants_.gapOpenPenalty());
    seq_al_.setGapExtPenalty(preparatory_->constants_.gapExtPenalty());
    
    double identity_score_by_p = 0.0;
    double identity_score_by_q = 0.0;
   
    alignment_ = std::make_pair( seq_al_.getAlignedSequencePNumerally(preparatory_->p_protein_.getSequence(), identity_score_by_p), seq_al_.getAlignedSequenceQNumerally(preparatory_->q_protein_.getSequence(), identity_score_by_q));
    
    identity_ = std::make_pair(identity_score_by_p, identity_score_by_q);
}

void CalculationPhase::run(){
    this->setMean();
    this->setStandardDeviation();
    
    DistanceMatrix distance_matrix = this->calculateDistanceScoreMatrix();

    preparatory_->p_protein_.setSubunitChain(block_distance_calc_.getNumberedSubunitChain(1));
    preparatory_->q_protein_.setSubunitChain(block_distance_calc_.getNumberedSubunitChain(2));
    
    DirectionMatrix direction_matrix = this->algorithmNeedlemanWunsch(distance_matrix);
    this->align(direction_matrix, &distance_matrix);
}

