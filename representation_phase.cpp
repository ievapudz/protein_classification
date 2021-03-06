#include "representation_phase.h"

RepresentationPhase::RepresentationPhase(CalculationPhase* calculations) : calculations_(calculations), p_protein_(calculations_->getPreparatory()->p_protein_),
    q_protein_(calculations_->getPreparatory()->q_protein_){
    
}

void RepresentationPhase::representAlignment(){
    std::string alignment_file_name = std::to_string(calculations_->getSubstructureLength());
    alignment_file_name.append("_"+p_protein_.getName()+"_"+q_protein_.getName()+"_traditional_alignment_file.txt");
    TXTFile alignment_file(alignment_file_name);
    alignment_file.writeData("./alignment_results_traditional/", calculations_->getAlignment().first, calculations_->getAlignment().second);
}

void RepresentationPhase::representNumeralAlignment(){
    std::string alignment_file_name = std::to_string(calculations_->getSubstructureLength());
    alignment_file_name.append("_"+p_protein_.getName()+"_"+q_protein_.getName()+"_numeral_alignment_file.txt");
    
    TXTFile alignment_file(alignment_file_name, calculations_->getAlignment2());
    
    alignment_file.writeData("./alignment_results_numeral/");
}

std::pair<double, double> RepresentationPhase::representIdentityScore(){
    return calculations_->getIdentity();
}
