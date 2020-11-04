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
    TXTFile alignment_file(alignment_file_name);
    alignment_file.writeData("./alignment_results_numeral/", calculations_->getAlignment().first, calculations_->getAlignment().second);
}

void RepresentationPhase::representIdentityScore(){
    std::vector<std::string> identity_figure;
    identity_figure.push_back(std::to_string(calculations_->getIdentity().first));
    identity_figure.push_back(std::to_string(calculations_->getIdentity().second));
    std::string identity_score_file_name = std::to_string(calculations_->getSubstructureLength());
    identity_score_file_name.append("_"+p_protein_.getName()+"_"+q_protein_.getName()+"_identity_score.txt");
    TXTFile identity_score_file(identity_score_file_name, identity_figure);
    identity_score_file.writeData("./alignment_results_traditional/");
}
