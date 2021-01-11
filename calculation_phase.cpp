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

std::vector<std::string>& CalculationPhase::getAlignment2(){
    return alignment_2_;
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

std::vector<std::string> CalculationPhase::algorithmNeedlemanWunsch(DistanceMatrix& match){

    std::vector< std::pair<int, int> > alignment;
    std::vector< std::string > final_alignment;
    
    std::vector<int> seq1 = preparatory_->p_protein_.getSubunitChain();
    std::vector<int> seq2 = preparatory_->q_protein_.getSubunitChain();

    if((seq1.size() != 0)&&(seq2.size() != 0)){
        double scores_matrix [seq1.size()+1][seq2.size()+1];
        int directions_matrix [seq1.size()+1][seq2.size()+1];

        for(int i = 1; i<=seq1.size(); i++){
            scores_matrix[i][0] = ( (i==1) ? (preparatory_->constants_.gapOpenPenalty() * (-1)) : (scores_matrix[i-1][0]-preparatory_->constants_.gapExtPenalty()));
            directions_matrix[i][0] = 1;
        }
        for(int j = 1; j<=seq2.size(); j++){
            scores_matrix[0][j] = ( (j==1) ? (preparatory_->constants_.gapOpenPenalty() * (-1)) : (scores_matrix[0][j-1]-preparatory_->constants_.gapExtPenalty()));
            directions_matrix[0][j] = 2;
        }
        
        std::pair<int, int> result_score_pos = std::make_pair(0, 0);
        
        for(int i = 1; i<=seq1.size(); i++){
            for(int j = 1; j<=seq2.size(); j++){
                int v1 = seq1[i-1];
                int v2 = seq2[j-1];
                
                std::pair<int, int> pair = std::make_pair(v1, v2);

                double match_score = scores_matrix[i-1][j-1] + match.getScore(pair) + ( ((i==1 && j==1) || (i==seq1.size() && j==seq2.size())) ? -1 : 0 );
                double deletion_score = scores_matrix[i-1][j] + ( directions_matrix[i-1][j] != 1 ? preparatory_->constants_.gapOpenPenalty() * (-1) : preparatory_->constants_.gapExtPenalty() * (-1));
                double insertion_score = scores_matrix[i][j-1] + ( directions_matrix[i][j-1] != 2 ? preparatory_->constants_.gapOpenPenalty() * (-1) : preparatory_->constants_.gapExtPenalty() * (-1));
                double max_score = std::max( match_score, std::max(deletion_score, insertion_score));
                directions_matrix[i][j] = ( max_score == insertion_score ? 2 : ( max_score == deletion_score ? 1 : 0));
                scores_matrix[i][j] = max_score;
            }
        }

        std::pair<int, int> result_score_pos_buffer = std::make_pair(seq1.size(), seq2.size());
        result_score_pos = result_score_pos_buffer;

      // constructing alignment

        int i = seq1.size();
        int j = seq2.size();
        
        double identity_score = 0.0;
        bool is_gap = false;

        while(i>0 && j>0){
            int dir = directions_matrix[i][j];
            if(dir == 0){
                i--;
                j--;
                is_gap = false;
                identity_score += match.getScore(std::make_pair(seq1[i], seq2[j]));
                std::pair<int, int> pair = std::make_pair(i, j);
                alignment.push_back(pair);
            }else if(dir == 1){
                i--;
                if(is_gap){
                  identity_score -= preparatory_->constants_.gapExtPenalty();
                }else{
                  is_gap = true;
                  identity_score -= preparatory_->constants_.gapOpenPenalty();
                }
                std::pair<int, int> pair = std::make_pair(i, -1);
                alignment.push_back(pair);
            }else{
                j--;
                if(is_gap){
                  identity_score -= preparatory_->constants_.gapExtPenalty();
                }else{
                  is_gap = true;
                  identity_score -= preparatory_->constants_.gapOpenPenalty();
                }
                std::pair<int, int> pair = std::make_pair(-1, j);
                alignment.push_back(pair);
            }
        }
        
        identity_ = std::make_pair(identity_score / seq1.size(), identity_score / seq2.size());
        
        alignment.pop_back();
        std::reverse(alignment.begin(), alignment.end());

        for(int k = 0; k < alignment.size(); k++){
            if((alignment[k].first >= 0)&&(alignment[k].second >= 0)){
                std::string alignment_cell = std::to_string(seq1[alignment[k].first]);
                alignment_cell.append(" ");
                alignment_cell.append(std::to_string(seq2[alignment[k].second]));
                final_alignment.push_back(alignment_cell);
            }else if((alignment[k].first < 0)&&(alignment[k].second > 0)){
                std::string alignment_cell = "- ";
                alignment_cell.append( std::to_string(seq2[alignment[k].second]));
                final_alignment.push_back(alignment_cell);
            }else{
                std::string alignment_cell = std::to_string( seq1[alignment[k].first] );
                alignment_cell.append(" -");
                final_alignment.push_back(alignment_cell);
            }
        }
    }
    return final_alignment;
  }

void CalculationPhase::run(){
    this->setMean();
    this->setStandardDeviation();
    
    DistanceMatrix distance_matrix = this->calculateDistanceScoreMatrix();

    preparatory_->p_protein_.setSubunitChain(block_distance_calc_.getNumberedSubunitChain(1));
    preparatory_->q_protein_.setSubunitChain(block_distance_calc_.getNumberedSubunitChain(2));
    
    alignment_2_ = this->algorithmNeedlemanWunsch(distance_matrix);

}

