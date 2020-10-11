#ifndef _CONSTANTS_
#define _CONSTANTS_
#include<string>

class Constants{
    private:
        std::string cif_file_path_;
        int z_score_max_;
        int z_score_min_;
        double gap_open_penalty_;
        double gap_ext_penalty_;
        int min_substructure_length_;
        int max_substructure_length_;
    public:
        Constants(std::string cif_file_path, int z_score_min, int z_score_max, double gap_open_penalty, double gap_ext_penalty, int min_substructure_length, int max_substructure_length) :
            cif_file_path_(cif_file_path),
            z_score_min_(z_score_min),
            z_score_max_(z_score_max),
            gap_open_penalty_(gap_open_penalty),
            gap_ext_penalty_(gap_ext_penalty),
            min_substructure_length_(min_substructure_length),
            max_substructure_length_(max_substructure_length){
            
        }
        std::string cifFilePath() const{
            return cif_file_path_;
        }
        int zScoreMax() const{
            return z_score_max_;
        }
        int zScoreMin() const{
            return z_score_min_;
        }
        double gapOpenPenalty() const{
            return gap_open_penalty_;
        }
        double gapExtPenalty() const{
            return gap_ext_penalty_;
        }
        int minSubstructureLength() const{
            return min_substructure_length_;
        }
        int maxSubstructureLength() const{
            return max_substructure_length_;
        }
};

#endif
