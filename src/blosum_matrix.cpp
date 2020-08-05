#include "blosum_matrix.h"

void BLOSUMMatrix::fillMatrix(){
    std::vector<char> aminoacids;
    std::vector<char>::iterator it;
    std::vector<char>::iterator it2;
    
    std::ifstream input;
    input.open("BLOSUM_matrix.txt");
    char x;
    do{
        input >> x;
        aminoacids.push_back(x);
    }while(x != '*');
    
    std::string s;
    int y;
    do{
        input >> x;
        for(it = aminoacids.begin(); it != aminoacids.end(); ++it){
            if(x == *it){
                for(it2 = aminoacids.begin(); it2 != aminoacids.end(); ++it2){
                    input >> y;
                    s.push_back(*it);
                    s.push_back(*it2);
                    matrix_.insert( std::pair<std::string, int> (s, y));
                    s.pop_back();
                    s.pop_back();
                }
            }
        }
    }while(!input.eof());
}

void BLOSUMMatrix::displayMatrix(){
    for(it_ = matrix_.begin(); it_ != matrix_.end(); ++it_){
        std::cout << it_->first << "-" << it_->second << " ";
    }
}

int BLOSUMMatrix::returnScore(const std::string& pair) const{
    if(this->matrix_.find(pair) == matrix_.end()){
        std::out_of_range oor("Reached pair that is not present in the BLOSUM matrix.");
        throw oor;
    }
    else return this->matrix_.find(pair)->second;
}

