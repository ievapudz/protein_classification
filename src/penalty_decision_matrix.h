#ifndef _PENALTY_DECISION_MATRIX_
#define _PENALTY_DECISION_MATRIX_
#include <iostream>
#include <vector>
#include <iomanip>

class PenaltyDecisionMatrix{
private:
    int rows_;
    int columns_;
    std::vector< std::vector<int> > matrix_;
public:
    PenaltyDecisionMatrix(int rows, int columns);
    void fillMatrix();
    void displayMatrix() const;
    std::vector< std::vector<int> > getMatrix();
};

#endif
