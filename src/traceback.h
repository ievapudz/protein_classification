#ifndef _TRACEBACK_
#define _TRACEBACK_
#include <iostream>
#include <vector>
#include "direction_matrix.h"
#define NORTH 360
#define NORTH_WEST 315
#define WEST 270

class Traceback{
    private:
        int rows_;
        int columns_;
        std::vector< std::vector<double> > scoring_matrix_;
        
        char algorithm_choice_;
        
        double max_value_;
        int max_value_row_coord_;
        int max_value_col_coord_;
    
        double trace_;
        
        double max_oper_;
        int max_oper_row_coord_;
        int max_oper_col_coord_;
        
        double north_value_;
        int north_row_coord_;
        int north_col_coord_;
        
        double north_west_value_;
        int north_west_row_coord_;
        int north_west_col_coord_;
        
        double west_value_;
        int west_row_coord_;
        int west_col_coord_;
        
        bool zeroth_element_;
        char zeroth_coordinate_;
    
        char direction_indicator_;
    public:
        Traceback(int rows, int columns, std::vector< std::vector<double> > scoring_matrix);
        void setAlgorithmChoice(char algorithm_choice);
        char getAlgorithmChoice() const;
    
        void setStart();
        void setZerothElement(bool is_element_zeroth);
        void setZerothCoordinate(char zeroth_coordinate);
        void setInitialMaxOper();
        void northOperation();
        void northWestOperation();
        void westOperation();
        void zerothElementCheck();
        void setProceedingMaxOper(double value, int row_coord, int col_coord);
        void setDirectionIndicator();
        void setDirection(DirectionMatrix& direction_matrix);
        void setTrace();
    
        DirectionMatrix getDirectionMatrix();
};

#endif
