#include "traceback.h"

Traceback::Traceback(int rows, int columns, std::vector< std::vector<double> > scoring_matrix) : rows_(rows), columns_(columns), scoring_matrix_(scoring_matrix){
    
}

void Traceback::setAlgorithmChoice(char algorithm_choice){
    algorithm_choice_ = algorithm_choice;
}

char Traceback::getAlgorithmChoice() const{
    return algorithm_choice_;
}

void Traceback::setStart(){
    // Method that sets the begining of the traceback process respectively to the chosen algorithm.
    if(algorithm_choice_ == '1'){
        max_value_ = 0;
        
        for(int i = 0; i < rows_; i++){
            for(int j = 0; j < columns_; j++){
                if(scoring_matrix_[i][j] > max_value_){
                    max_value_ = scoring_matrix_[i][j];
                    max_value_row_coord_ = i;
                    max_value_col_coord_ = j;
                }
            }
        }
    }
    else if(algorithm_choice_ == '2'){
        max_value_row_coord_ = rows_-1;
        max_value_col_coord_ = columns_-1;
        max_value_ = scoring_matrix_[max_value_row_coord_][max_value_col_coord_];
    }
}

void Traceback::setZerothElement(bool is_element_zeroth){
    zeroth_element_ = is_element_zeroth;
}

void Traceback::setZerothCoordinate(char zeroth_coordinate){
    zeroth_coordinate_ = zeroth_coordinate;
}

void Traceback::setInitialMaxOper(){
    // Method that prepares the max_oper_ value for the first set of north, north_west, west operations.
    if(algorithm_choice_ == '1'){
        max_oper_ = 0;
    }
    if(algorithm_choice_ == '2'){
        if(rows_ >= columns_) max_oper_ = scoring_matrix_[rows_-1][0];
        if(rows_ < columns_) max_oper_ = scoring_matrix_[0][columns_-1];
    }
}

void Traceback::northOperation(){
    if(max_value_row_coord_ == 0){
        north_row_coord_ = max_value_row_coord_;
        zeroth_coordinate_ = 'r'; // marking that row coordinate reached zero in the traceback process
    }
    else{
        north_row_coord_ = max_value_row_coord_-1;
    }
    
    north_col_coord_ = max_value_col_coord_;
    north_value_ = scoring_matrix_[north_row_coord_][north_col_coord_];
}

void Traceback::northWestOperation(){
    if(max_value_row_coord_ == 0){
        north_west_row_coord_ = max_value_row_coord_;
    }
    else{
        north_west_row_coord_ = max_value_row_coord_-1;
    }
    
    if(max_value_col_coord_ == 0){
        north_west_col_coord_ = max_value_col_coord_;
    }
    else{
        north_west_col_coord_ = max_value_col_coord_-1;
    }
    
    north_west_value_ = scoring_matrix_[north_west_row_coord_][north_west_col_coord_];
}

void Traceback::westOperation(){
    west_row_coord_ = max_value_row_coord_;
    if(max_value_col_coord_ == 0){
        west_col_coord_ = max_value_col_coord_;
        zeroth_coordinate_ = 'c'; // marking that column coordinate reached zero in the traceback process
    }
    else{
        west_col_coord_ = max_value_col_coord_-1;
    }
    
    west_value_ = scoring_matrix_[west_row_coord_][west_col_coord_];
}

void Traceback::zerothElementCheck(){
    if((max_value_col_coord_ == 0) && (max_value_row_coord_ == 0)){
        zeroth_element_ = 1;
        zeroth_coordinate_ = 'b'; // marking that both coordinates reached zero in the traceback process
    }
}

void Traceback::setProceedingMaxOper(double value, int row_coord, int col_coord){
    max_oper_ = value;
    max_value_row_coord_ = row_coord;
    max_value_col_coord_ = col_coord;
}

void Traceback::setDirectionIndicator(){
    direction_indicator_ = '0';
    
    if(west_value_ > max_oper_){
        this->setProceedingMaxOper(west_value_, west_row_coord_, west_col_coord_);
        direction_indicator_ = '1';
    }
    if(north_value_ > max_oper_){
        this->setProceedingMaxOper(north_value_, north_row_coord_, north_col_coord_);
        direction_indicator_ = '2';
    }
    if(north_west_value_ >= max_oper_){
        this->setProceedingMaxOper(north_west_value_, north_west_row_coord_, north_west_col_coord_);
        direction_indicator_ = '3';
    }
}

void Traceback::setDirection(DirectionMatrix& direction_matrix){
    if(zeroth_coordinate_ == '-'){
        if(direction_indicator_ == '1'){
            direction_matrix.setDirection(west_row_coord_, west_col_coord_+1, WEST);
        }
        if(direction_indicator_ == '2'){
            direction_matrix.setDirection(north_row_coord_+1, north_col_coord_, NORTH);
        }
        if(direction_indicator_ == '3'){
            direction_matrix.setDirection(north_west_row_coord_+1, north_west_col_coord_+1, NORTH_WEST);
        }
    }
    else{
        switch(zeroth_coordinate_){
            case 'r':
                direction_matrix.setDirection(0, north_col_coord_, WEST);
                break;
            case 'c':
                direction_matrix.setDirection(west_row_coord_, 0, NORTH);
                break;
            case 'b':
                direction_matrix.setDirection(0, 0, 0);
        }
    }
}

void Traceback::setTrace(){
    trace_ = max_oper_;
}

DirectionMatrix Traceback::getDirectionMatrix(){
    DirectionMatrix direction_matrix(rows_, columns_);
    
    this->setStart();
    this->setZerothElement(0);
    this->setZerothCoordinate('-');
    
    do{
        this->setInitialMaxOper();
        this->northOperation();
        this->northWestOperation();
        this->westOperation();
        this->zerothElementCheck();
        this->setDirectionIndicator();
        this->setDirection(direction_matrix);
        this->setTrace();
        
    }while(!zeroth_element_);
    
    return direction_matrix;
}

