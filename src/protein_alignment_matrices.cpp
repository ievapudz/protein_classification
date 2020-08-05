#include "protein_alignment_matrices.h"

//---BLOSUMMatrix class

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

//---DistanceMatrix class

DistanceMatrix::DistanceMatrix(std::vector< std::vector<double> > matrix) : matrix_(matrix){
    
    rows_ = matrix.size();
    columns_ = matrix[0].size();
    
    if(rows_ == 0){
        std::length_error le("Zero length of rows in distance matrix.");
        throw le;
    }
    if(columns_ == 0){
        std::length_error le("Zero length of columns in distance matrix.");
        throw le;
    }
    
}

void DistanceMatrix::displayMatrix(int neighbour_number_left) const{
    //Method that displays the matrix of distances between chain subunits (numbers represent accurate subunit number in the aminoacid chain).
    for(int h = -1; h < columns_; h++){
        if(h == -1)
            std::cout << std::setw(10) << std::right << " ";
        else
            std::cout << std::setw(10) << std::right << h + neighbour_number_left + 1;
    }
    
    std::cout << std::endl;
    
    for(int i = 0; i < rows_; i++){
        std::cout << std::setw(10) << std::right << i + neighbour_number_left + 1;
        for(int j = 0; j < columns_; j++){
            std::cout << std::setw(10) << std::right << matrix_[i][j];
        }
        std::cout << std::endl;
    }
}


void DistanceMatrix::setMatrix(std::vector< std::vector<double> > matrix){
    matrix_ = matrix;
}

std::vector< std::vector<double> > DistanceMatrix::getMatrix() const{
    return matrix_;
}

double DistanceMatrix::getElement(int index_row, int index_column) const{
    return matrix_[index_row][index_column];
}

double DistanceMatrix::getMaxValue(){
    
    double max_value = matrix_[0][0];
    
    for(int i = 0; i < rows_; i++){
        for(int j = 0; j < columns_; j++){
            if(matrix_[i][j] > max_value)
                max_value = matrix_[i][j];
        }
    }
    
    return max_value;
}

double DistanceMatrix::getMinValue(){
    
    double min_value = matrix_[0][0];
    
    for(int i = 0; i < rows_; i++){
        for(int j = 0; j < columns_; j++){
            if(matrix_[i][j] < min_value)
                min_value = matrix_[i][j];
        }
    }
    
    return min_value;
}

double DistanceMatrix::getMean(StatisticCalculator& statistic_calculator){
    
    double mean = statistic_calculator.getMean();
    
    return mean;
}

double DistanceMatrix::getStandardDeviation(StatisticCalculator& statistic_calculator){
    
    double standard_deviation = statistic_calculator.getStandardDeviation();
    
    return standard_deviation;
}

std::vector< std::vector<double> > DistanceMatrix::calculateScoreMatrix(double mean, double standard_deviation){
    /* Method that calculates z-scores for each distance in the distance matrix. Formula: https://en.wikipedia.org/wiki/Standard_score */
    
    std::vector< std::vector<double> > score_matrix(rows_, std::vector<double>(columns_, 0.0));
    
    for(int i = 0; i < rows_; i++){
        for(int j = 0; j < columns_; j++){
            double z_score = ( matrix_[i][j] - mean ) / standard_deviation;
            score_matrix[i][j] = z_score;
        }
    }
    
    return score_matrix;
}

//---DistanceScoreMatrix class

DistanceScoreMatrix::DistanceScoreMatrix(int rows, int columns) : rows_(rows), columns_(columns),
    matrix_(rows, std::vector<double>(columns, 0.0)){
        
    if(rows_ == 0){
        std::length_error le("Zero length of rows in distance score matrix.");
        throw le;
    }
    if(columns_ == 0){
        std::length_error le("Zero length of columns in distance score matrix.");
        throw le;
    }
    
}

DistanceScoreMatrix::DistanceScoreMatrix(std::vector< std::vector<double> > matrix) : matrix_(matrix){
    
    rows_ = matrix.size();
    columns_ = matrix[0].size();
    
    if(rows_ == 0){
        std::length_error le("Zero length of rows in distance score matrix.");
        throw le;
    }
    if(columns_ == 0){
        std::length_error le("Zero length of columns in distance score matrix.");
        throw le;
    }
}

void DistanceScoreMatrix::displayMatrix(int neighbour_number_left) const{
    //Method that displays the matrix of distances between chain subunits (numbers represent accurate subunit number in the aminoacid chain).
    
    for(int h = -1; h < columns_; h++){
        if(h == -1)
            std::cout << std::setw(10) << std::right << " ";
        else
            std::cout << std::setw(10) << std::right << h + neighbour_number_left + 1;
    }
    
    std::cout << std::endl;
    
    for(int i = 0; i < rows_; i++){
        std::cout << std::setw(10) << std::right << i + neighbour_number_left + 1;
        for(int j = 0; j < columns_; j++){
            std::cout << std::setw(10) << std::right << matrix_[i][j];
        }
        std::cout << std::endl;
    }
}

void DistanceScoreMatrix::setMatrix(std::vector< std::vector<double> >& score_matrix){
    matrix_ = score_matrix;
}

double DistanceScoreMatrix::getMaxValue(){
    
    double max_value = matrix_[0][0];
    
    for(int i = 0; i < rows_; i++){
        for(int j = 0; j < columns_; j++){
            if(matrix_[i][j] > max_value)
                max_value = matrix_[i][j];
        }
    }
    
    return max_value;
}

double DistanceScoreMatrix::getMinValue(){
    
    double min_value = matrix_[0][0];
    
    for(int i = 0; i < rows_; i++){
        for(int j = 0; j < columns_; j++){
            if(matrix_[i][j] < min_value)
                min_value = matrix_[i][j];
        }
    }
    
    return min_value;
}

std::vector< std::vector<double> > DistanceScoreMatrix::getMatrix() const{
    return matrix_;
}

double DistanceScoreMatrix::getScore(std::pair<int, int> pair) const{
    // Method that returns score of the corresponding pair of subunits.
    
    double score;
    
    for(int i = 0; i < rows_; i++){
        for(int j = 0; j < columns_; j++){
            if((i+1 == pair.first) && (j+1 == pair.second)){
                score = matrix_[i][j];
                break;
            }
        }
    }
    return score;
}

void DistanceScoreMatrix::limitScores(int z_score_min, int z_score_max){
    for(int i = 0; i < rows_; i++){
        for(int j = 0; j < columns_; j++){
            if(matrix_[i][j] < z_score_min)
                matrix_[i][j] = z_score_min;
            if(matrix_[i][j] > z_score_max)
                matrix_[i][j] = z_score_max;
        }
    }
}

void DistanceScoreMatrix::reverseScoreSigns(){
    for(int i = 0; i < rows_; i++){
        for(int j = 0; j < columns_; j++){
            matrix_[i][j] = matrix_[i][j] * (-1);
        }
    }
}

//---DirectionMatrix class

DirectionMatrix::DirectionMatrix(int rows, int columns) :
    rows_(rows),
    columns_(columns),
    direction_matrix_(rows, std::vector<int>(columns, 0)){
    
}

void DirectionMatrix::setDirection(int row_coordinate, int column_coordinate, int direction){
    direction_matrix_[row_coordinate][column_coordinate] = direction;
}

std::vector< std::vector<int> > DirectionMatrix::returnMatrix() const{
    return this->direction_matrix_;
}

void DirectionMatrix::displayMatrix() const{
    int i, j;
    for (i = 0; i < rows_; i++){
        for (j = 0; j < columns_; j++){
            std::cout << std::setw(4) << std::right << direction_matrix_[i][j];
        }
        std::cout << std::endl;
    }
}

std::vector< std::pair<int, int> > DirectionMatrix::returnNonZeroCoords() const{
    int i, j;
    std::vector< std::pair<int, int> > non_zero_coords;
    for (i = 0; i < rows_; i++){
        for (j = 0; j < columns_; j++){
            if(direction_matrix_[i][j] != 0){
                std::pair<int, int> coordinates;
                coordinates = std::make_pair (i, j);
                non_zero_coords.push_back(coordinates);
            }
        }
    }
    return non_zero_coords;
}

std::vector<int> DirectionMatrix::returnDirections() const{
    int i, j;
    std::vector<int> directions;
    for (i = 0; i < rows_; i++){
        for (j = 0; j < columns_; j++){
            if(direction_matrix_[i][j] != 0){
                directions.push_back(direction_matrix_[i][j]);
            }
        }
    }
    return directions;
}

//---ScoringMatrix class

ScoringMatrix::ScoringMatrix(int rows, int columns, double gap_open_penalty, double gap_ext_penalty) :
    rows_(rows),
    columns_(columns),
    gap_open_penalty_(gap_open_penalty),
    gap_ext_penalty_(gap_ext_penalty),
    scoring_matrix_(rows, std::vector<double>(columns, 0)){
    
}

void ScoringMatrix::setGapOpenPenalty(double gap_open_penalty){
    gap_open_penalty_ = gap_open_penalty;
}

void ScoringMatrix::setGapExtPenalty(double gap_ext_penalty){
    gap_ext_penalty_ = gap_ext_penalty;
}

int ScoringMatrix::getRows() const{
    return rows_;
}

int ScoringMatrix::getColumns() const{
    return columns_;
}

std::vector< std::vector<double> > ScoringMatrix::getScoringMatrix() const{
    return scoring_matrix_;
}

double ScoringMatrix::getGapOpenPenalty() const{
    return gap_open_penalty_;
}

double ScoringMatrix::getGapExtPenalty() const{
    return gap_ext_penalty_;
}

void ScoringMatrix::printScoringMatrix() const{
    int i, j;
    for (i = 0; i < rows_; i++){
        for (j = 0; j < columns_; j++){
            std::cout << std::setw(4) << std::right << scoring_matrix_[i][j];
        }
        std::cout << std::endl;
    }
}

void ScoringMatrix::fillWithGapPenalties(){
    // Method that fills "0th" row and "0th" column with gap extension penalty scores.
    int i, j;
    for(i = 0; i < rows_; i++){
        scoring_matrix_[i][0] = gap_ext_penalty_ * i * (-1);
    }
    for(j = 0; j < columns_; j++){
        scoring_matrix_[0][j] = gap_ext_penalty_ * j * (-1);
    }
}

std::pair<int, int> ScoringMatrix::getPair(int p_subunit, int q_subunit){
    // Method that returns a pair made of subunit numbers in P and Q chains.
    std::pair<int, int> pair = std::make_pair(p_subunit, q_subunit);
    return pair;
}

double ScoringMatrix::getNorthResult(double north_arg, double penalty){
    double north_result = north_arg - penalty;
    return north_result;
}

double ScoringMatrix::getNorthWestResult(double north_west_arg, double score){
    double north_west_result = north_west_arg + score;
    return north_west_result;
}

double ScoringMatrix::getWestResult(double west_arg, double penalty){
    double west_result = west_arg - penalty;
    return west_result;
}

double ScoringMatrix::getInitialMaxOper(double north_arg, double north_west_arg, double west_arg){
    double initial_max_oper = 0;
    
    if(north_arg < initial_max_oper)
       initial_max_oper = north_arg;
    if(north_west_arg < initial_max_oper)
        initial_max_oper = north_west_arg;
    if(west_arg < initial_max_oper)
        initial_max_oper = west_arg;
    
    return initial_max_oper;
}

double ScoringMatrix::getMaxOperationValue(double max_oper, double north, double north_west, double west){
    
    if(north > max_oper)
        max_oper = north;
    if(north_west > max_oper)
        max_oper = north_west;
    if(west > max_oper)
        max_oper = west;
    
    return max_oper;
}

void ScoringMatrix::fillWithScores(char algorithm_choice, std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q, const DistanceScoreMatrix& matrix){
    
    double north, north_west, west;
    int gap_ext_count = 0;
    double max_oper;
    
    for(int i = 1; i < rows_; i++){
        for(int j = 1; j < columns_; j++){
            std::pair<int, int> pair = this->getPair(subunit_chain_P[i-1], subunit_chain_Q[j-1]);
            
            north_west = this->getNorthWestResult(scoring_matrix_[i-1][j-1], matrix.getScore(pair));
            
            if(pair.first != pair.second){
                gap_ext_count++;
                if(gap_ext_count == 1){
                    north = this->getNorthResult(scoring_matrix_[i-1][j], gap_open_penalty_);
                    west = this->getWestResult(scoring_matrix_[i][j-1], gap_open_penalty_);
                }
                if(gap_ext_count > 1){
                    north = this->getNorthResult(scoring_matrix_[i-1][j], gap_ext_penalty_);
                    west = this->getWestResult(scoring_matrix_[i][j-1], gap_ext_penalty_);
                }
            }
            else{
                gap_ext_count = 0;
                north = this->getNorthResult(scoring_matrix_[i-1][j], gap_open_penalty_);
                west = this->getWestResult(scoring_matrix_[i][j-1], gap_open_penalty_);
            }
            
            if(algorithm_choice == '1'){
                max_oper = 0;
            }
            else if(algorithm_choice == '2'){
                max_oper = this->getInitialMaxOper(scoring_matrix_[i][j-1],scoring_matrix_[i-1][j-1], scoring_matrix_[i-1][j]);
            }
            else{
                std::invalid_argument ia("Exception in Scoring_matrix::fillWithScores() invalid algorithm choice.");
                throw ia;
            }
            // Finding and setting the max value from operations.
            scoring_matrix_[i][j] = this->getMaxOperationValue(max_oper, north, north_west, west);
        }
    }
}

void ScoringMatrix::setTracebackStart(char algorithm_choice, double& start_value, int& start_row_coord, int& start_col_coord){
    
    if(algorithm_choice == '1'){
        start_value = 0;
        
        for(int i = 0; i < rows_; i++){
            for(int j = 0; j < columns_; j++){
                if(scoring_matrix_[i][j] > start_value){
                    start_value = scoring_matrix_[i][j];
                    start_row_coord = i;
                    start_col_coord = j;
                }
            }
        }
    }
    if(algorithm_choice == '2'){
        start_row_coord = rows_-1;
        start_col_coord = columns_-1;
        start_value = scoring_matrix_[start_row_coord][start_col_coord];
    }
}

void ScoringMatrix::setTracebackInitialMaxOper(char algorithm_choice, double& max_oper){
    if(algorithm_choice == '1'){
        max_oper = 0;
    }
    if(algorithm_choice == '2'){
        if(rows_ >= columns_) max_oper = scoring_matrix_[rows_-1][0];
        if(rows_ < columns_) max_oper = scoring_matrix_[0][columns_-1];
    }
}

void ScoringMatrix::tracebackNorthOperation(int row_coord, int col_coord, double& north_value, int& north_row_coord, int& north_col_coord, char& zero_coordinate){
    
    if(row_coord == 0){
        north_row_coord = row_coord;
        zero_coordinate = 'r';
    }
    else{
        north_row_coord = row_coord-1;
    }
    
    north_col_coord = col_coord;
    north_value = scoring_matrix_[north_row_coord][north_col_coord];
}

void ScoringMatrix::tracebackNorthWestOperation(int row_coord, int col_coord, double& north_west_value, int& north_west_row_coord, int& north_west_col_coord){
    
    if(row_coord == 0){
        north_west_row_coord = row_coord;
    }
    else{
        north_west_row_coord = row_coord-1;
    }
    
    if(col_coord == 0){
        north_west_col_coord = col_coord;
    }
    else{
        north_west_col_coord = col_coord-1;
    }
    
    north_west_value = scoring_matrix_[north_west_row_coord][north_west_col_coord];
}

void ScoringMatrix::tracebackWestOperation(int row_coord, int col_coord, double& west_value, int& west_row_coord, int& west_col_coord, char& zero_coordinate){
    
    west_row_coord = row_coord;
    if(col_coord == 0){
        west_col_coord = col_coord;
        zero_coordinate = 'c';
    }
    else{
        west_col_coord = col_coord-1;
    }
    
    west_value = scoring_matrix_[west_row_coord][west_col_coord];
}

void ScoringMatrix::tracebackMaxOper(double& max_oper, int& max_value_row_coord, int& max_value_col_coord, double value, int row_coord, int col_coord){
    
    max_oper = value;
    max_value_row_coord = row_coord;
    max_value_col_coord = col_coord;
}

void ScoringMatrix::setTracebackDirections(char zero_coordinate, char indicator, DirectionMatrix& direction_matrix, int west_row_coord, int west_col_coord, int north_row_coord, int north_col_coord, int north_west_row_coord, int north_west_col_coord){
    
    if(zero_coordinate == '-'){
        if(indicator == '1'){
            direction_matrix.setDirection(west_row_coord, west_col_coord+1, WEST);
        }
        if(indicator == '2'){
            direction_matrix.setDirection(north_row_coord+1, north_col_coord, NORTH);
        }
        if(indicator == '3'){
            direction_matrix.setDirection(north_west_row_coord+1, north_west_col_coord+1, NORTH_WEST);
        }
    }
    else{
        switch(zero_coordinate){
            case 'r':
                direction_matrix.setDirection(0, north_col_coord, WEST);
                break;
            case 'c':
                direction_matrix.setDirection(west_row_coord, 0, NORTH);
                break;
            case 'b':
                direction_matrix.setDirection(0, 0, 0);
        }
    }
}

char ScoringMatrix::getIndicator(double& max_oper, int& max_value_row_coord, int& max_value_col_coord, double west_value, int west_row_coord, int west_col_coord, double north_value, int north_row_coord, int north_col_coord, double north_west_value, int north_west_row_coord, int north_west_col_coord){
    
    char indicator = '0';
    
    if(west_value > max_oper){
        this->tracebackMaxOper(max_oper, max_value_row_coord, max_value_col_coord, west_value, west_row_coord, west_col_coord);
        indicator = '1';
    }
    if(north_value > max_oper){
        this->tracebackMaxOper(max_oper, max_value_row_coord, max_value_col_coord, north_value, north_row_coord, north_col_coord);
        indicator = '2';
    }
    if(north_west_value >= max_oper){
        this->tracebackMaxOper(max_oper, max_value_row_coord, max_value_col_coord, north_west_value, north_west_row_coord, north_west_col_coord);
        indicator = '3';
    }
    
    return indicator;
}

DirectionMatrix ScoringMatrix::getDirectionMatrix(char algorithm_choice, std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q){
    
    DirectionMatrix direction_matrix(rows_, columns_);
    
    double max_value;
    int max_value_row_coord, max_value_col_coord;
    double trace;
    
    this->setTracebackStart(algorithm_choice, max_value, max_value_row_coord, max_value_col_coord);
    
    // Continuation of the traceback process
    double north, north_west, west;
    double max_oper;
    int north_row_coord, north_col_coord;
    int north_west_row_coord, north_west_col_coord;
    int west_row_coord, west_col_coord;
    
    bool zero_element = 0;
    char zero_coordinate = '-';
    
    do{
        this->setTracebackInitialMaxOper(algorithm_choice, max_oper);
            
        // North operation
        this->tracebackNorthOperation(max_value_row_coord, max_value_col_coord, north, north_row_coord, north_col_coord, zero_coordinate);
        
        // North-West operation
        this->tracebackNorthWestOperation(max_value_row_coord, max_value_col_coord, north_west, north_west_row_coord, north_west_col_coord);
        
        // West operation
        this->tracebackWestOperation(max_value_row_coord, max_value_col_coord, west, west_row_coord, west_col_coord, zero_coordinate);
        
        // Checking whether the element (0; 0) in the matrix was reached
        if((max_value_col_coord == 0) && (max_value_row_coord == 0)){
            zero_element = 1;
            zero_coordinate = 'b';
        }
        
        // Get indicator to set directions in the direction matrix.
        char indicator = this->getIndicator(max_oper, max_value_row_coord, max_value_col_coord, west, west_row_coord, west_col_coord, north, north_row_coord, north_col_coord, north_west, north_west_row_coord, north_west_col_coord);

        // Directions (degrees) in the matrix show where the cursor for sequence aligning is going to (the character before determines the next one).
        this->setTracebackDirections(zero_coordinate, indicator, direction_matrix, west_row_coord, west_col_coord, north_row_coord, north_col_coord, north_west_row_coord, north_west_col_coord);
        
        trace = max_oper;
        
    }while(!zero_element);
    
    return direction_matrix;
}

void ScoringMatrix::algorithmSmithWaterman(std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q, DistanceScoreMatrix& matrix){
    try{
        this->fillWithScores('1', subunit_chain_P, subunit_chain_Q, matrix);
        this->printScoringMatrix();
        
        DirectionMatrix direction_matrix = this->getDirectionMatrix('1', subunit_chain_P, subunit_chain_Q);
        direction_matrix.displayMatrix();
        
        SequenceAligner seq_al(direction_matrix.returnDirections(), direction_matrix.returnNonZeroCoords(), subunit_chain_P, subunit_chain_Q);
        seq_al.alignSequences();
        
    }catch(std::out_of_range& oor){
        std::cerr << "Out of range: " << oor.what() << std::endl;
    }
    catch(std::invalid_argument& ia){
        std::cerr << ia.what() << std::endl;
    }
}

void ScoringMatrix::algorithmNeedlemanWunsch(std::vector<int> subunit_chain_P, std::vector<int> subunit_chain_Q, DistanceScoreMatrix& matrix){
    this->fillWithGapPenalties();
    try{
        this->fillWithScores('2', subunit_chain_P, subunit_chain_Q, matrix);
        this->printScoringMatrix();
        
        DirectionMatrix direction_matrix = this->getDirectionMatrix('2', subunit_chain_P, subunit_chain_Q);
        direction_matrix.displayMatrix();
            
        SequenceAligner seq_al(direction_matrix.returnDirections(), direction_matrix.returnNonZeroCoords(), subunit_chain_P, subunit_chain_Q);
        seq_al.alignSequences();
            
    }catch(std::out_of_range& oor){
        std::cerr << "Out of range: " << oor.what() << std::endl;
    }
    catch(std::invalid_argument& ia){
        std::cerr << ia.what() << std::endl;
    }
}


