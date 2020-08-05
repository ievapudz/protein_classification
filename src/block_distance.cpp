#include "block_distance.h"

BlockDistanceCalculator::BlockDistanceCalculator() : substructure_length_(1){
    
}

BlockDistanceCalculator::BlockDistanceCalculator(Protein protein_1, Protein protein_2) : substructure_length_(1),
    protein_1_(protein_1),
    protein_2_(protein_2),
    phi_angles_1_(protein_1.getPhiAngles()),
    psi_angles_1_(protein_1.getPsiAngles()),
    phi_angles_2_(protein_2.getPhiAngles()),
    psi_angles_2_(protein_2.getPsiAngles()){
        
    
}

void BlockDistanceCalculator::setSubstructureLength(int substructure_length){
    // Method that determines the length of blocks, between which the distance calculations will be made.
    
    if(substructure_length == 0){
        std::invalid_argument ia("Invalid argument for substructure_length in setSubstructureLength.");
        throw ia;
    }
    substructure_length_ = substructure_length;
}

void BlockDistanceCalculator::setPhiAngles(int protein_choice){
    switch (protein_choice) {
        case 1:
            phi_angles_1_ = protein_1_.getPhiAngles();
            break;
        case 2:
            phi_angles_2_ = protein_2_.getPhiAngles();
            break;
        default:
            std::invalid_argument ia("Exception in BlockDistanceCalculator::setPhiAngles: invalid protein choice.");
            throw ia;
            break;
    }
}

void BlockDistanceCalculator::setPsiAngles(int protein_choice){
    switch (protein_choice) {
        case 1:
            psi_angles_1_ = protein_1_.getPsiAngles();
            break;
        case 2:
            psi_angles_2_ = protein_2_.getPsiAngles();
            break;
        default:
            std::invalid_argument ia("Exception in BlockDistanceCalculator::setPhiAngles: invalid protein choice.");
            throw ia;
            break;
    }
}

void BlockDistanceCalculator::setAngles(){
    this->setPhiAngles(1);
    this->setPsiAngles(1);
    this->setPhiAngles(2);
    this->setPsiAngles(2);
}

void BlockDistanceCalculator::setSubunitChain(int chain_choice){
    // A single C_alfa atom is considered as a single subunit.
    // A single aminoacid has got a single C_alfa atom.
    // This method converts the polypeptide chain to the chain of subunits.
    
    std::vector< std::pair<double, double> > subunit_chain;
    
    switch(chain_choice){
        case 1:{
            for(int i = 0; i < this->getNumberOfSubunits(1); i++){
                std::pair<double, double> subunit(phi_angles_1_[i], psi_angles_1_[i]);
                subunit_chain.push_back(subunit);
            }
            subunit_chain_1_ = subunit_chain;
            
            break;
        }
        case 2:{
            for(int i = 0; i < this->getNumberOfSubunits(2); i++){
                std::pair<double, double> subunit(phi_angles_2_[i], psi_angles_2_[i]);
                subunit_chain.push_back(subunit);
            }
            subunit_chain_2_ = subunit_chain;
            
            break;
        }
        default:{
            std::invalid_argument ia("Invalid argument for chain number in setSubunitChain.");
            throw ia;
            break;
        }
    }
}

std::vector<int> BlockDistanceCalculator::getNumberedSubunitChain(int chain_choice){
    // Method that converts the subunit_chain_? to numbered (by the seqID in the corresponding proteins) chain of subunits.
    // Return result will be used for ScoringMatrix::fillWithScores and SequenceAligner::alignSequences methods.
    
    std::vector<int> numbered_subunit_chain;
    
    if(chain_choice == 1){
        for(int i = this->getNeighbourNumberLeft(); i < (subunit_chain_1_.size() - this->getNeighbourNumberRight()); i++){
            numbered_subunit_chain.push_back(i+1);
        }
    }
    else if(chain_choice == 2){
        for(int i = this->getNeighbourNumberLeft(); i < (subunit_chain_2_.size() - this->getNeighbourNumberRight()); i++){
            numbered_subunit_chain.push_back(i+1);
        }
    }
    
    return numbered_subunit_chain;
}

int BlockDistanceCalculator::getSubstructureLength() const{
    return substructure_length_;
}

int BlockDistanceCalculator::getNumberOfSubunits(int chain_choice){
    // A subunit is a C_alfa atom that has a pair of dihedral angles.
    
    int number_of_pairs;
    
    switch (chain_choice) {
        case 1:{
            number_of_pairs = phi_angles_1_.size();
            if(phi_angles_1_.size() > psi_angles_1_.size())
                number_of_pairs = psi_angles_1_.size();
            break;
        }
        case 2:{
            number_of_pairs = phi_angles_2_.size();
            if(phi_angles_2_.size() > psi_angles_2_.size())
                number_of_pairs = psi_angles_2_.size();
            break;
        }
        default:{
            std::invalid_argument ia("Invalid argument for chain number in getNumberOfSubunits.");
            throw ia;
            number_of_pairs = 0;
            break;
        }
    }
    return number_of_pairs;
}

bool BlockDistanceCalculator::isSubstructureLengthEven(){
    if(substructure_length_ % 2 == 0)
        return true;
    else
        return false;
}

int BlockDistanceCalculator::getNeighbourNumberLeft(){
    int neighbour_number_left;
    
    if(this->isSubstructureLengthEven()){
        neighbour_number_left = (substructure_length_ / 2) - 1;
    }else{
        neighbour_number_left = substructure_length_ / 2;
    }
    
    return neighbour_number_left;
}

int BlockDistanceCalculator::getNeighbourNumberRight(){
    int neighbour_number_right;
    
    if(this->isSubstructureLengthEven()){
        neighbour_number_right = substructure_length_ / 2;
    }else{
        neighbour_number_right = this->getNeighbourNumberLeft();
    }
    
    return neighbour_number_right;
}

int BlockDistanceCalculator::getRandomIndex(int min, int max){
    int random_index;
    std::random_device rd;
    do{
        random_index = rd() % max;
    }while(random_index < min);
    return random_index;
}

double BlockDistanceCalculator::getPairDistance(std::pair<double, double> subunit_1, std::pair<double, double> subunit_2, int distance_calculation_choice){
    
    double pair_distance;
    
    switch(distance_calculation_choice){
        case 1:
        {
            // Manhattan distance calculation.
            pair_distance = abs(subunit_1.first - subunit_2.first) + abs(subunit_1.second - subunit_2.second);
            break;
        }
        case 2:
        {
            // Euclidean distance (before taking the square root from the sum) calculation.
            pair_distance = (subunit_1.first - subunit_2.first)*(subunit_1.first - subunit_2.first) + (subunit_1.second - subunit_2.second)*(subunit_1.second - subunit_2.second);
            break;
        }
    }
    return pair_distance;
}

std::vector<double> BlockDistanceCalculator::getRandomDistances(bool both_random){
    
    std::vector<double> distances;
    
    try{
        
        int number_of_subunits_1 = this->getNumberOfSubunits(1);
        int number_of_subunits_2 = this->getNumberOfSubunits(2);
        
        if(number_of_subunits_1 < substructure_length_){
            std::string error_message = "Exception in ";
            error_message.append(protein_1_.getName());
            error_message.append(" BlockDistanceCalculator::getDistances: the chain is shorter than the substructure length.");
            std::length_error le(error_message);
            throw le;
        }
        if(number_of_subunits_2 < substructure_length_){
            std::string error_message = "Exception in ";
            error_message.append(protein_2_.getName());
            error_message.append(" BlockDistanceCalculator::getDistances: the chain is shorter than the substructure length.");
            std::length_error le(error_message);
            throw le;
        }
        
        int neighbour_number_left = this->getNeighbourNumberLeft();
        int neighbour_number_right = this->getNeighbourNumberRight();
        
        this->setSubunitChain(1);
        this->setSubunitChain(2);
        
        if(both_random){
        // Blocks from both chains are taken randomly.
            
            // Index to access a random block in the first chain.
            int l;
            // Index to access a random block in the second chain.
            int k;
            
            for(int i = neighbour_number_left; i < (number_of_subunits_1 - neighbour_number_right); i++){
                
                // Generating random index to identify a block in the first chain.
                l = this->getRandomIndex(neighbour_number_left, number_of_subunits_1 - neighbour_number_right);
                
                // Generating random index to identify a block in the second chain.
                k = this->getRandomIndex(neighbour_number_left, number_of_subunits_2 - neighbour_number_right);
                
                double block_distance = 0;
                
                for(int j = -1 * neighbour_number_left; j <= neighbour_number_right; j++){
                    // Calculating distance between the subunits making up the pair.
                    double pair_distance = this->getPairDistance(subunit_chain_1_[l+j], subunit_chain_2_[k+j], 1);
                    
                    // Calculating the distance between the blocks.
                    block_distance += pair_distance;
                }
                distances.push_back(block_distance);
            }
        }
        else{
        // Blocks only from second chain are taken randomly. Blocks from the first chain are taken orderly.
            
            // Index to access a random block in the second chain.
            int k;
            
            for(int i = neighbour_number_left; i < (number_of_subunits_1 - neighbour_number_right); i++){
                
                // Generating the random index to identify a block in the second chain.
                k = this->getRandomIndex(neighbour_number_left, number_of_subunits_2 - neighbour_number_right);
                
                double block_distance = 0;
                
                for(int j = -1 * neighbour_number_left; j <= neighbour_number_right; j++){
                    // Calculating distance between the subunits making up the pair.
                    double pair_distance = this->getPairDistance(subunit_chain_1_[i+j], subunit_chain_2_[k+j], 1);
                    
                    // Calculating the distance between the blocks.
                    block_distance += pair_distance;
                }
                distances.push_back(block_distance);
            }
        }
    }catch(const std::invalid_argument& ia){
        std::cerr << ia.what() << std::endl;
    }
    catch(const std::length_error& le){
        std::cerr << le.what() << std::endl;
    }
    
    return distances;
}

std::vector< std::vector<double> > BlockDistanceCalculator::getSuccessiveDistances(int distance_calculation_choice){
    
    int number_of_subunits_1 = this->getNumberOfSubunits(1);
    int number_of_subunits_2 = this->getNumberOfSubunits(2);
    
    int neighbour_number_left = this->getNeighbourNumberLeft();
    int neighbour_number_right = this->getNeighbourNumberRight();
    
    int modified_rows = number_of_subunits_1 - (neighbour_number_left + neighbour_number_right);
    int modified_columns = number_of_subunits_2 - (neighbour_number_left + neighbour_number_right);
    
    std::vector< std::vector<double> > distances(modified_rows, std::vector<double>(modified_columns, 0.0));
    
    try{
        if(number_of_subunits_1 < substructure_length_){
            std::string error_message = "Exception in ";
            error_message.append(protein_1_.getName());
            error_message.append(" BlockDistanceCalculator::getDistances: the chain is shorter than the substructure length.");
            std::length_error le(error_message);
            throw le;
        }
        if(number_of_subunits_2 < substructure_length_){
            std::string error_message = "Exception in ";
            error_message.append(protein_2_.getName());
            error_message.append(" BlockDistanceCalculator::getDistances: the chain is shorter than the substructure length.");
            std::length_error le(error_message);
            throw le;
        }
        
        this->setSubunitChain(1);
        this->setSubunitChain(2);
        
        for(int i = neighbour_number_left; i < (number_of_subunits_1 - neighbour_number_right); i++){
            for(int j = neighbour_number_left; j < (number_of_subunits_2 - neighbour_number_right); j++){
                
                double block_distance = 0;
                    
                switch(distance_calculation_choice){
                    case(1):
                    {
                        for(int l = -1 * neighbour_number_left; l <= neighbour_number_right; l++){
                            // Calculating distance between the subunits making up the pair.
                            double pair_distance = this->getPairDistance(subunit_chain_1_[i+l], subunit_chain_2_[j+l], distance_calculation_choice);
                                
                            // Calculating Manhattan distance between the blocks.
                            block_distance += pair_distance;
                        }
                            
                        distances[i-neighbour_number_left][j-neighbour_number_left] = block_distance;
                        break;
                    }
                    case(2):
                    {
                        for(int l = -1 * neighbour_number_left; l <= neighbour_number_right; l++){
                            // Calculating squared coordinate difference between the subunits making up the pair.
                                
                            double pair_distance = this->getPairDistance(subunit_chain_1_[i+l], subunit_chain_2_[j+l], distance_calculation_choice);
                                    
                            // Calculating Euclidean distance between the blocks.
                            block_distance += pair_distance;
                        }
                        distances[i-neighbour_number_left][j-neighbour_number_left] = sqrt(block_distance);
                        break;
                    }
                    default:
                    {
                        std::cerr << "invalid calculation choice\n";
                    }
                }
            }
        }
    }catch(const std::invalid_argument& ia){
        std::cerr << ia.what() << std::endl;
    }
    catch(const std::length_error& le){
        std::cerr << le.what() << std::endl;
    }
    
    return distances;
}
