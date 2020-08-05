#ifndef _BLOCK_DISTANCE_
#define _BLOCK_DISTANCE_
#include <iostream>
#include <utility>
#include <vector>
#include <stdexcept>
#include <math.h>
#include <random>
#include "biological_structures.h"

class BlockDistanceCalculator{
        int substructure_length_;
    
        Protein protein_1_;
        Protein protein_2_;
    
        std::vector<double> phi_angles_1_;
        std::vector<double> psi_angles_1_;
        std::vector<double> phi_angles_2_;
        std::vector<double> psi_angles_2_;
    
        std::vector< std::pair<double, double> > subunit_chain_1_;
        std::vector< std::pair<double, double> > subunit_chain_2_;
    public:
        BlockDistanceCalculator();
        BlockDistanceCalculator(Protein protein_1, Protein protein_2);
    
        void setSubstructureLength(int substructure_length);
        void setPhiAngles(int protein_choice);
        void setPsiAngles(int protein_choice);
        void setAngles();
        void setSubunitChain(int chain_choice);
    
        std::vector<int> getNumberedSubunitChain(int chain_choice);
        int getSubstructureLength() const;
        int getNumberOfSubunits(int chain_choice);
        bool isSubstructureLengthEven();
        int getNeighbourNumberLeft();
        int getNeighbourNumberRight();
        int getRandomIndex(int min, int max);
    
        double getPairDistance(std::pair<double, double> subunit_1, std::pair<double, double> subunit_2, int distance_calculation_choice);
        
        std::vector<double> getRandomDistances(bool both_random);
    
        std::vector< std::vector<double> > getSuccessiveDistances(int distance_calculation_choice);
       
};

#endif
