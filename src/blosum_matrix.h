#ifndef _BLOSUM_MATRIX_
#define _BLOSUM_MATRIX_
#include <iostream>
#include <string>
#include <stdexcept> 
#include <map>

class BLOSUMMatrix{
    private:
        std::map<std::string, int> matrix_;
        std::map<std::string, int>::iterator it_;
    public:
        void fillMatrix();
        void displayMatrix();
        int returnScore(const std::string& pair) const;
};

#endif
