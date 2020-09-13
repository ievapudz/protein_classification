#ifndef _CIFPARSER_
#define _CIFPARSER_
#include<stdlib.h>
#include<fstream>
#include<string>
#include<vector>
#include "biological_structures.h"

class CIFParser{
        std::ifstream input_;
        std::string reading_;
        std::string file_path_;
        std::vector<std::string> atom_site_columns_;
    public:
        void setFilePath(std::string file_path);
        std::string getFilePath() const;
        void parseAtomSiteColumns();
        std::vector<Atom> parseAtoms();
        std::vector<std::string> parseAminoacidSequence();
};

#endif
