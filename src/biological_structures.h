#ifndef _BIOLOGICAL_STRUCTURES_
#define _BIOLOGICAL_STRUCTURES_
#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>
#include "geometry.h"

class Atom{
        std::string type_id_;
        std::string label_alt_id_;
        std::string label_asym_id_;
        std::string label_entity_id_;
        Vector3D location_;
        std::string occupancy_;
        std::string auth_seq_id_;
        std::string auth_comp_id_;
        std::string auth_asym_id_;
        std::string pdbx_PDB_model_num_;
    public:
        Atom();
        Atom(std::string type_id, Vector3D location);
    
        void setTypeID(std::string type_id);
        void setLabelAltID(std::string label_alt_id);
        void setLabelAsymID(std::string label_asym_id);
        void setLabelEntityID(std::string label_entity_id);
        void setLocation(Vector3D& location);
        void setOccupancy(std::string occupancy);
        void setAuthSeqID(std::string auth_seq_id);
        void setAuthCompID(std::string auth_comp_id);
        void setAuthAsymID(std::string auth_asym_id);
        void setPdbxPDBModelNum(std::string pdbx_PDB_model_num);
    
        std::string getTypeID() const;
        std::string getLabelAltID() const;
        std::string getLabelAsymID() const;
        std::string getLabelEntityID() const;
        Vector3D getLocation() const;
        std::string getOccupancy() const;
        std::string getAuthSeqID() const;
        std::string getAuthCompID() const;
        std::string getAuthAsymID() const;
        std::string getPdbxPDBModelNum() const;
    
        bool equalTypeID(std::string type_id) const;
        bool validLabelAltID() const;
        bool isOccupancyZero() const;
        bool validAuthAsymID() const;
        bool equalAuthAsymID(std::string auth_asym_id);
        bool equalPdbxPDBModelNum(std::string pdbx_PDB_model_num) const;
};

class Protein{
        std::string name_;
        std::vector<Atom> all_atoms_;
        std::vector<Atom> n_atoms_;
        std::vector<Atom> ca_atoms_;
        std::vector<Atom> c_atoms_;
        std::vector< std::string > sequence_;
    
    /*  subunit_chain is chain made of aminoacid numbers.
        If substructure length (defined in BlockDistanceCalculator) is one, the subunit_chain is analogous to aminoacid_sequence */
        std::vector<int> subunit_chain_;
    public:
        Protein();
        Protein(std::string name, std::vector<Atom> all_atoms);
        Protein(std::string name, std::vector<Atom> n_atoms, std::vector<Atom> ca_atoms, std::vector<Atom> c_atoms);
        void setName(std::string name);
        void setAllAtoms(std::vector<Atom> all_atoms);
        void setAtoms(std::vector<Atom>& atoms);
        void setSubunitChain(std::vector<int> subunit_chain);
        void setAminoacidSequence(std::vector<std::pair<std::string,std::string> > aminoacid_codes);
    
        void filterAtoms(std::string atom_type, std::string auth_asym_id);
    
        std::vector<std::string> getAminoacidSequence(std::vector<std::pair<std::string,std::string> > aminoacid_codes);
    
        std::string getName() const;
        std::vector<Atom> getAllAtoms() const;
        std::vector<Atom> getAtoms(std::string atom_type) const;
        unsigned int getResidueNumber() const;
        std::vector<std::string> getSequence() const;
        std::vector<int> getSubunitChain() const;
    
        std::vector< std::vector<Atom> > getNCACTriplets();
    
        std::vector<double> getPhiAngles();
        std::vector<double> getPsiAngles();
    
        std::vector<double> getAngleDistances(unsigned int substructure_length, char choice);
    
        void reset();
};

#endif
