#include "cif_parser.h"

//---class CIFParser

void CIFParser::setFilePath(std::string file_path){
    file_path_ = file_path;
}

std::string CIFParser::getFilePath() const{
    return file_path_;
}

void CIFParser::parseAtomSiteColumns(){
    input_.open(file_path_);
    while(!input_.eof()){
        input_ >> reading_;
        if(reading_ == "_atom_site.group_PDB"){
            atom_site_columns_.push_back(reading_);
            do{
                input_ >> reading_;
                if(reading_ != "ATOM"){
                    atom_site_columns_.push_back(reading_);
                }
            }while(reading_ != "_atom_site.pdbx_PDB_model_num");
            break;
        }
    }
}

std::vector<Atom> CIFParser::parseAtoms(){
    std::vector<Atom> all_atoms;
    std::vector<std::string> atom_row(atom_site_columns_.size());
    
    do{
        Atom atom;
        Vector3D location;
        for(int i = 0; i < atom_site_columns_.size(); i++){
            input_ >> reading_;
            if(reading_ == "#")
                break;
            atom_row[i] = reading_;
            if(atom_site_columns_[i] == "_atom_site.label_atom_id"){
                atom.setTypeID(atom_row[i]);
            }
            if(atom_site_columns_[i] == "_atom_site.label_alt_id"){
                atom.setLabelAltID(atom_row[i]);
            }
            if(atom_site_columns_[i] == "_atom_site.label_asym_id"){
                atom.setLabelAsymID(atom_row[i]);
            }
            if(atom_site_columns_[i] == "_atom_site.label_entity_id"){
                atom.setLabelEntityID(atom_row[i]);
            }
            double coordinate_x, coordinate_y, coordinate_z;
            if(atom_site_columns_[i] == "_atom_site.Cartn_x"){
                const char * c = atom_row[i].c_str();
                coordinate_x = atof(c);
                location.setX(coordinate_x);
            }
            if(atom_site_columns_[i] == "_atom_site.Cartn_y"){
                const char * c = atom_row[i].c_str();
                coordinate_y = atof(c);
                location.setY(coordinate_y);
            }
            if(atom_site_columns_[i] == "_atom_site.Cartn_z"){
                const char * c = atom_row[i].c_str();
                coordinate_z = atof(c);
                location.setZ(coordinate_z);
            }
            atom.setLocation(location);
            if(atom_site_columns_[i] == "_atom_site.occupancy"){
                atom.setOccupancy(atom_row[i]);
            }
            if(atom_site_columns_[i] == "_atom_site.auth_seq_id"){
                atom.setAuthSeqID(atom_row[i]);
            }
            if(atom_site_columns_[i] == "_atom_site.auth_comp_id"){
                atom.setAuthCompID(atom_row[i]);
            }
            if(atom_site_columns_[i] == "_atom_site.auth_asym_id"){
                atom.setAuthAsymID(atom_row[i]);
            }
            if(atom_site_columns_[i] == "_atom_site.pdbx_PDB_model_num"){
                atom.setPdbxPDBModelNum(atom_row[i]);
            }
        }
        if(atom.getTypeID() != "-")
            all_atoms.push_back(atom);
    }while(reading_ != "#");
    
    return all_atoms;
}

