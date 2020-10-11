#include"biological_structures.h"

//---class Atom

Atom::Atom() : type_id_("-"), location_(){
    
}

Atom::Atom(std::string type_id, Vector3D location) : type_id_(type_id), location_(location){
    
}

void Atom::setTypeID(std::string type_id){
    type_id_ = type_id;
}

void Atom::setLabelAltID(std::string label_alt_id){
    label_alt_id_ = label_alt_id;
}

void Atom::setLabelAsymID(std::string label_asym_id){
    label_asym_id_ = label_asym_id;
}

void Atom::setLabelEntityID(std::string label_entity_id){
    label_entity_id_ = label_entity_id;
}

void Atom::setLocation(Vector3D& location){
    location_ = location;
}

void Atom::setOccupancy(std::string occupancy){
    occupancy_ = occupancy;
}

void Atom::setAuthSeqID(std::string auth_seq_id){
    auth_seq_id_ = auth_seq_id;
}

void Atom::setAuthCompID(std::string auth_comp_id){
    auth_comp_id_ = auth_comp_id;
}

void Atom::setAuthAsymID(std::string auth_asym_id){
    auth_asym_id_ = auth_asym_id;
}

void Atom::setPdbxPDBModelNum(std::string pdbx_PDB_model_num){
    pdbx_PDB_model_num_ = pdbx_PDB_model_num;
}

std::string Atom::getTypeID() const{
    return type_id_;
}

std::string Atom::getLabelAltID() const{
    return label_alt_id_;
}

std::string Atom::getLabelAsymID() const{
    return label_asym_id_;
}

std::string Atom::getLabelEntityID() const{
    return label_entity_id_;
}

Vector3D Atom::getLocation() const{
    return location_;
}

std::string Atom::getOccupancy() const{
    return occupancy_;
}

std::string Atom::getAuthSeqID() const{
    return auth_seq_id_;
}

std::string Atom::getAuthCompID() const{
    return auth_comp_id_;
}

std::string Atom::getAuthAsymID() const{
    return auth_asym_id_;
}

std::string Atom::getPdbxPDBModelNum() const{
    return pdbx_PDB_model_num_;
}

bool Atom::equalTypeID(std::string type_id) const{
    if(type_id_ == type_id)
        return true;
    else
        return false;
}

bool Atom::validLabelAltID() const{
    if(label_alt_id_ == ".")
        return true;
    else if(label_alt_id_ == " ")
        return true;
    else if(label_alt_id_ == "A")
        return true;
    else
        return false;
}

bool Atom::isOccupancyZero() const{
    int n = occupancy_.length();
    char char_array[n + 1];
    std::strcpy(char_array, occupancy_.c_str());
    double occupancy_double = atof(char_array);
    if(occupancy_double == 0.0){
        return true;
    }
    else
        return false;
}

bool Atom::validAuthAsymID() const{
    if(auth_asym_id_ == "A")
        return true;
    else if(auth_asym_id_ == "AA")
        return true;
    else if(auth_asym_id_ == "I")
        return true;
    else
        return false;
}

bool Atom::equalAuthAsymID(std::string auth_asym_id){
    if(auth_asym_id_ == auth_asym_id)
        return true;
    else
        return false;
}

bool Atom::equalPdbxPDBModelNum(std::string pdbx_PDB_model_num) const{
    if(pdbx_PDB_model_num_ == pdbx_PDB_model_num)
        return true;
    else
        return false;
}

//---class Protein

Protein::Protein() : name_("-"), n_atoms_(), ca_atoms_(), c_atoms_(), sequence_(), subunit_chain_(){
    
}

Protein::Protein(std::string name, std::vector<Atom> all_atoms) : name_(name), all_atoms_(all_atoms), sequence_(), subunit_chain_(){
    
}

Protein::Protein(std::string name, std::vector<Atom> n_atoms, std::vector<Atom> ca_atoms, std::vector<Atom> c_atoms) : name_(name), n_atoms_(n_atoms), ca_atoms_(ca_atoms), c_atoms_(c_atoms), sequence_(), subunit_chain_(){
    
}

void Protein::setName(std::string name){
    name_ = name;
}

void Protein::setAtoms(std::vector<Atom>& atoms){
    try{
        if(atoms[0].getTypeID() == "N")
            n_atoms_ = atoms;
        if(atoms[0].getTypeID() == "CA")
            ca_atoms_ = atoms;
        if(atoms[0].getTypeID() == "C")
            c_atoms_ = atoms;
        else{
            std::string exception_text = "Data is not accepted from this type of atom: ";
            exception_text.append(atoms[0].getTypeID());
            std::out_of_range oor(exception_text);
            throw oor;
        }
    }catch(std::out_of_range& oor){
        std::cerr << "Out of range: " << oor.what() << std::endl;
    }
}

void Protein::setSubunitChain(std::vector<int> subunit_chain){
    subunit_chain_ = subunit_chain;
}

void Protein::setAminoacidSequence(std::vector<std::pair<std::string,std::string> > aminoacid_codes){
    // Method should be invoked after atom filtering.
    
    std::vector<std::string> aminoacid_sequence;
    
    for(int i = 0; i < ca_atoms_.size(); i++){
        for(int j = 0; j < aminoacid_codes.size(); j++){
            if(ca_atoms_[i].getAuthCompID() == aminoacid_codes[j].first)
                aminoacid_sequence.push_back(aminoacid_codes[j].second);
        }
    }
    sequence_ = aminoacid_sequence;
}

void Protein::filterAtoms(std::string filter_atom_type, std::string auth_asym_id){
    
    std::string current_seq_id = all_atoms_.begin()->getAuthSeqID();
    
    std::vector<Atom>::iterator it = all_atoms_.begin();
    
    bool triplet_n_present = 0;
    std::vector<Atom>::iterator n_position = all_atoms_.end();
    
    bool triplet_ca_present = 0;
    std::vector<Atom>::iterator ca_position = all_atoms_.end();
    
    bool triplet_c_present = 0;
    std::vector<Atom>::iterator c_position = all_atoms_.end();
    
    while(it->getAuthSeqID() == current_seq_id){
        if(!(it->isOccupancyZero())){
            if(it->getTypeID() == "N"){
                triplet_n_present = 1;
                n_position = it;
            }
            if(it->getTypeID() == "CA"){
                triplet_ca_present = 1;
                ca_position = it;
            }
            if(it->getTypeID() == "C"){
                triplet_c_present = 1;
                c_position = it;
            }
        }
        it++;
    }
    
    if((!triplet_n_present) || (!triplet_ca_present) || (!triplet_c_present)){
        if((n_position != all_atoms_.end()) && (triplet_n_present)){
            all_atoms_.erase(n_position);
        }
        if((ca_position != all_atoms_.end()) && (triplet_ca_present)){
            all_atoms_.erase(ca_position);
        }
        if((c_position != all_atoms_.end()) && (triplet_c_present)){
            all_atoms_.erase(c_position);
        }
    }
    
    for(int i = 0; i < all_atoms_.size(); i++){
        if((all_atoms_[i].equalTypeID(filter_atom_type)) && (all_atoms_[i].validLabelAltID()) && (all_atoms_[i].equalAuthAsymID(auth_asym_id)) && (all_atoms_[i].equalPdbxPDBModelNum("1")) &&
            (!(all_atoms_[i].isOccupancyZero()))){
            
            if(filter_atom_type == "N")
                n_atoms_.push_back(all_atoms_[i]);
            if(filter_atom_type == "CA")
                ca_atoms_.push_back(all_atoms_[i]);
            if(filter_atom_type == "C")
                c_atoms_.push_back(all_atoms_[i]);
                 
        }
    }
}

std::vector<std::string> Protein::getAminoacidSequence(std::vector<std::pair<std::string,std::string> > aminoacid_codes){
    std::vector<std::string> aminoacid_sequence;
    
    for(int i = 0; i < ca_atoms_.size(); i++){
        for(int j = 0; j < aminoacid_codes.size(); j++){
            if(ca_atoms_[i].getAuthCompID() == aminoacid_codes[j].first)
                aminoacid_sequence.push_back(aminoacid_codes[j].second);
        }
    }
    sequence_ = aminoacid_sequence;
    return aminoacid_sequence;
}

std::string Protein::getName() const{
    return name_;
}

std::vector<Atom> Protein::getAllAtoms() const{
    return all_atoms_;
}

std::vector<Atom> Protein::getAtoms(std::string atom_type) const{
    try{
        if(atom_type == "N")
            return n_atoms_;
        if(atom_type == "CA")
            return ca_atoms_;
        if(atom_type == "C")
            return c_atoms_;
        else{
            std::string exception_text = "There is no data about atom type: ";
            exception_text.append(atom_type);
            std::out_of_range oor(exception_text);
            throw oor;
        }
    }catch(std::out_of_range& oor){
        std::cerr << "Out of range: " << oor.what() << std::endl;
        std::vector<Atom> empty;
        return empty;
    }
}

unsigned int Protein::getResidueNumber() const{
    int residue_number = n_atoms_.size();
    
    if(ca_atoms_.size() < residue_number)
        residue_number = ca_atoms_.size();
    if(c_atoms_.size() < residue_number)
        residue_number = c_atoms_.size();
    
    return residue_number;
}

std::vector<std::string> Protein::getSequence() const{
    return sequence_;
}

std::vector<int> Protein::getSubunitChain() const{
    return subunit_chain_;
}

std::vector< std::vector<Atom> > Protein::getNCACTriplets(){
    std::vector< std::vector<Atom> > n_ca_c_triplets;
    
    unsigned int residue_number = this->getResidueNumber();
    
    // N_i atom - n_ca_c_triplets[i][0]
    // C_alfa_i_atom - n_ca_c_triplets[i][1]
    // C_i_atom - n_ca_c_triplets[i][2]
    
    for(int i = 0; i < residue_number; i++){
        std::vector<Atom> n_ca_c(3);
        n_ca_c[0] = n_atoms_[i];
        n_ca_c[1] = ca_atoms_[i];
        n_ca_c[2] = c_atoms_[i];
        n_ca_c_triplets.push_back(n_ca_c);
    }
    
    return n_ca_c_triplets;
}

std::vector<double> Protein::getPhiAngles(){
    
    // Method implemented with atan2 function (vector coordinate explanation: https://en.wikipedia.org/wiki/Dihedral_angle)
    
    std::vector<double> phi_angles;
    try{
        phi_angles.push_back(360.0);
        
        std::vector< std::vector<Atom> > n_ca_c_triplets = this->getNCACTriplets();
        
        if(n_ca_c_triplets.size() == 0){
            std::string error_message = "Exception in ";
            error_message.append(name_);
            error_message.append(" Protein::getPhiAngles: the number of N, CA, C triplets is zero.");
            std::runtime_error re(error_message);
            throw re;
        }
        
        for(int i = 1; i < n_ca_c_triplets.size(); i++){
            
            Vector3D b1(n_ca_c_triplets[i][0].getLocation(), n_ca_c_triplets[i-1][2].getLocation());
            Vector3D b2(n_ca_c_triplets[i][1].getLocation(), n_ca_c_triplets[i][0].getLocation());
            Vector3D b3(n_ca_c_triplets[i][1].getLocation(), n_ca_c_triplets[i][2].getLocation());
             
            Vector3D cross_b3_b2 = b3.crossProduct(b2);
            Vector3D cross_b1_b2 = b1.crossProduct(b2);
            Vector3D cross_final = cross_b3_b2.crossProduct(cross_b1_b2);
            
            // y is the dot product of b2 and cross_final.
            double y = b2 * cross_final;
            
            Vector3D first_of_x = cross_b1_b2 * b2.getLength();
            
            // x is the dot product of first_of_x and cross_b3_b2.
            double x = first_of_x * cross_b3_b2;
            
            double phi_angle = atan2(y, x) * 180.00 / PI * (-1);
            
            phi_angles.push_back(phi_angle);
        }
    }catch(const std::runtime_error& re){
        std::cerr << re.what() << std::endl;
    }
    
    return phi_angles;
    
}

std::vector<double> Protein::getPsiAngles(){
    
    // Method implemented with atan2 function (vector coordinate explanation: https://en.wikipedia.org/wiki/Dihedral_angle)
    
    std::vector<double> psi_angles;
       
    try{
        std::vector< std::vector<Atom> > n_ca_c_triplets = this->getNCACTriplets();
        
        if(n_ca_c_triplets.size() == 0){
            std::string error_message = "Exception in ";
            error_message.append(name_);
            error_message.append(" Protein::getPsiAngles: the number of N, CA, C triplets is zero.");
            std::runtime_error re(error_message);
            throw re;
        }
           
        for(int i = 0; i < n_ca_c_triplets.size()-1; i++){
            
            Vector3D b1(n_ca_c_triplets[i][1].getLocation(), n_ca_c_triplets[i][0].getLocation());
            Vector3D b2(n_ca_c_triplets[i][2].getLocation(), n_ca_c_triplets[i][1].getLocation());
            Vector3D b3(n_ca_c_triplets[i][2].getLocation(), n_ca_c_triplets[i+1][0].getLocation());
            
            Vector3D cross_b3_b2 = b3.crossProduct(b2);
            Vector3D cross_b1_b2 = b1.crossProduct(b2);
            Vector3D cross_final = cross_b3_b2.crossProduct(cross_b1_b2);
            
            // y is the dot product of b2 and cross_final.
            double y = b2 * cross_final;
            
            Vector3D first_of_x = cross_b1_b2 * b2.getLength();
            
            // x is the dot product of first_of_x and cross_b3_b2.
            double x = first_of_x * cross_b3_b2;
            
            double psi_angle = atan2(y, x) * 180.00 / PI * (-1);
            
            psi_angles.push_back(psi_angle);
        }
        
       if(all_atoms_[0].getTypeID() == "N")
            psi_angles.push_back(360.0);
        
    }catch(const std::runtime_error& re){
        std::cerr << re.what() << std::endl;
    }
    return psi_angles;
}

std::vector<double> Protein::getAngleDistances(unsigned int substructure_length, char choice){
    
    std::vector<double> phi_angles = this->getPhiAngles();
    std::vector<double> psi_angles = this->getPsiAngles();
    
    std::vector<double> angle_distances;
    
    int neighbour_number = substructure_length / 2;
    
    if(choice == 'e'){
        // Euclidean distance
        std::cout << "---Euclidean distance---\n";
        for(int i = 3; i < (this->getResidueNumber() - 3); i++){
            double block_distance = 0;
            for(int j = -1 * neighbour_number; j <= neighbour_number; j++){
                double d_squared = (phi_angles[i+j] - psi_angles[i+j])*(phi_angles[i+j] - psi_angles[i+j]);
                block_distance = block_distance + d_squared;
            }
            std::cout << i+1 << " " << sqrt(block_distance) << std::endl;
            angle_distances.push_back(sqrt(block_distance));
        }
    }
    
    if(choice == 'm'){
        // Manhattan distance
        std::cout << "---Manhattan distance---\n";
        for(int i = 3; i < (this->getResidueNumber() - 3); i++){
            double block_distance = 0;
            for(int j = -1 * neighbour_number; j <= neighbour_number; j++){
                double d_absolute = abs(phi_angles[i+j] - psi_angles[i+j]);
                block_distance = block_distance + d_absolute;
            }
            std::cout << i+1 << " " << block_distance  << std::endl;
            angle_distances.push_back(block_distance);
        }
    }
    
    return angle_distances;
}

