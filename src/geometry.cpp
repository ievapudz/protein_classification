#include "geometry.h"

//---class Vector3D

Vector3D::Vector3D() : coordinate_x_(0.0), coordinate_y_(0.0), coordinate_z_(0.0){
    
}

Vector3D::Vector3D(double coordinate_x, double coordinate_y, double coordinate_z) : coordinate_x_(coordinate_x), coordinate_y_(coordinate_y), coordinate_z_(coordinate_z){
    
}

Vector3D::Vector3D(Vector3D vector_1, Vector3D vector_2) : coordinate_x_(vector_2.coordinate_x_ - vector_1.coordinate_x_),
    coordinate_y_(vector_2.coordinate_y_ - vector_1.coordinate_y_),
    coordinate_z_(vector_2.coordinate_z_ - vector_1.coordinate_z_){
    
}

void Vector3D::setX(double coordinate_x){
    coordinate_x_ = coordinate_x;
}

void Vector3D::setY(double coordinate_y){
    coordinate_y_ = coordinate_y;
}

void Vector3D::setZ(double coordinate_z){
    coordinate_z_ = coordinate_z;
}

double Vector3D::getX() const{
    return coordinate_x_;
}

double Vector3D::getY() const{
    return coordinate_y_;
}

double Vector3D::getZ() const{
    return coordinate_z_;
}

void Vector3D::displayCoordinates() const{
    std::cout << "(" << coordinate_x_ << "; " << coordinate_y_ << "; " << coordinate_z_ << ")\n";
}

double Vector3D::getLength(){
    double squared_length = this->getX()*this->getX() + this->getY()*this->getY() + this->getZ()*this->getZ();
    return sqrt(squared_length);
}

Vector3D Vector3D::crossProduct(Vector3D& other){
    // Method that performs cross product calculations of the two given vectors and returns a vector that is perpendicular to each of the given two.
    
    // The other vector is the second one in the cross product operation.
    
    double perp_x = coordinate_y_ * other.getZ() - coordinate_z_ * other.getY();
    double perp_y = coordinate_x_ * other.getZ() - coordinate_z_ * other.getX();
    double perp_z = coordinate_x_ * other.getY() - coordinate_y_ * other.getX();
    
    Vector3D perpendicular_vector(perp_x, perp_y, perp_z);
    return perpendicular_vector;
}

double Vector3D::calculateAngle(Vector3D& other){
    double dot_product = *(this) * other;
    double mult_lengths = this->getLength() * other.getLength();
    double cos_angle;
    cos_angle = dot_product / mult_lengths;
    return acos(cos_angle) * 180.0 / PI;
}

Vector3D Vector3D::getUnitVector(){
    double length = this->getLength();
    double unit_x = coordinate_x_ / length;
    double unit_y = coordinate_y_ / length;
    double unit_z = coordinate_z_ / length;
    
    Vector3D unit_vector(unit_x, unit_y, unit_z);
    return unit_vector;
}

Vector3D Vector3D::operator*(double number){
    Vector3D result;
    result.coordinate_x_ = coordinate_x_ * number;
    result.coordinate_y_ = coordinate_y_ * number;
    result.coordinate_z_ = coordinate_z_ * number;
    return result;
}

double Vector3D::operator*(Vector3D& other){
    double result = coordinate_x_ * other.coordinate_x_ + coordinate_y_ * other.coordinate_y_ + coordinate_z_ * other.coordinate_z_;
    return result;
}
