#ifndef _GEOMETRY_
#define _GEOMETRY_
#include<iostream>
#include<math.h>
#define PI 3.14159265

class Vector3D{
        double coordinate_x_;
        double coordinate_y_;
        double coordinate_z_;
    public:
        Vector3D();
        Vector3D(double coordinate_x, double coordinate_y, double coordinate_z);
        Vector3D(Vector3D vector_1, Vector3D vector_2);
        void setX(double coordinate_x);
        void setY(double coordinate_y);
        void setZ(double coordinate_z);
        double getX() const;
        double getY() const;
        double getZ() const;
    
        void displayCoordinates() const;
    
        double getLength();
        Vector3D crossProduct(Vector3D& other);
        double calculateAngle(Vector3D& other);
        Vector3D getUnitVector();
    
        Vector3D operator*(double number);
        double operator*(Vector3D& other);
};

#endif
