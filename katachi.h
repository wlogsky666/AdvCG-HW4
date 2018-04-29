#ifndef KATACHI_H_INCLUDED
#define KATACHI_H_INCLUDED

#include <iostream>
#include <vector>

#include "algebra3.h"
#include "ppm.h"

// COUT & CIN METHOD : vec3
std::istream& operator >> (std::istream& is, vec3& v){
    double x, y, z;
    is >> x >> y >> z;
    v.set(x, y, z);
    return is;
}
std::ostream& operator << (std::ostream& os,  vec3& v){
    os << v[0] << " " << v[1] << " " << v[2] << std::endl;
    return os;
}

// Material
struct Material{
    vec3 rgb;
    double Ka, Kd, Ks;
    double exp;
    double reflect, refract;
    double nr;
}tmpm;
std::istream& operator >> (std::istream& is, Material& m){
    is >> m.rgb >> m.Ka >> m.Kd >> m.Ks >> m.exp >> m.reflect >> m.refract >> m.nr;
    return is;
}

// Grid : Regular Grids -> Speed up ray Tracing
struct Grid{
    vec3 p0 = {-1000, 1000, 1000};
    vec3 p1 = {1000, -1000, -1000};
};

// *************** Shape ***************

// Triangle
struct Triangle{
    vec3 t1, t2, t3;
    vec3 normalVector;
    Material m;

    Grid getBounding(){
        Grid All;
        // p0
        if(t1[0] > All.p0[0]){
            All.p0[0] = t1[0];
        }
        if(t1[1] < All.p0[1]){
            All.p0[1] = t1[1];
        }
        if(t1[2] < All.p0[2]){
            All.p0[2] = t1[2];
        }
        if(t2[0] > All.p0[0]){
            All.p0[0] = t2[0];
        }
        if(t2[1] < All.p0[1]){
            All.p0[1] = t2[1];
        }
        if(t2[2] < All.p0[2]){
            All.p0[2] = t2[2];
        }
        if(t3[0] > All.p0[0]){
            All.p0[0] = t3[0];
        }
        if(t3[1] < All.p0[1]){
            All.p0[1] = t3[1];
        }
        if(t3[2] < All.p0[2]){
            All.p0[2] = t3[2];
        }

        // p1
        if(t1[0] < All.p1[0]){
            All.p1[0] = t1[0];
        }
        if(t1[1] > All.p1[1]){
            All.p1[1] = t1[1];
        }
        if(t1[2] > All.p1[2]){
            All.p1[2] = t1[2];
        }
        if(t2[0] < All.p1[0]){
            All.p1[0] = t2[0];
        }
        if(t2[1] > All.p1[1]){
            All.p1[1] = t2[1];
        }
        if(t2[2] > All.p1[2]){
            All.p1[2] = t2[2];
        }
        if(t3[0] < All.p1[0]){
            All.p1[0] = t3[0];
        }
        if(t3[1] > All.p1[1]){
            All.p1[1] = t3[1];
        }
        if(t3[2] > All.p1[2]){
            All.p1[2] = t3[2];
        }

        if(All.p0[0] == All.p1[0]){
            All.p0[0] += 0.05;
            All.p1[0] -= 0.05;
        }
        if(All.p0[1] == All.p1[1]){
            All.p0[1] += 0.05;
            All.p1[1] -= 0.05;
        }
        return All;
    }

}tmpt;
std::istream& operator >> (std::istream& is, Triangle& triangle){
    is >> triangle.t1 >> triangle.t2 >> triangle.t3 >> triangle.normalVector;
    return is;
}

// Sphere
struct Sphere{
    vec3 o;
    double radius;
    Material m;

    Grid getBounding(){
        Grid All;
        // p0
        if(o[0] + radius > All.p0[0]){
            All.p0[0] = o[0] + radius;
        }
        if(o[1] - radius < All.p0[1]){
            All.p0[1] = o[1] - radius;
        }
        if(o[2] - radius < All.p0[2]){
            All.p0[2] = o[2] - radius;
        }
        // p1
        if(o[0] - radius < All.p1[0]){
            All.p1[0] = o[0] - radius;
        }
        if(o[1] + radius > All.p1[1]){
            All.p1[1] = o[1] + radius;
        }
        if(o[2] + radius > All.p1[2]){
            All.p1[2] = o[2] + radius;
        }

        if(All.p0[0] == All.p1[0]){
            All.p0[0] += 0.05;
            All.p1[0] -= 0.05;
        }
        if(All.p0[1] == All.p1[1]){
            All.p0[1] += 0.05;
            All.p1[1] -= 0.05;
        }
        return All;
    }
}tmps;
std::istream& operator >> (std::istream& is, Sphere& sphere){
    is >> sphere.o >> sphere.radius;
    return is;
}


struct Shader{

    vec3 eyePosition;
    vec3 viewDirection;
    double viewAngle;
    int r_width, r_height;
    vec3 Corners[5];        // center, leftUp, rightUp, leftDown, rightDown

    // Store multi shapes
    std::vector<Triangle> tris;
    std::vector<Sphere> sphs;
    std::vector<vec3> lights;
    std::vector<Sphere> cellSTable[4];
    std::vector<Triangle> cellTTable[4];
};


struct Draw{
    Pixel rgb = {0, 0, 0};
    bool hit = false;
    bool operator = (vec3 v){
        if(v[0] < 0) v[0] = 0;
        else if(v[0] > 255) v[0] = 255;
        if(v[1] < 0) v[1] = 0;
        else if(v[1] > 255) v[1] = 255;
        if(v[2] < 0) v[2] = 0;
        else if(v[2] > 255) v[2] = 255;

        rgb.R = v[0], rgb.G = v[1], rgb.B = v[2];
        return true;
    }


};


bool readFile(Shader& myShader, char filename[]){
    std::fstream fp(filename, std::ios::in);
    if( !fp ){
        std::cerr << "[" << filename << "] File not exists\n";
        return false;
    }

    std::string line;
    while(fp >> line){
        if(line == "E"){
            fp >> myShader.eyePosition;
        }
        else if(line == "V"){
            fp >> myShader.viewDirection;
            myShader.viewDirection.normalize();
        }
        else if(line == "F"){
            fp >> myShader.viewAngle;
        }
        else if(line == "R"){
            fp >> myShader.r_width >> myShader.r_height ;
        }

        else if(line == "S"){
            fp >> tmps;
            tmps.m = tmpm;
            myShader.sphs.push_back(tmps);
        }
        else if(line == "T"){
            fp >> tmpt;
            tmpt.m = tmpm;
            myShader.tris.push_back(tmpt);
        }
        else if(line == "M"){
            fp >> tmpm;
            tmpm.rgb *= 255;
        }
        else if(line == "L"){
            vec3 light;
            fp >> light;
            myShader.lights.push_back(light);
        }
        else{
            std::cerr << "Invalid char\n";
        }
    }
    return true;
}

#endif // KATACHI_H_INCLUDED
