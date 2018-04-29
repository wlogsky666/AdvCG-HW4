#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <limits>
#include <omp.h>
#include <cstdio>
#include <atomic>
#include <windows.h>

#include "katachi.h"
#include "ppm.h"
#include "algebra3.h"

using namespace std;

struct Config{
    char filename[30] = "in.txt";
    char outputFilename[30] = "output.ppm";
    float eye2Canvas = 0.5;
    int tracingDepth = 2;
    int threadNum = 0;
}myconfig;

int tra;

void readConfig(char* name, Config& config ){
    fstream fp(name, ios::in);
    char tline[128];
    cout << "************** Config **************\n";
    while( fp.getline(tline, 100, '\n') ){
        string line = tline;
        int idx = line.find('=');
        if( idx == std::string::npos){
            cout << "Error config line: " + line << "\n";
        }
        else{
            string parameter = line.substr(0, idx);
            string value = line.substr(idx+1, line.length()-1-idx);
            parameter.erase(parameter.find_last_not_of(" \n\r\t")+1);
            value.erase(value.find_last_not_of(" \n\r\t")+1);
            parameter.erase(0, parameter.find_first_not_of(" \n\r\t"));
            value.erase(0, value.find_first_not_of(" \n\r\t"));

            if(parameter == "filename") {
                strcpy(config.filename, value.c_str());
            }
            else if(parameter == "outputFilename"){
                strcpy(config.outputFilename, value.c_str());
            }
            else if(parameter == "eye2Canvas"){
                float f = stof(value);
                config.eye2Canvas = f;

            }
            else if(parameter == "tracingDepth"){
                int i  = stoi(value);
                config.tracingDepth = i;
            }
            else if(parameter == "threadNum")
            {
                int i  = stoi(value);
                config.threadNum = i;
            }
            else{
                cout << "Error config key: " << parameter << "\n" ;
                continue;
            }
            cout << parameter << ": " << value << endl;
        }
    }
    cout << "\n************* PARALLEL *************\n";
    int max_threads = omp_get_max_threads();
    if(config.threadNum == 0) config.threadNum = max_threads-1;
    cout << "Max number of threads: " << max_threads << endl;
    cout << "Number of thread used: " << config.threadNum << endl;
    omp_set_num_threads(config.threadNum);
}

// mutex printMutex; // MinGW does not have mutex supported
std::atomic_flag lock_stream = ATOMIC_FLAG_INIT;
atomic_int complete(0);
float step;
HANDLE ScreenOut  = GetStdHandle(STD_OUTPUT_HANDLE) ;

void printProgess(){
    while (lock_stream.test_and_set()) {}
    char bar[4] = {'\\', '|', '/', '-'};

    COORD dot = { 14, 25 } ;
    SetConsoleCursorPosition( ScreenOut, dot );

    //printf(" Progress: [");

    SetConsoleTextAttribute(ScreenOut, FOREGROUND_GREEN);
    //for( int i = 0; i <= 100 ; i+=2){
        //if( complete*step >= i )
            //printf("|");
        //else
            //printf(" ");
    //}

    printf("%c", bar[complete%4]);

    SetConsoleTextAttribute(ScreenOut, FOREGROUND_GREEN|FOREGROUND_RED|FOREGROUND_BLUE);
    //printf("] %.2f%     \n", complete * step);

    lock_stream.clear();
}


void output2PPM(Draw image[],
                int width, int height,
                char filename[]);

bool rayTracing(Shader& myShader, Draw& image, vec3 curCenter, vec3 curViewDirection, int depth, double percentage, Grid& all, bool on);

vec3 phongModel(const Shader& myshader,
                const Material& m,
                const vec3& curCenter,
                const vec3& minDisPoint,
                const vec3& normal_vec);


int main(int argc, char* argv[]){
    readConfig("config.txt", myconfig);
    if( argc > 1 ) strcpy(myconfig.filename, argv[1]);
    Shader myShader;
    cout << "\nLoading graphic file...\n";
    readFile(myShader, myconfig.filename);
    tra = myconfig.tracingDepth;

    cout << "Initializing...\n";
    // Find corners of image
    vec3 upSide(0, 1, 0);   // Set y as upside

    vec3 left = upSide ^ myShader.viewDirection;    // Find canvas left
    vec3   up = myShader.viewDirection ^ left;      // Find canvas up

    double eye2canvas = myconfig.eye2Canvas;     // Assume that distance is 0.5 from eye to Canvas

    double width = eye2canvas * tanh(myShader.viewAngle) ;
    double height = width*myShader.r_height/myShader.r_width;

    myShader.Corners[0] = myShader.eyePosition + eye2canvas * myShader.viewDirection;   // Center
    myShader.Corners[1] = myShader.Corners[0] + left * width + up * height;             // Left-Up
    myShader.Corners[2] = myShader.Corners[0] - left * width + up * height;             // Right-Up
    myShader.Corners[3] = myShader.Corners[0] + left * width - up * height;             // Left-Down
    myShader.Corners[4] = myShader.Corners[0] - left * width - up * height;             // Right-Down
    // == END

    vec3 stepWidth  = (myShader.Corners[2]-myShader.Corners[1]) / myShader.r_width ;        // x step for each pixel
    vec3 stepHeight = (myShader.Corners[3]-myShader.Corners[1]) / myShader.r_height;        // y step for each pixel

    vec3 pixelCenter = myShader.Corners[1] + stepWidth/2 + stepHeight/2;                     // center position of each pixel

    Draw* image = new Draw[myShader.r_height*myShader.r_width];
    const int tracingDepth = myconfig.tracingDepth;


    // Calculate p0 and p1
    Grid All;
    for(Sphere& s : myShader.sphs){
        Grid sall = s.getBounding();
        if(sall.p0[0] > All.p0[0]){
            All.p0[0] = sall.p0[0];
        }
        if(sall.p0[1] < All.p0[1]){
            All.p0[1] = sall.p0[1];
        }
        if(sall.p0[2] < All.p0[2]){
            All.p0[2] = sall.p0[2];
        }

        if(sall.p1[0] < All.p1[0]){
            All.p1[0] = sall.p1[0];
        }
        if(sall.p1[1] > All.p1[1]){
            All.p1[1] = sall.p1[1];
        }
        if(sall.p1[2] > All.p1[2]){
            All.p1[2] = sall.p1[2];
        }
    }
    for(Triangle& t : myShader.tris){
        Grid sall = t.getBounding();

        if(sall.p0[0] > All.p0[0]){
            All.p0[0] = sall.p0[0];
        }
        if(sall.p0[1] < All.p0[1]){
            All.p0[1] = sall.p0[1];
        }
        if(sall.p0[2] < All.p0[2]){
            All.p0[2] = sall.p0[2];
        }

        if(sall.p1[0] < All.p1[0]){
            All.p1[0] = sall.p1[0];
        }
        if(sall.p1[1] > All.p1[1]){
            All.p1[1] = sall.p1[1];
        }
        if(sall.p1[2] > All.p1[2]){
            All.p1[2] = sall.p1[2];
        }
    }

    //All.p0[0] += 0.5, All.p0[1] -= 0.5, All.p0[2] -= 0.5;
    //All.p1[0] -= 0.5, All.p1[1] += 0.5, All.p1[2] += 2;
    //cout << All.p0 << All.p1;

    //float splitx = abs(All.p0[0]-All.p1[0]);
    //float splity = abs(All.p0[1]-All.p1[1]);
    //float cx = width / 2;
    //float cy = splity / 2;

    float midx = myShader.eyePosition[0];
    float midy = myShader.eyePosition[1];

    for(Sphere& s : myShader.sphs){
        Grid gs = s.getBounding();
        if(gs.p0[0] >= midx && gs.p1[1] >= midy ){
            myShader.cellSTable[0].push_back(s);
        }
        if(gs.p1[0] <= midx && gs.p1[1] >= midy){
            myShader.cellSTable[1].push_back(s);
        }
        if(gs.p0[0] >= midx && gs.p0[1] <= midy){
            myShader.cellSTable[2].push_back(s);
        }
        if(gs.p1[0] <= midx && gs.p0[1] <= midy){
            myShader.cellSTable[3].push_back(s);
        }
    }
    for(Triangle& t : myShader.tris){
        Grid gs = t.getBounding();
        if(gs.p0[0] >= midx && gs.p1[1] >= midy ){
            myShader.cellTTable[0].push_back(t);
        }
        if(gs.p1[0] <= midx && gs.p1[1] >= midy){
            myShader.cellTTable[1].push_back(t);
        }
        if(gs.p0[0] >= midx && gs.p0[1] <= midy){
            myShader.cellTTable[2].push_back(t);
        }
        if(gs.p1[0] <= midx && gs.p0[1] <= midy){
            myShader.cellTTable[3].push_back(t);
        }
    }

    //cout << myShader.tris.size() << endl;
    //for( int i = 0 ; i < 4 ; ++i )
    //    cout << i << ": " << myShader.cellTTable[i].size() << endl;

    bool on = (myShader.sphs.size()+myShader.tris.size()) > 500;

    // Calculate progress of thread
    step = 100.0 / myShader.r_height;

    cout << "Calculating...\n" ;
    COORD dot = { 0, 25 } ;
    SetConsoleCursorPosition( ScreenOut, dot );
    cout << "Processing ...";
    printProgess();


    // Ray tracing for each Position
    #pragma omp parallel for
    for(int y = 0; y < myShader.r_height; ++y){
        for(int x = 0; x < myShader.r_width; ++x){
            vec3 curCenter = pixelCenter + x*stepWidth + y*stepHeight;
            vec3 curViewDirection = (curCenter-myShader.eyePosition).normalize();       // View Direction of current pixel

            rayTracing(myShader, image[y*myShader.r_width+x], curCenter, curViewDirection, tracingDepth, 1, All, on);
        }
        complete ++;
        printProgess();
    }

    complete = myShader.r_height;
    printProgess();
    cout << "\n[" << myconfig.outputFilename << "] Output to file... \n" ;
    output2PPM(image, myShader.r_width, myShader.r_height, myconfig.outputFilename);

}

bool rayTracing(Shader& myShader, Draw& image, vec3 curCenter, vec3 curViewDirection, int depth, double percentage, Grid& grid, bool on){

    if(depth <= 0) return 0;
    vector<Sphere>* ss = &myShader.sphs;
    vector<Triangle>* tt = &myShader.tris;

    if ( on && depth == tra )
    {
        vec4 line1 = {-curViewDirection[0], 0, grid.p1[0]-grid.p0[0], curCenter[0]-grid.p0[0]};
        vec4 line2 = {-curViewDirection[1], grid.p1[1]-grid.p0[1], 0, curCenter[1]-grid.p0[1]};
        vec4 line3 = {-curViewDirection[2], 0, 0, curCenter[2]-grid.p0[2]};

        double delta = line1[0]*line2[1]*line3[2]+line1[1]*line2[2]*line3[0]+line1[2]*line2[0]*line3[1];
        delta -= line1[0]*line2[2]*line3[1]+line1[2]*line2[1]*line3[0]+line1[1]*line2[0]*line3[2];

        double deltax = line1[3]*line2[1]*line3[2]+line1[1]*line2[2]*line3[3]+line1[2]*line2[3]*line3[1];
        deltax -= line1[3]*line2[2]*line3[1]+line1[1]*line2[3]*line3[2]+line1[2]*line2[1]*line3[3];
        double deltay = line1[0]*line2[3]*line3[2]+line1[3]*line2[2]*line3[0]+line1[2]*line2[0]*line3[3];
        deltay -= line1[0]*line2[2]*line3[3]+line1[3]*line2[0]*line3[2]+line1[2]*line2[3]*line3[0];
        double deltaz = line1[0]*line2[1]*line3[3]+line1[1]*line2[3]*line3[0]+line1[3]*line2[0]*line3[1];
        deltaz -= line1[0]*line2[3]*line3[1]+line1[1]*line2[0]*line3[3]+line1[3]*line2[1]*line3[0];


        double t1 = deltax/delta;
        double s1 = deltay/delta;
        double s2 = deltaz/delta;

        if( !(t1 > 0 && 0 <= s1 && s1 <= 1 && 0 <= s2 && s2 <= 1 )){
            //image.hit = true;
            //image = vec3(0, 255, 0);
            return 0;
        }
        else{

            vec3 node = curCenter + curViewDirection*t1;
            float midx = myShader.eyePosition[0], midy = myShader.eyePosition[1];

            if( node[0] >= midx && node[1] >= midy ){
                ss = &myShader.cellSTable[0];
                tt = &myShader.cellTTable[0];
            }
            else if( node[0] <= midx && node[1] >= midy ){
                ss = &myShader.cellSTable[1];
                tt = &myShader.cellTTable[1];
            }
            else if( node[0] >= midx && node[1] <= midy ){
                ss = &myShader.cellSTable[2];
                tt = &myShader.cellTTable[2];
            }
            else if( node[0] <= midx && node[1] <= midy ){
                ss = &myShader.cellSTable[3];
                tt = &myShader.cellTTable[3];
            }

        }
    }

    double minDis = std::numeric_limits<double>::max();
    vec3 minDisPoint, normal_vec;
    Material minDisM;

    // Draw Spheres
    for(Sphere& s : *ss){
        if( (curCenter-s.o).length() < s.radius ){
            //image[y*myShader.r_width+x] = 1;
            continue;
        }

        // (x-x0)^2 +(y-y0)^2 +(z-z0)^2 = r^2
        // (x+t*xd-x0)^2 + (y+t*yd-y0)^2 + (z+t*zd-z0)^2 = r^2;
        // Check if t is solvable, expand first
        // (x-x0)^2 + 2*t*xd*(x-x0) + (t*xd)^2 + ...... = r^2

        // (xd^2+yd^2+zd^2)*t^2 + (2*xd*(x-x0)+2*yd*(y-y0)+2*zd*(z-z0))*t + ((x-x0)^2+(y-y0)^2+(z-z0)^2-r^2)
        //         a                                  b                                      c
        // check if d >= 0, d = b^2 - 4ac

        double a, b, c, d;
        a = curViewDirection[0]*curViewDirection[0]
                + curViewDirection[1]*curViewDirection[1]
                    + curViewDirection[2]*curViewDirection[2];
        b = 2*curViewDirection[0]*(curCenter[0]-s.o[0])
                + 2*curViewDirection[1]*(curCenter[1]-s.o[1])
                    + 2*curViewDirection[2]*(curCenter[2]-s.o[2]);
        c = (curCenter[0]-s.o[0]) * (curCenter[0]-s.o[0])
                + (curCenter[1]-s.o[1]) * (curCenter[1]-s.o[1])
                    + (curCenter[2]-s.o[2]) * (curCenter[2]-s.o[2])
                        - (s.radius * s.radius);
        d = b*b - 4*a*c;

        // If hit shape
        if( d >= 0 ){
            // image[y*myShader.r_width+x] = 1;

            double t1 = (-b-sqrt(d))/(2*a);      // t1 < t2
            double t2 = (-b+sqrt(d))/(2*a);

            // Need to calculate distance to check the most front shape
            double tmin = t1 >= 0 ? t1 : t2 ;
            double tminDis = (curViewDirection*tmin).length();
            if( minDis > tminDis){
                minDis = tminDis;
                minDisPoint = curCenter + curViewDirection*tmin;
                minDisM = s.m;

                // Compute normal vector of shape for phone-model
                normal_vec = (minDisPoint - s.o).normalize();
            }
        }

    }
    // Draw Triangles
    for(Triangle& t : *tt){

        // x + t*xd = v0_x + s1 * (v1_x- v0_x) + s2 *(v2_x-v0_x)
        // a + t*b  =  c   + S1 *      d       + e  *    s2
        // -t*b + s1*d + s2* e = a-c
        // 解三元一次


        vec4 line1 = {-curViewDirection[0], t.t2[0]-t.t1[0], t.t3[0]-t.t1[0], curCenter[0]-t.t1[0]};
        vec4 line2 = {-curViewDirection[1], t.t2[1]-t.t1[1], t.t3[1]-t.t1[1], curCenter[1]-t.t1[1]};
        vec4 line3 = {-curViewDirection[2], t.t2[2]-t.t1[2], t.t3[2]-t.t1[2], curCenter[2]-t.t1[2]};

        double delta = line1[0]*line2[1]*line3[2]+line1[1]*line2[2]*line3[0]+line1[2]*line2[0]*line3[1];
        delta -= line1[0]*line2[2]*line3[1]+line1[2]*line2[1]*line3[0]+line1[1]*line2[0]*line3[2];

        double deltax = line1[3]*line2[1]*line3[2]+line1[1]*line2[2]*line3[3]+line1[2]*line2[3]*line3[1];
        deltax -= line1[3]*line2[2]*line3[1]+line1[1]*line2[3]*line3[2]+line1[2]*line2[1]*line3[3];
        double deltay = line1[0]*line2[3]*line3[2]+line1[3]*line2[2]*line3[0]+line1[2]*line2[0]*line3[3];
        deltay -= line1[0]*line2[2]*line3[3]+line1[3]*line2[0]*line3[2]+line1[2]*line2[3]*line3[0];
        double deltaz = line1[0]*line2[1]*line3[3]+line1[1]*line2[3]*line3[0]+line1[3]*line2[0]*line3[1];
        deltaz -= line1[0]*line2[3]*line3[1]+line1[1]*line2[0]*line3[3]+line1[3]*line2[1]*line3[0];

        if( delta == 0 ) continue;

        double t1 = deltax/delta;
        double s1 = deltay/delta;
        double s2 = deltaz/delta;

        if( t1 > 0 && s1+s2 <= 1 && 0 <= s1 && s1 <= 1 && 0 <= s2 && s2 <= 1 ){
            //image[y*myShader.r_width+x] = 1;
            double tminDis = (curViewDirection*t1).length();
            if( minDis > tminDis){
                minDis = tminDis;
                minDisPoint = curCenter + curViewDirection*t1;
                minDisM = t.m;

                // Compute normal vector of shape for phone-model
                vec3 N = t.normalVector;

                // Check direction is right, same direction is wrong
                normal_vec = (N * curViewDirection) > 0 ? -N : N;
            }
        }
    }

    // use Phone-Model compute rgb
    if(minDis != std::numeric_limits<double>::max()){
        image.hit = true;
        // Mix color Complementary
        image = vec3(image.rgb.R, image.rgb.G, image.rgb.B) * (1-percentage) + phongModel(myShader, minDisM, curCenter, minDisPoint, normal_vec) * percentage;
        //image = phongModel(myShader, minDisM, curCenter, minDisPoint, normal_vec) * percentage;

        if(minDisM.reflect > 0){
            vec3 newViewDirection = curViewDirection - 2 * (curViewDirection * normal_vec) * normal_vec;
            rayTracing(myShader, image, minDisPoint, newViewDirection, depth-1, minDisM.reflect, grid, on);
        }
    }
    else{
        // Nothing hit by ray
        return false;
    }

    return image.hit;
}

vec3 phongModel(const Shader& myshader,
                const Material& m,
                const vec3& curCenter,
                const vec3& minDisPoint,
                const vec3& normal_vec){

    vec3 Ia = m.rgb;
    vec3 Ii = Ia; // *** Set Ii as material's rgb ***

    vec3 ambient = m.Ka * m.rgb;
    vec3 diffuse = {0,0,0}, specular = {0,0,0};


    for(const auto& light : myshader.lights){
        vec3 lightDirection = (light - minDisPoint).normalize();

        // Diffuse, Kd * Id, Id = Ii * N.L
        if(m.Kd > 0){
            vec3 Id = Ii * (lightDirection * normal_vec);
            if(lightDirection * normal_vec > 0)
                diffuse += m.Kd * Id;
        }

        // Specular, Ks * Is, Is = Ii * (N.H)^n
        if(m.Ks > 0){
            vec3 canvasDirection = (curCenter - minDisPoint).normalize();

            vec3 H = (lightDirection + canvasDirection).normalize();

            vec3 Is = Ii * pow((normal_vec*H), m.exp);
            if(normal_vec*H > 0)
                specular += m.Ks * Is;
        }
    }

    vec3 I = ambient + diffuse + specular;
    return I;
}

void output2PPM(Draw image[],
                int width, int height,
                char filename[]){

    ColorImage output;

    output.init(width, height);
    Pixel white = {255, 255, 255}, black = {0, 0, 0}, red = {255, 0, 0}, green = {0, 255, 0};

    for(int y = 0; y < height; ++y){
        for(int x = 0; x < width; ++x){
            if( image[y*width+x].hit ){
                output.writePixel(x, y, image[y*width+x].rgb);
            }
            else
                {

                }
        }
    }
    output.outputPPM(filename);

}


