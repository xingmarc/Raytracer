//
// template-rt.cpp
//

#define _CRT_SECURE_NO_WARNINGS
#include "matm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

int g_width;
int g_height;

struct Ray
{
    vec4 origin;
    vec4 dir;
};

// TODO: add structs for spheres, lights and anything else you may need.

//SPHERE <name> <pos x> <pos y> <pos z> <scl x> <scl y> <scl z> <r> <g> <b> <Ka> <Kd> <Ks> <Kr> <n>
//LIGHT <name> <pos x> <pos y> <pos z> <Ir> <Ig> <Ib>
struct sphere
{
    string name;
    float pos_x,pos_y,pos_z;
    float scl_x,scl_y,scl_z;
    float r,g,b;
    float Ka,Kd,Ks,Kr;
    float n;
};
vector<sphere> spheres;
struct light
{
    string name;
    float pos_x,pos_y,pos_z;
    float Ir,Ig,Ib;
};
vector<light> lights;

vector<vec4> g_colors;

float g_left;
float g_right;
float g_top;
float g_bottom;
float g_near;

struct backGroundColor{
    float r,g,b;
}back;
struct ambientIntensity{
    float Ir,Ig,Ib;
}ambient;

//output file name:
string outputFileName;



// -------------------------------------------------------------------
// Input file parsing

vec4 toVec4(const string& s1, const string& s2, const string& s3)
{
    stringstream ss(s1 + " " + s2 + " " + s3);
    vec4 result;
    ss >> result.x >> result.y >> result.z;
    result.w = 1.0f;
    return result;
}

float toFloat(const string& s)
{
    stringstream ss(s);
    float f;
    ss >> f;
    return f;
}

void parseLine(const vector<string>& vs)
{
    //TODO: add parsing of NEAR, LEFT, RIGHT, BOTTOM, TOP, SPHERE, LIGHT, BACK, AMBIENT, OUTPUT.
    
    if (vs[0] == "NEAR") {g_near = toFloat(vs[1]);}
    else if(vs[0] == "LEFT"){g_left = toFloat(vs[1]);}
    else if(vs[0] == "RIGHT"){g_right = toFloat(vs[1]);}
    else if(vs[0] == "TOP"){g_top = toFloat(vs[1]);}
    else if(vs[0] == "BOTTOM"){g_bottom = toFloat(vs[1]);}
    
    //RES (Resolution)
    else
    if (vs[0] == "RES")
    {
        g_width = (int)toFloat(vs[1]);
        g_height = (int)toFloat(vs[2]);
        g_colors.resize(g_width * g_height);
    }
    else if(vs[0] == "BACK"){
        back.r = toFloat(vs[1]);
        back.g = toFloat(vs[2]);
        back.b = toFloat(vs[3]);
    }
    else if(vs[0] == "AMBIENT"){
        ambient.Ir = toFloat(vs[1]);
        ambient.Ig = toFloat(vs[2]);
        ambient.Ib = toFloat(vs[3]);
    }
    else if(vs[0] == "SPHERE"){
        sphere s;
        s.name = vs[1];
        s.pos_x = toFloat(vs[2]);
        s.pos_y = toFloat(vs[3]);
        s.pos_z = toFloat(vs[4]);
        s.scl_x = toFloat(vs[5]);
        s.scl_y = toFloat(vs[6]);
        s.scl_z = toFloat(vs[7]);
        s.r = toFloat(vs[8]);
        s.g = toFloat(vs[9]);
        s.b = toFloat(vs[10]);
        s.Ka = toFloat(vs[11]);
        s.Kd = toFloat(vs[12]);
        s.Ks = toFloat(vs[13]);
        s.Kr = toFloat(vs[14]);
        s.n = toFloat(vs[15]);
        
        spheres.push_back(s);
    }
    else if(vs[0] == "LIGHT"){
        light l;
        l.name = vs[1];
        l.pos_x = toFloat(vs[2]);
        l.pos_y = toFloat(vs[3]);
        l.pos_z = toFloat(vs[4]);
        l.Ir = toFloat(vs[5]);
        l.Ig = toFloat(vs[6]);
        l.Ib = toFloat(vs[7]);
        
        lights.push_back(l);
    }
    else if(vs[0] == "OUTPUT"){
        outputFileName = vs[1];
    }
}

void loadFile(const char* filename)
{
    ifstream is(filename);
    if (is.fail())
    {
        cout << "Could not open file " << filename << endl;
        exit(1);
    }
    string s;
    vector<string> vs;
    while(!is.eof())
    {
        vs.clear();
        getline(is, s);
        istringstream iss(s);
        while (!iss.eof())
        {
            string sub;
            iss >> sub;
            vs.push_back(sub);
        }
        parseLine(vs);
    }
}


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec4& color)
{
    int iy2 = g_height - iy - 1; // Invert iy coordinate.
    g_colors[iy2 * g_width + ix] = color;
}


// -------------------------------------------------------------------
// Intersection routine

// TODO: add your ray-sphere intersection routine here.


// -------------------------------------------------------------------
// Ray tracing

vec4 trace(const Ray& ray)
{
    // TODO: implement your ray tracing routine here.
    return vec4(0.0f, 0.0f, 0.0f, 1.0f);
}

vec4 getDir(int ix, int iy)
{
    // TODO: modify this. This should return the direction from the origin
    // to pixel (ix, iy), normalized.
    vec4 dir;
    dir = vec4(0.0f, 0.0f, -1.0f, 0.0f);
    return dir;
}

void renderPixel(int ix, int iy)
{
    Ray ray;
    ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    ray.dir = getDir(ix, iy);
    vec4 color = trace(ray);
    setColor(ix, iy, color);
}

void render()
{
    for (int iy = 0; iy < g_height; iy++)
        for (int ix = 0; ix < g_width; ix++)
            renderPixel(ix, iy);
}


// -------------------------------------------------------------------
// PPM saving

void savePPM(int Width, int Height, char* fname, unsigned char* pixels) 
{
    FILE *fp;
    const int maxVal=255;

    printf("Saving image %s: %d x %d\n", fname, Width, Height);
    fp = fopen(fname,"wb");
    if (!fp) {
        printf("Unable to open file '%s'\n", fname);
        return;
    }
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);

    for(int j = 0; j < Height; j++) {
        fwrite(&pixels[j*Width*3], 3, Width, fp);
    }

    fclose(fp);
}

void saveFile()
{
    // Convert color components from floats to unsigned chars.
    // TODO: clamp values if out of range.
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
            for (int i = 0; i < 3; i++)
                buf[y*g_width*3+x*3+i] = (unsigned char)(((float*)g_colors[y*g_width+x])[i] * 255.9f);
    
    // TODO: change file name based on input file name.
    savePPM(g_width, g_height, "output.ppm", buf);
    delete[] buf;
}



//
void testForDebug ()
{
    //This is the test for my debugging
    //Outputing the parsing result
    ofstream outfile("testForDebug.txt");
    outfile << "NEAR: " <<g_near << endl;
    outfile << "LEFT: " <<g_left << endl;
    outfile << "RIGHT: " << g_right << endl;
    outfile << "BOTTOM: "<<g_bottom<<endl;
    outfile << "TOP: " << g_top<<endl;
    outfile << "RES: "<< g_width <<" "<< g_height<<endl;
    vector<sphere>::iterator sph;
    vector<light>::iterator li;
    
    for (sph = spheres.begin(); sph != spheres.end(); sph++) {
        outfile << "SPHERES: ";
        outfile << sph->name <<" "
        <<sph->pos_x<<" "<<sph->pos_y<<" "<<sph->pos_z<<" "
        <<sph->scl_x<<" "<<sph->scl_y<<" "<<sph->scl_z<<" "
        <<sph->r<<" "<<sph->g<<" "<<sph->b<<" "
        <<sph->Ka<<" "<<sph->Kd<<" "<<sph->Ks<<" "<<sph->Kr<<" "
        <<sph->n<<endl;
    
    }
    for (li = lights.begin(); li != lights.end(); li++) {
        outfile << "LIGHTS: ";
        outfile << li->name <<" "
        <<li->pos_x<<" "<<li->pos_y<<" "<<li->pos_z<<" "
        <<li->Ir << " " <<li->Ig<<" "<<li->Ib<<" "<<endl;
        
        
    }
    outfile <<"BACK: " <<back.r<<" "<<back.g<<" "<<back.b<<endl;
    outfile <<"AMBIENT: "<<ambient.Ir<<" "<<ambient.Ig<<" "<<ambient.Ib<<endl;
    outfile <<"OUTPUT: "<<outputFileName<<endl;
    
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout << "Usage: template-rt <input_file.txt>" << endl;
        exit(1);
    }
    loadFile(argv[1]);
    
    
    //
    testForDebug();
    //
    
    //REMEMBER TO UNCOMMENT THESE TWO LINES:
    //render();
    //saveFile();
	return 0;
}

