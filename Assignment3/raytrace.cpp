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


const float EPSILON = 0.0001;



// -------------------------------------------------------------------
//The data structures listed here:

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
    
    // you can also set this: Garrett's in the discussion:
    // mat4 transform // Or make it faster: mat4 inverse_transform
    // vec4 color
    // float etc
    // sphere (vec4 color, float etc,...)// the construction function
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




vec4 addIllumination(float t, Ray ray, sphere sph,light li);

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
        
        //the size of the vector g_colors = the number of the pixels
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
//
bool calculateTheRaysphereIntersection(vec3 S, vec3 c,float &t_near, float &t_far){
    //for these two parameters: vec3
    float t1, t2;
    float B = dot(S, c);
    float A = length(c)*length(c);
    float C = length(S)*length(S) - 1;
    
    float delta = B*B - A*C;
    
    if (delta < 0) {
        
        return false;// no solution case, return -1;
    }
    else{
        t1 = -1*B/A + sqrt(delta)/A;
        t2 = -1*B/A - sqrt(delta)/A;
        t_near = t2;
        t_far = t1;
        //
        return true;
    }

}

// -------------------------------------------------------------------
// Ray tracing

vec3 truncateVec4ToVec3(vec4 input){
    return vec3(input.x,input.y,input.z);
}

mat4 calculateInverseM (sphere s){
    mat4 inv_M;
    
    bool result = InvertMatrix(mat4(s.scl_x,0,0,0,
                                    0,s.scl_y,0,0,
                                    0,0,s.scl_z,0,
                                    s.pos_x,s.pos_y,s.pos_z,1),
                               inv_M);
    
    if (!result) {
        cerr << "Matrix un-invertable!"<<endl;
        return mat4(0.0f);//return 0
    }
    return inv_M;
}



// if the return value is true, there is a shadow,
// if the return value is false, there is not a shadow
bool addShadow(float t, Ray r, light li){
    
    vec4 S = r.origin;
    vec4 c_ray = r.dir;
    
    vec4 Start = S + t*c_ray;
    
    
    //vector<light> ::iterator itl;
    vector<sphere>::iterator isph;
    
    
    
    vec4 c = vec4( li.pos_x - Start.x,
                      li.pos_y - Start.y,
                      li.pos_z - Start.z,0);
        
        
    for (isph = spheres.begin(); isph != spheres.end(); isph++) {
        vec4 S_prime_vec4 = calculateInverseM(*isph) * Start;
        vec4 c_prime_vec4 = calculateInverseM(*isph) * c;
            
        vec3 S_prime = truncateVec4ToVec3(S_prime_vec4);
        vec3 c_prime = truncateVec4ToVec3(c_prime_vec4);
        
        float t1, t2;
        //float t;
        if(calculateTheRaysphereIntersection(S_prime, c_prime, t1, t2))
        {
            //
            if ((t1 > EPSILON && t1 < 1)||(t2 > EPSILON && t2 < 1 )) {
                return true;
            }
        }
    }
    return false;
}


//This function get the nearest sphere the ray encounters
bool rayIntersectWithSphere(sphere &sph, Ray r, float &t, bool eye)
// this eye checks if the ray starts from the eye point
{
    
    // bool result = false;
    //the iterator for spheres
    vector<float> vecForT;
    vector<float>:: iterator it;
    vector<sphere> vecForSphere;
    
    
    vector<sphere>::iterator s_iterator;
    
    vec4 S_prime_vec4, c_prime_vec4;
    vec3 S_prime,c_prime;
    mat4 inv_M;
    
    
    
    sphere temp_sph;
    float temp_t;
    //int flag = 0;
    
    // iterate all the spheres
    for (s_iterator = spheres.begin(); s_iterator != spheres.end(); s_iterator++) {
        
        temp_sph = *s_iterator;
        
        inv_M = calculateInverseM(*s_iterator);
        
        S_prime_vec4 = inv_M * r.origin;
        c_prime_vec4 = inv_M * r.dir;
        
        //truncate to vec3
        S_prime = truncateVec4ToVec3(S_prime_vec4);
        c_prime = truncateVec4ToVec3(c_prime_vec4);
        
        
        float tt_near,tt_far;
        // there is an intersection
        if (calculateTheRaysphereIntersection(S_prime, c_prime, tt_near, tt_far))
        {
            if (eye) {
                if (tt_near >= 1) {
                    vecForT.push_back(tt_near);
                    vecForSphere.push_back(*s_iterator);
                }
                else{
                    if (tt_far >= 1) {
                        vecForT.push_back(tt_far);
                        vecForSphere.push_back(*s_iterator);
                    }
                }
            }
            else
            {
                if (tt_near > 0.0001f) {
                    vecForT.push_back(tt_near);
                    vecForSphere.push_back(*s_iterator);
                }
                
            }
            
        }
        //else {flag ++;}
    }
    
    
    if (vecForT.empty()) {
        return false;
    }
    else
    {
        
        int i = 0;
        
        int num = 0;
        
        temp_t = vecForT[0];
        
        for (it = vecForT.begin(); it != vecForT.end(); it++,i++) {
            if (*it < temp_t) {
                temp_t = *it;
                num = i;
            }
        }
        
        sph = vecForSphere[num];
        t = vecForT[num];
        //cout << t<<" ";
        return  true;
        
        
    }
}




vec4 reflection (Ray ray, float t, sphere sph, int depth)
{
    mat4 inv_M = calculateInverseM(sph);
    //vec4 color = vec4(0,0,0,1);
    
    vec3 S = truncateVec4ToVec3(ray.origin);
    vec3 c = truncateVec4ToVec3(ray.dir);
    
    vec3 S_prime = truncateVec4ToVec3(inv_M * ray.origin);
    vec3 c_prime = truncateVec4ToVec3(inv_M * ray.dir);
    
    //get the hit point
    vec3 hitPoint = S + t*c;
    
    //get the normal vector:
    vec3 hitPointOnTheUnitSphere = S_prime + c_prime *t;
    
    vec3 normalOnTheUnitSphere = normalize(hitPointOnTheUnitSphere);
    vec4 normalOnTheUnitSphere_vec4 =
    vec4(normalOnTheUnitSphere.x,normalOnTheUnitSphere.y,
         normalOnTheUnitSphere.z,0);
    
    //the normal vector = inv_M * the normal vector on the unit untransformed sphere
    
    vec4 normal_vec4 = normalize(inv_M * normalOnTheUnitSphere_vec4);
    //the last component is 0, so I can truncate it into vec3.
    vec3 normal = vec3(normal_vec4.x,normal_vec4.y,normal_vec4.z);
    
    
    
    
    vec3 in = normalize(vec3(ray.origin.x - hitPoint.x,
                            ray.origin.y - hitPoint.y,
                            ray.origin.z - hitPoint.z));
    
    vec3 R = normalize( 2*normal*dot(normal, in) - in );
    
    Ray newRay;
    newRay.origin = vec4(hitPoint.x,hitPoint.y,hitPoint.z,1);
    newRay.dir = vec4(R.x,R.y,R.z,0);
    
    sphere intersectedSphere;
    float intersectedT;
    
    //if the ray intersected with other spheres:
    if (!rayIntersectWithSphere(intersectedSphere, newRay, intersectedT, false)) {
        return vec4 (0,0,0,1);
    }
    else{

    if (depth == 0)
    {
        
        //cout << intersectedT <<"newRay :"<< newRay.origin <<" "<<newRay.dir<<endl;
        //cout <<"vec4:"<< vec4(sph.r, sph.g, sph.b, 1.0)<<endl;
        vec4 temp_color(0,0,0,1);
        vector<light> ::iterator itl;
        for (itl = lights.begin(); itl != lights.end(); itl++)
        {
            temp_color += addIllumination(intersectedT, newRay, intersectedSphere,*itl);
        }
        
        
        
        return vec4(temp_color.x + ambient.Ir * intersectedSphere.Ka * intersectedSphere.r,
                    temp_color.y + ambient.Ig * intersectedSphere.Ka * intersectedSphere.g,
                    temp_color.z + ambient.Ib * intersectedSphere.Ka * intersectedSphere.b,
                    1.0);
        //return vec4(sph.r * sph.Ka, sph.g  * sph.Ka, sph.b * sph.Ka, 1.0);//*Ka????
        
    }
    else{
        depth++;
        return reflection(newRay, intersectedT, intersectedSphere, depth);
    }
    }
}




//int count1 = 0;
//Add Illumination for one light source
vec4 addIllumination(float t, Ray ray, sphere sph,light li)
                     //vec3 S, vec3 c, vec3 S_prime, vec3 c_prime,
                     //sphere sph, light li)
{
    mat4 inv_M = calculateInverseM(sph);
    vec4 color = vec4(0,0,0,1);
    
    vec3 S = truncateVec4ToVec3(ray.origin);
    vec3 c = truncateVec4ToVec3(ray.dir);
    
    vec3 S_prime = truncateVec4ToVec3(inv_M * ray.origin);
    vec3 c_prime = truncateVec4ToVec3(inv_M * ray.dir);
    
    //get the hit point
    vec3 hitPoint = S + t*c;
    
    //get the normal vector:
    vec3 hitPointOnTheUnitSphere = S_prime + c_prime *t;
    
    vec3 normalOnTheUnitSphere = normalize(hitPointOnTheUnitSphere);
    vec4 normalOnTheUnitSphere_vec4 =
           vec4(normalOnTheUnitSphere.x,normalOnTheUnitSphere.y,
                normalOnTheUnitSphere.z,0);
    
    //the normal vector = inv_M * the normal vector on the unit untransformed sphere
    
    vec4 normal_vec4 = normalize(inv_M * normalOnTheUnitSphere_vec4);
    //the last component is 0, so I can truncate it into vec3.
    vec3 normal = vec3(normal_vec4.x,normal_vec4.y,normal_vec4.z);
    
    
    //vector<light> ::iterator itl;
    //for (itl = lights.begin(); itl != lights.end(); itl++) {
        vec3 L = normalize(vec3(li.pos_x - hitPoint.x,
                                li.pos_y - hitPoint.y,
                                li.pos_z - hitPoint.z));
        
        vec3 R = normalize( 2*normal*dot(normal, L) - L );
        
        vec3 v = normalize(vec3(S.x - hitPoint.x,
                                S.y - hitPoint.y,
                                S.z - hitPoint.z));
        
        float cosTheta = dot(normal,L);
        float cosPhi_n = pow(dot(R, v),sph.n);
    
    
   // cout << " RV : "<< dot(R, v) << endl;
    
    
        if (cosTheta > 0) {//if we can't see the light
            //color.x += ambient.Ir * sph.Ka;
            //color.y += ambient.Ig * sph.Ka;
            //color.z += ambient.Ib * sph.Ka;

            if (dot(R,v) > 0) {
                color.x += li.Ir * sph.Kd * cosTheta * sph.r + li.Ir * sph.Ks * cosPhi_n;
                //ambient.Ir*sph.Ka;
                color.y += li.Ig * sph.Kd * cosTheta * sph.g + li.Ig * sph.Ks * cosPhi_n ;
                //ambient.Ig * sph.Ka;
                color.z += li.Ib * sph.Kd * cosTheta * sph.b + li.Ib * sph.Ks * cosPhi_n;
                //ambient.Ib * sph.Ka;

                
            }
            else{
                color.x += li.Ir * sph.Kd * cosTheta * sph.r ;
                //ambient.Ir*sph.Ka;
                color.y += li.Ig * sph.Kd * cosTheta * sph.g ;
                //ambient.Ig * sph.Ka;
                color.z += li.Ib * sph.Kd * cosTheta * sph.b ;
                //ambient.Ib * sph.Ka;
            }
            
        }
        
    //}
    
    
    /*
    vec4 tempColor = reflection(ray, t, sph, 0);
    //cout << tempColor <<endl;
    color.x += tempColor.x * sph.Kr ;
    color.y += tempColor.y * sph.Kr ;
    color.z += tempColor.z * sph.Kr ;*/
    
    
    return color;
}






/*
// the temperary struct
struct tempShperes{
    float t;
    string name;
    vec3 S_prime,c_prime;
    vec3 S,c;
}t_sphere;*/


vec4 trace(const Ray& ray)
{
    // TODO: implement your ray tracing routine here.
    //the color
    vec4 color = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    
    sphere sph;
    float t;

     vec4 temp_color = vec4(0,0,0,1);
    
    if (rayIntersectWithSphere(sph, ray, t,true)) {
        //cout << t << endl;
        //cout << sph.name << endl;
        
        vector<light> ::iterator itl;
        for (itl = lights.begin(); itl != lights.end(); itl++)
        {
            if (!addShadow(t, ray, *itl)) {
            
                temp_color += addIllumination(t, ray, sph ,*itl);
            }
        }
        
        temp_color += reflection(ray, t, sph, 0) * sph.Kr;
        
        
        //cout << "tmepcolor: "<<temp_color.x <<" "<<temp_color.y<<" "<<temp_color.z<<endl;
        color.x +=  temp_color.x + ambient.Ir * sph.Ka * sph.r;
        color.y +=  temp_color.y + ambient.Ig * sph.Ka * sph.g;
        color.z +=  temp_color.z + ambient.Ib * sph.Ka * sph.b;
        
        /*
        color.x += sph.r;
        color.y += sph.g;
        color.z += sph.b;
        */
        //cout << "color: "<<color.x <<" "<<color.y<<" "<<color.z<<endl;
    
    }
    else
    {
        color.x = back.r;color.y = back.g; color.z = back.b;
    }

    return color;
}

vec4 getDir(int ix, int iy)
{
    // TODO: modify this. This should return the direction from the origin
    // to pixel (ix, iy), normalized.
    vec4 dir;
    //
    dir.x = g_left   + (ix/float(g_width))*(g_right - g_left);
    dir.y = g_bottom + (iy/float(g_height))*(g_top - g_bottom);
    dir.z = -1 * g_near;
    //dir = vec4(0.0f, 0.0f, -1.0f, 0.0f);
    dir.w = 0.0f;
    //cout <<dir.x <<" "<<dir.y <<" "<<dir.z <<" "<<endl;
    
    //dir = normalize(dir);
    
    return dir;
}

void renderPixel(int ix, int iy)
{
    Ray ray;
    ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    ray.dir = getDir(ix, iy);
    vec4 color = trace(ray);//trace function returning the color in the pic
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

//to clamp the color, I have to do this first
void clamp(vector<vec4> &input){
    vector<vec4>:: iterator i;
    
    for (i = input.begin(); i != input.end(); i++) {
        if (i->x > 1.0f) i->x = 1.0;
        if (i->y > 1.0f) i->y = 1.0;
        if (i->z > 1.0f) i->z = 1.0;
    }
}

void saveFile()
{
    // Convert color components from floats to unsigned chars.
    // TODO: clamp values if out of range.
    
    clamp(g_colors);
    
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
            for (int i = 0; i < 3; i++){
                buf[y*g_width*3+x*3+i] = (unsigned char)(((float*)g_colors[y*g_width+x])[i] * 255.9f);
                //cout <<g_colors[y*g_width+x][i]<< " " <<endl;
            }
    
    // TODO: change file name based on input file name.
    
    string toCertainDirectory;
    string directory = "/Users/xing/Documents/cs174/Assignment3/Assignment3/testFilesResults/";
    
    toCertainDirectory = directory + outputFileName;
    
    //cout << toCertainDirectory << endl;
    
    
    char *cstr = new char [toCertainDirectory.length() + 1];
    strcpy(cstr, toCertainDirectory.c_str());
    
    //char *cstr = new char [outputFileName.length() + 1];
    //strcpy(cstr, outputFileName.c_str());
    savePPM(g_width, g_height, cstr, buf);
    
    
    delete[] buf;
    delete [] cstr;
}



//*********************************
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
    //testForDebug();
    render();
    saveFile();
    
    

	return 0;
}

