/*
    * Abudureheman Adila
    * ECS 175 Project #2 
    * UC Davis, Fall 2019
    */

#ifdef WIN32
#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include "tk.h"
#endif

#if defined(__APPLE__) || defined(MACOSX)
#include <OpenGL/gl.h>
//#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#else //linux
#include <GL/gl.h>
#include <GL/glut.h>
#endif

//other includes
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
//other files


typedef int OutCode;
constexpr int INSIDE = 0; // 0000
constexpr int LEFT = 1;   // 0001
constexpr int RIGHT = 2;  // 0010
constexpr int BOTTOM = 4; // 0100
constexpr int TOP = 8;    // 1000
/****set in main()****/
//the number of pixels in the grid
char *inputFileName;
int grid_width;
int grid_height;

float xMin;
float xMax;
float yMin;
float yMax;
const double pi = 3.14159265358979323846;


//the size of pixels sets the inital window height and width
//don't make the pixels too large or the screen size will be larger than
//your display size
float pixel_size;

/*Window information*/
int win_height;
int win_width;

void init();
void idle();
void display();
void draw_pix(int x, int y);
void draw();
void reshape(int width, int height);
void key(unsigned char ch, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void check();


class Coord {
	public:
		float x, y, z;
    public:
    	Coord(){};
    	Coord(float x, float y, float z) {
    		this->x = x;
    		this->y = y;
            this->z = z;
    	}
};
class ECoord{
	public:
		int x, y;
    public:
    	ECoord(){};
    	ECoord(int x, int y) {
    		this->x = x;
    		this->y = y;
    	}
};


class Polygon {
    //Center of mass is 0 respect to polygon itself
    //But position vector is the centroid from the viewPort
    public:
		int count;
		std::vector<float> position;
        int edgeCount;
		Coord transVec;
		float angle;
		float scale;
		std::vector<std::vector<float>> vertices;
        std::vector<std::vector<int>> edges;

    public:
    	void printPolygon();
    	void updateCentroid();
    	Polygon(){};
    	
    	Polygon(std::vector<Coord> &vert, std::vector <ECoord> &edges) {
		    float xtotal = 0, ytotal = 0;
		    this -> count = vert.size(); 
            this ->edgeCount = edges.size();
		    for(int i = 0; i < vert.size(); i++){
		        vertices.push_back(std::vector<float>{vert[i].x, vert[i].y, vert[i].z});
		    }
            for(int j=0; j<edges.size();j++){
                this->edges.push_back(std::vector<int>{edges[j].x, edges[j].y});
            }
		    this -> updateCentroid();

		    transVec.x = 0.0;
		    transVec.y = 0.0;
            transVec.z = 0.0;
		    angle = 0.0;
		    scale = 1.0;
		}
};

void Polygon::printPolygon() {
        std::cout<<"numVertices: "<<this->count<<std::endl;
        for(int i = 0; i < this->count; i++) {
            std::cout<<"x: "<<(this->vertices[i])[0]<<" y: "<<(this->vertices[i])[1]<<" z: "<<(this->vertices[i])[2]<<std::endl;
        }
        //  for(int j = 0; j < this->edgeCount; j++) {
        //     std::cout<<"Edges "<<(this->edges[j])[0]<<" and "<<(this->edges[j])[1]<<std::endl;
        // }
        std::cout<<std::endl;
    }

void Polygon::updateCentroid() {
        float xtotal = 0, ytotal = 0, ztotal=0;
        for(int i = 0; i < this->count; i++){
            xtotal += this->vertices[i][0];
            ytotal += this->vertices[i][1];
            ztotal += this->vertices[i][2];        
        }
        this->position = {xtotal/(float)(this->count), ytotal/(float)(this->count),  ztotal/(float)(this->count)};
}

void readinput(char *filename, std::vector<Polygon> &polygons);
void writeFile(char *filename, std::vector<Polygon> &polygons);
void drawLineDDA(std::vector<float> start, std::vector<float> end);
void drawLineBresenham(Coord start, Coord end);
//std::vector<float> ComputeIntersection(std::vector<float> a, std::vector<float> b);
//OutCode computeoutbound(std::vector<float> point);
void rasterization(Polygon &p);
bool* loadBuffer;
std::vector<Polygon> polygonList;
std::vector<Polygon> cPolygonList;
void translation(Coord transl, Polygon &poly);
void rotation(float angle, Polygon &poly);
void scaling(float scal, Polygon &poly);
char lineMode;
bool rasterswitch;
float angleG;
float scaleG;
int iD;
float translationXG, translationYG,
sFactor, cliponeX,cliponeY, cliptwoX, cliptwoY;
bool sortVert(const std::vector<float> &a, const std::vector<float> &b);
void polyClip(Polygon &poly);
void polyClipLeft(Polygon &poly);
void polyClipRight(Polygon &poly); 
void polyClipBottom(Polygon &poly);
void polyClipTop(Polygon &poly);
void exchangeV(std::vector<float> &vA, std::vector<float> &vB);  
void copyList(std::vector<std::vector<float>> &s, std::vector<std::vector<float>> &t, int n);
void copyVertex(std::vector<float> &s, std::vector<float> &t);
int main(int argc, char **argv)
{
    inputFileName = "bunny_Scene.txt";
    pixel_size = .5;

    /*Window information*/
    // win_height = grid_height * pixel_size;
    // win_width = grid_width * pixel_size;

    /*Set up glut functions*/
    /** See https://www.opengl.org/resources/libraries/glut/spec3/spec3.html ***/
    
    float translationX=0, translationY=0 , sFactor=1, cliponeX=0,cliponeY=0, cliptwoX=0, cliptwoY=0;
    grid_width = 1000;
    grid_height = 1000;

    xMin = 0;
    xMax = grid_width;
    yMin = 0;
    yMax = grid_height;

    lineMode = 'd';
    rasterswitch = false;
    angleG = 0;
    scaleG = 1;
    iD = -1;
    translationXG = 0;
    translationYG = 0;
    win_height = grid_height * pixel_size;
    win_width = grid_width * pixel_size;
    loadBuffer = new bool[grid_height* grid_width];
    for(int i = 0; i < grid_width; i++){
       for(int j=0; j < grid_height; j++){
           loadBuffer[j*grid_width + i] = false;
       }
    }
    readinput(inputFileName, polygonList);
    polygonList[0].printPolygon();
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    /*initialize variables, allocate memory, create buffers, etc. */
    //create window of size (win_width x win_height
    glutInitWindowSize(win_width, win_height);
    //windown title is "glut demo"
    glutCreateWindow("Two Dimentional Drawing");
   
    /*defined glut callback functions*/
    glutDisplayFunc(display); //rendering calls here
    glutReshapeFunc(reshape); //update GL on window size change
    glutMouseFunc(mouse);     //mouse button events
    glutMotionFunc(motion);   //mouse movement events
    glutKeyboardFunc(key);    //Keyboard events
    glutIdleFunc(idle);       //Function called while program is sitting "idle"

    //initialize opengl variables
    init();
    //start glut event loop
    glutMainLoop();
    return 0;
}

/*initialize gl stufff*/
void init()
{
    //set clear color (Default background to white)
    glClearColor(1.0, 1.0, 1.0, 1.0);
    //checks for OpenGL errors
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	glOrtho(0, grid_width, 0, grid_height, -1, 1);
    check();
}

//called repeatedly when glut isn't doing anything else
void idle()
{
    //redraw the scene over and over again
    glutPostRedisplay();
    iD=0;
    angleG=0;
    translationXG = 0; 
    translationYG=0;
    
    // int choice;
    //     char line;
    //     char ra;
    //     std::cout << "Bresemham or DDA [b/d]: "; 
    //     std::cin>>line;
    //     lineMode = line;
    //     std::cout << "Rasterize all polygons? [y/n]: "; 
    //     std::cin>>ra;
    //     if(ra=='y'){rasterswitch=true;}else{rasterswitch=false;}
    //     std::cout << "1. Rotation \n";  
    //     std::cout << "2. Translation\n";  
    //     std::cout << "3. Scalling \n";  
    //     std::cout << "4. Clipping \n";
    //     std::cout << "5. Exit \n";
    //     std::cout << "Please select one of options above for your operation: ";
    //     std::cin>> choice;
    //     switch (choice) 
    //     { 
    //         case 1:  
    //             std::cout << "Please enter rotation angle: ";
    //             std::cin>> angleG;
    //             std::cout << "Please enter Polygon ID, -1 for all together, 1,2,3.. for specific ones: ";
    //             std::cin>> iD;
    //             break;
    //         case 2:
    //             std::cout << "Please enter translation in x and y direction: ";
    //             std::cin>> translationXG >> translationYG;
    //             std::cout << "Please enter Polygon ID, -1 for all together, 1,2,3.. for specific ones: ";
    //             std::cin>> iD;
    //             break; 
    //         case 3:  
    //             std::cout << "Please enter scalling factor: " ;
    //             std::cin>> scaleG;
    //             std::cout << "Please enter Polygon ID, -1 for all together, 0,1,2, for specific ones: ";
    //             std::cin>> iD;
    //             break; 
    //         case 4:  
    //             std::cout << "Please enter left clipping coordinates seperated by a space: ";
    //             std::cin>> xMin >> yMin;
    //             std::cout << "Please enter the right clipping coordinates seperated by a space: " ;
    //             std::cin>> xMax >> xMax;
    //             break; 
    //         case 5: 
    //             writeFile(inputFileName, polygonList);
    //             exit(0);
    //             break; 
    //         default:  
    //             break;
    //     }
}

void swapCor(Coord start, Coord end)
{
    int tempx = start.x;
    int tempy = start.y;
    start.x = end.x;
    start.y = end.y;
    end.x = tempx;
    end.y = tempy;
}
void readinput(char *filename, std::vector<Polygon> &polygons){
    std::ifstream inputFile;
    inputFile.open(filename);
    std::string line;
    int count;
    int edgeCount;
    inputFile >> count;
    getline(inputFile, line); //empty
    getline(inputFile, line); //point count
    for (int i=0; i< count; i++){
        int num;
        std::vector <Coord> vertices;
        std::vector <ECoord> edges;
        inputFile >> num;
        getline(inputFile, line);
        for (int j=0; j<num; j++){
            float x, y, z;
            std::string inputX, inputY, inputZ;
            getline(inputFile, line);
            std :: istringstream record(line);
            getline(record, inputX, ' ');
            getline(record, inputY, ' ');
            getline(record, inputZ);
            x = std::stof(inputX);
            y = std::stof(inputY);
            z = std::stof(inputZ);
            Coord point;
            point.x = x;
            point.y = y;
            point.z = z;
            vertices.push_back(point);
        }
        std :: string input;
        getline(inputFile, line);
        edgeCount = std::stoi(line);
        for (int j=0; j<edgeCount; j++){
            int x, y;
            std :: string inputX, inputY;
            getline(inputFile, line);
            std :: istringstream record(line);
            getline(record, inputX, ' ');
            getline(record, inputY);
            x = std::stoi(inputX);
            y = std::stoi(inputY);
            ECoord point;
            point.x = x;
            point.y = y;
            edges.push_back(point);
        }
        Polygon polygon(vertices, edges);
        polygons.push_back(polygon);
        getline(inputFile, line);
        if(line == ""){
            getline(inputFile, line);
        }
    }
    inputFile.close();
}
void writeFile(char *filename,std::vector<Polygon> &polygons){
    std::ofstream outputFile(filename);
    outputFile << polygons.size() << "\n\n";
    for (int i = 0; i<polygons.size();i++) {
        outputFile << polygons[i].count << std::endl;
        for(int j = 0; j<polygons[i].count;j++) {
            outputFile << polygons[i].vertices[j][0] << ' ' << polygons[i].vertices[j][1] << polygons[i].vertices[j][2]<< std::endl;
        }
        outputFile << polygons[i].edgeCount<< std::endl;
        for(int k = 0; k<polygons[i].count;k++) {
            outputFile << polygons[i].edges[k][0] << ' ' << polygons[i].edges[k][1] << polygons[i].edges[k][2]<< std::endl;
        }
        outputFile << std::endl;
    }
}

void drawLineDDA(std::vector<float> start, std::vector<float> end)
{
    int startX = (int)(start[0] + 0.5);
    int startY = (int)(start[1] + 0.5);
    int endX = (int)(end[0] + 0.5);
    int endY = (int)(end[1] + 0.5);
    
    int dX = endX - startX;
    int dY = endY - startY;
    int steps, k;
    float incX, incY;
    float x = (float)startX;
    float y = (float)startY;
    
    if(fabs(dX) > fabs(dY)) {
    	steps = fabs(dX);
    } else {
    	steps = fabs(dY);
    }
    incX = (float)dX / (float)steps;
    incY = (float)dY / (float)steps;
    
    loadBuffer[(int)(x + grid_width * y + 0.5)] = true;
    for(int i = 0; i < steps; i++) {
    	x += incX;
    	y += incY;
    	int roundX = (int)(x+0.5);
    	int roundY = (int)(y+0.5);
    	loadBuffer[roundX + grid_width * roundY] = true;
    }
}



//Algorithm from class notes & textbook 
void drawLineBresenham(std::vector<float> start, std::vector<float> end)
{
    float m = (end[1] - start[1]) / (end[0] - start[0]);
    int x = 0, y = 0;
    if(m == 1){
        if(start[0]<end[0]){
            y = start[0];
        }else{
            y = end[0];
        }
        for (int x = fmin(round(start[0]), round(end[0])); x <= fmax(round(start[0]), round(end[0])); x++) {
            //Coord point(x, y++);
            loadBuffer[y*grid_width+x]=true;
        }
    }else if(m == -1){
        int y;
        for (int x = fmin(round(start[0]), round(end[0])), y = fmax(round(start[1]), round(end[1])); x <= fmax(round(start[0]), round(end[0])); x++, y--) {
            //Coord point(x, y);
            loadBuffer[y*grid_width+x]=true;
        }

    }
    if (fabs(m) < 1) {
        int dx = fabs(end[0] - start[0]),
            dy = fabs(end[1] - start[1]),
            p = 2 * dy - dx;
        if (start[0] > end[0])
        {
            x = end[0];
            y = end[1];
            end[0] = start[0];
        }else{
            x = start[0];
            y = start[1];
        }
        while (x < end[0])
        {
            x++;
            if (p < 0)
            {
                p = p + 2 * dy;
            }
            else
            {
                if (m> 0) {y++;} else {y--;}
                p = p + 2 * dy - 2 * dx;
            }
            //draw_pix(x, y);
            loadBuffer[y*grid_width+x]=true;
        }
    } else if (fabs(m) >= 1) {                    
        int dx = fabs(end[0] - start[0]),
            dy = fabs(end[1] - start[1]),
            p = 2 * dx - dy;
        if (start[1] > end[1])
        {
            x = end[0];
            y = end[1];
            end[1] = start[1];
        }
        else
        {
            x = start[0];
            y = start[1];
        }
        while (y < end[1])
        {
            y++;
            if (p < 0)
            {
                p = p + (2 * dx);
            }else{
                if(m>0){
                    x++;
                }else{
                    x--;
                }
                p = p + (2 * dx) - (2 * dy);
            }
             //draw_pix(x, y);
             loadBuffer[y*grid_width+x]=true;
        }
    } 
}

bool sortVert(const std::vector<float> &a, const std::vector<float> &b) {
	return (a[0] < b[0]);
}
   
void rasterization(Polygon &p)
{
	int numLine = p.count;
	for(int y = 0; y < grid_height; y++) {
		std::vector<std::vector<float>> arr;
		//for each line of poly
		for(int n = 0; n < numLine; n++) {
			std::vector<float> vA = p.vertices[n];
			std::vector<float> vB = p.vertices[(n + 1) % numLine];
			
			float x1 = vA[0];
			float y1 = vA[1];
			float x2 = vB[0];
			float y2 = vB[1];		
			
			if(y1 == y && y2 == y) {
				continue;
			}
			
			if((y1 > y && y2 > y) || (y1 < y && y2 < y)) {
				continue;
			}
			
			//swap y2 >> y1
			if(y2 < y1) {
				float tempx = x1;
				float tempy = y1;
				x1 = x2;
				y1 = y2;
				x2 = tempx;
				y2 = tempy;
			}
			
			//min
			if(y1 == y) {
				arr.push_back(std::vector<float>{x1, (float)y, 1.0});
				continue;
			}
			
			//max
			if(y2 == y) {
				continue;
			}
			
			//vertical
			if(x1 == x2) {
				arr.push_back(std::vector<float>{x1, (float)y, 1.0});
				continue;
			}
			
			//intx
			float k = (y2-y1) / (x2-x1);
            float intersectVal = x1 + (y-y1) / k;
            arr.push_back(std::vector<float>{intersectVal, (float)y, 1.0});
		}
				
		std::sort(arr.begin(), arr.end(), sortVert);
		
		int numPair = arr.size() / 2;
		for(int k = 0; k < numPair; k++) {
			std::vector<float> p1 = arr[2 * k];
			std::vector<float> p2 = arr[2 * k + 1];
			
			float xL = p1[0];
			float xH = p2[0];
			
			if(floor(xH) == xH) {
				xH--;
			}
			
			int xLi = ceil(xL);
			int xHi = floor(xH);
			
			for(int posX = xLi; posX < xHi + 1; posX++) {
				loadBuffer[posX + y * grid_width] = true;
			}
		}
	}
}
void translation(Coord transl, Polygon &poly){
    
    for(int i = 0; i < poly.count; i++) {
        std::vector<float> temp = poly.vertices[i];
        temp[0] = temp[0] + transl.x;
        temp[1] = temp[1] + transl.y;
        temp[2] = temp[2] + transl.z;
        poly.vertices[i] = temp;
    } 
    poly.updateCentroid();
}

//需要改～
void rotation(float angle, Polygon &poly){
    
    Coord trans2Ori;
    trans2Ori.x = -poly.position[0];
    trans2Ori.y = -poly.position[1];
    trans2Ori.z = -poly.position[2];
    Coord trans2Cen;
    trans2Cen.x = poly.position[0]; 
    trans2Cen.y = poly.position[1];
    trans2Cen.z = poly.position[2];
    float rotateAngleInRad = (angle / 180.0) * 3.14;
    float cosA = cos(rotateAngleInRad);
    float sinA = sin(rotateAngleInRad);

    translation(trans2Ori, poly);   
    for(int i = 0; i < poly.count; i++) {
        float currX = poly.vertices[i][0];
        float currY = poly.vertices[i][1];
        float currZ = poly.vertices[i][2];
        float nextX = currX * cosA + currY * (-sinA);
        float nextY = currX * sinA + currY * cosA;
        
        std::vector<float> temp = {nextX,nextY,1} ;
        poly.vertices[i] = temp;     
    }
    translation(trans2Cen, poly);    
}


void scaling(float scal, Polygon &poly){
    Coord trans2Ori; 
    trans2Ori.x = -poly.position[0];
    trans2Ori.y = -poly.position[1];
    trans2Ori.z = -poly.position[2];
    Coord trans2Cen;
    trans2Cen.x = poly.position[0];
    trans2Cen.y = poly.position[1];
    trans2Cen.z = poly.position[2];
    translation(trans2Ori, poly);
    for(int i = 0; i < poly.count; i++) {
        float currX = poly.vertices[i][0];
        float currY = poly.vertices[i][1];
        float currZ = poly.vertices[i][2];
        //std::cout<<"currx: "<<currX<<" currY: "<<currY<<std::endl;
        float nextX = currX * scal;
        float nextY = currY * scal;
        float nextZ = currZ * scal;
        //std::cout<<"nextx: "<<nextX<<" nextY: "<<nextY<<std::endl;

        poly.vertices[i][0] = nextX;
        poly.vertices[i][1] = nextY;
        poly.vertices[i][2] = nextZ;
    }
    translation(trans2Cen, poly);  
}

/////////////////////////////////////////////////////////////////////////////////////////////
typedef float Matrix4x4 [4][4];
Matrix4x4 matComposite;
void matrix4x4SetIdentity (Matrix4x4 matIdent4x4)
{
   int row, col;
   for (row = 0; row < 4; row++)
      for (col = 0; col < 4 ; col++)
        matIdent4x4 [row][col] = (row == col);
}

/* Premultiply matrix m1 by matrix m2, store result in m2. */
void matrix4x4PreMultiply (Matrix4x4 m1, Matrix4x4 m2)
{
   int row, col;
   Matrix4x4 matTemp;
   for (row = 0; row < 4; row++)
      for (col = 0; col < 4 ; col++)
         matTemp [row][col] = m1 [row][0] * m2 [0][col] + m1 [row][1] *
                            m2 [1][col] + m1 [row][2] * m2 [2][col] +
                            m1 [row][3] * m2 [3][col];
   for (row = 0; row < 4; row++)
      for (col = 0; col < 4; col++)
         m2 [row][col] = matTemp [row][col];
}
/*  Procedure for generating 3-D translation matrix.  */
void translate3D (float tx, float ty, float tz)
{
   Matrix4x4 matTransl3D;
   /*  Initialize translation matrix to identity.  */
   matrix4x4SetIdentity (matTransl3D);
   matTransl3D [0][3] = tx;
   matTransl3D [1][3] = ty;
   matTransl3D [2][3] = tz;
   /*  Concatenate matTransl3D with composite matrix.  */
   matrix4x4PreMultiply (matTransl3D, matComposite);
}
/*  Procedure for generating a quaternion rotation matrix.  */
void rotate3D (Coord p1, Coord p2, float angle)
{
   Matrix4x4 matQuatRot;
   float radianAngle= (angle / 180.0) * 3.14;
   float axisVectLength = sqrt ((p2.x - p1.x) * (p2.x - p1.x) +
                        (p2.y - p1.y) * (p2.y - p1.y) +
                        (p2.z - p1.z) * (p2.z - p1.z));
   float cosA = cosf (radianAngle);
   float oneC = 1 - cosA;
   float sinA = sinf (radianAngle);
   float ux = (p2.x - p1.x) / axisVectLength;
   float uy = (p2.y - p1.y) / axisVectLength;
   float uz = (p2.z - p1.z) / axisVectLength;
   /*  Set up translation matrix for moving p1 to origin,
    *  and concatenate translation matrix with matComposite.
    */
   translate3D (-p1.x, -p1.y, -p1.z);
   /*  Initialize matQuatRot to identity matrix.  */
   matrix4x4SetIdentity (matQuatRot);
   matQuatRot [0][0] = ux*ux*oneC + cosA;
   matQuatRot [0][1] = ux*uy*oneC - uz*sinA;
   matQuatRot [0][2] = ux*uz*oneC + uy*sinA;
   matQuatRot [1][0] = uy*ux*oneC + uz*sinA;
   matQuatRot [1][1] = uy*uy*oneC + cosA;
   matQuatRot [1][2] = uy*uz*oneC - ux*sinA;
   matQuatRot [2][0] = uz*ux*oneC - uy*sinA;
   matQuatRot [2][1] = uz*uy*oneC + ux*sinA;
   matQuatRot [2][2] = uz*uz*oneC + cosA;
   /*  Concatenate matQuatRot with composite matrix.  */
   matrix4x4PreMultiply (matQuatRot, matComposite);
   /*  Construct inverse translation matrix for p1 and
    *  concatenate with composite matrix.
    */
   translate3D (p1.x, p1.y, p1.z);
}
/*  Procedure for generating a 3-D scaling matrix.  */
void scale3D (float sx, float sy, float sz, std::vector<float>fixedPt)
{
   Matrix4x4 matScale3D;
   /*  Initialize scaling matrix to identity.  */
   matrix4x4SetIdentity (matScale3D);
   matScale3D [0][0] = sx;
   matScale3D [0][3] = (1 - sx) * fixedPt[0]*grid_width;
   matScale3D [1][1] = sy;
   matScale3D [1][3] = (1 - sy) * fixedPt[1]*grid_width;
   matScale3D [2][2] = sz;
   matScale3D [2][3] = (1 - sz) * fixedPt[2]*grid_width;
   /*  Concatenate matScale3D with composite matrix.  */
   matrix4x4PreMultiply (matScale3D, matComposite);
}
void drawLine(float x1 ,float y1, float x2, float y2){
    glBegin(GL_LINES);
    glVertex2f(x1,y1);
    glVertex2f(x2,y2);
    glEnd();
}
/////////////////////////////////////////////////////////////////////////////////////////////
//this is where we render the screen
void display()
{
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    glLoadIdentity();
    //Translate my NDC to WC
    //std::cout<<"X: "<<polygonList[0].vertices[0][0]<<"  Y: "<<polygonList[0].vertices[0][1]<<std::endl;
    //std::cout<<"X: "<<polygonList[0].vertices[0][1]<<"  Y: "<<polygonList[0].vertices[0][1]<<std::endl;
    //define my matrix with my coordinates
    //gltranslatef
    //glscale
    matrix4x4SetIdentity (matComposite);
    rotate3D (Coord((polygonList[0].position[0]*100),(polygonList[0].position[2]*100),0), Coord((polygonList[0].position[0]*100),(polygonList[0].position[2]*100)+5,0), 90);  //  First transformation: Rotate.
    scale3D (8.0, 8.0, 8.0, polygonList[0].position);   //  Second transformation: Scale.
    //translate3D (20, 20, 20);        //  Final transformation: Translate.
    //polygonList[0].printPolygon();
    cPolygonList.clear();
    for(int u=0;u<polygonList.size();u++){
        cPolygonList.push_back(polygonList[u]);
    }
    for(int i=0; i<polygonList[0].count; i++){
        cPolygonList[0].vertices[i][0] = matComposite[0][0]*polygonList[0].vertices[i][0]*grid_width +matComposite[0][1]* polygonList[0].vertices[i][1]*grid_width + matComposite[0][2]*polygonList[0].vertices[i][2]*grid_width+matComposite[0][3]*1;
        cPolygonList[0].vertices[i][1] = matComposite[1][0]*polygonList[0].vertices[i][0]*grid_width +matComposite[1][1]* polygonList[0].vertices[i][1]*grid_width+ matComposite[1][2]*polygonList[0].vertices[i][2]*grid_width+matComposite[1][3]*1; 
        cPolygonList[0].vertices[i][2] = matComposite[2][0]*polygonList[0].vertices[i][0]*grid_width +matComposite[2][1]* polygonList[0].vertices[i][1]*grid_width + matComposite[2][2]*polygonList[0].vertices[i][2]*grid_width+matComposite[2][3]*1;
    }
    //std::cout<<"a: "<<matComposite[0][3]<<"  b: "<<matComposite[1][3]<<" c: "<<matComposite[2][3]<<std::endl;
    for(int i = 0; i<polygonList[0].count; i++){
        draw_pix((cPolygonList[0].vertices[i][0]),(cPolygonList[0].vertices[i][1]));
        //std::cout<<"X: "<<(polygonList[0].vertices[i][0])*100<<"  Y: "<<(polygonList[0].vertices[i][1])*100<<std::endl;
    }
    for(int k = 0; k<cPolygonList[0].edgeCount;k++){
        int a = cPolygonList[0].edges[k][0]-1;
        int b = cPolygonList[0].edges[k][1]-1;
        drawLine(cPolygonList[0].vertices[a][0],cPolygonList[0].vertices[a][1],cPolygonList[0].vertices[b][0],cPolygonList[0].vertices[b][1]);
    }
    glutSwapBuffers();
    //checks for opengl errors
    check();
    
    /*
	cPolygonList.clear();
	delete(loadBuffer);
    loadBuffer = NULL;
	loadBuffer = new bool[grid_width * grid_height];
	for(int i = 0; i < grid_width * grid_height; i++) {
		loadBuffer[i] = false;
	}
	
    //clears the screen
    glClear(GL_COLOR_BUFFER_BIT);

    // make temp geometry 
    for(int u=0;u<polygonList.size();u++){
        cPolygonList.push_back(polygonList[u]);
    }
    if(iD == -1){
        for(int u=0;u<cPolygonList.size();u++){
        cPolygonList[u].angle = angleG;
        cPolygonList[u].scale = scaleG;
        cPolygonList[u].transVec.x = translationXG;
        cPolygonList[u].transVec.y = translationYG;

        if(!(cPolygonList[u].transVec.x == 0 && cPolygonList[u].transVec.y == 0 )) {
            translation(cPolygonList[u].transVec, cPolygonList[u]);
        }

        if(!(cPolygonList[u].angle == 0)) {
            rotation(cPolygonList[u].angle, cPolygonList[u]);
        }

        if(!(cPolygonList[u].scale == 1)) {
            scaling(cPolygonList[u].scale, cPolygonList[u]);
        }

        cPolygonList[u].transVec = Coord(0.0, 0.0, 0.0);
        cPolygonList[u].angle = 0.0;
        cPolygonList[u].scale = 1.0;
        }
    }else{
        cPolygonList[iD].angle = angleG;
        cPolygonList[iD].scale = scaleG;
        cPolygonList[iD].transVec.x = translationXG;
        cPolygonList[iD].transVec.y = translationYG;

        if(!(cPolygonList[iD].transVec.x == 0 && cPolygonList[iD].transVec.y == 0 )) {
            translation(cPolygonList[iD].transVec, cPolygonList[iD]);
        }

        if(!(cPolygonList[iD].angle == 0)) {
            rotation(cPolygonList[iD].angle, cPolygonList[iD]);
        }

        if(!(cPolygonList[iD].scale == 1)) {
            scaling(cPolygonList[iD].scale, cPolygonList[iD]);
        }
        cPolygonList[iD].transVec = Coord(0.0, 0.0, 0.0);
        cPolygonList[iD].angle = 0.0;
        cPolygonList[iD].scale = 1.0;
    }
    
    
    //copy cpoly to poly then work from there next round
    polygonList.clear();
    for(int u=0;u<cPolygonList.size();u++){
        polygonList.push_back(cPolygonList[u]);
        
    }
    for(int a = 0; a<cPolygonList.size(); a++){
        polyClip(cPolygonList[a]);
    }

    for(auto p : cPolygonList){
        for(int i = 0; i<p.count;i++){
            std::vector<float> cur = p.vertices[i];
            std::vector<float> prev = p.vertices[(i + p.count - 1) % p.count];
            
           	if(lineMode == 'd'){
                drawLineDDA(cur, prev);
           	}else{
                drawLineBresenham(cur, prev);
           	}
        }
        
        if(rasterswitch){
        	for(auto p : cPolygonList) {
        		rasterization(p);
        	}
        }
    }    
    
    // draw();
    delete(loadBuffer);
    loadBuffer = NULL;*/
}

//Draws a single "pixel" given the current grid size
//don't change anything in this for project 1
void draw_pix(int x, int y)
{
    glBegin(GL_POINTS);
    glColor3f(0.6, 0.5, 0.0);
    glVertex3f(x + .5, y + .5, 0);
    glEnd();
}

void draw() {
	glColor3f(.2, .2, 1.0);
	
	for(int i = 0; i < grid_width * grid_height; i++) {
		if(loadBuffer[i]) {
			glBegin(GL_POINTS);
			glVertex2i(i % grid_width, i / grid_width);
			glEnd();
		}
	}	
	glFlush();
}

/*Gets called when display size changes, including initial craetion of the display*/
void reshape(int width, int height)
{
    /*set up projection matrix to define the view port*/
    //update the ne window width and height
    win_width = width;
    win_height = height;

    //creates a rendering area across the window
    glViewport(0, 0, width, height);
    // up an orthogonal projection matrix so that
    // the pixel space is mapped to the grid space
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, grid_width, 0, grid_height, -10, 10);

    //clear the modelview matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    //set pixel size based on width, if the aspect ratio
    //changes this hack won't work as well
    pixel_size = width / (float)grid_width;

    //set pixel size relative to the grid cell size
    glPointSize(pixel_size);
    //check for opengl errors
    check();
}

//gets called when a key is pressed on the keyboard
void key(unsigned char ch, int x, int y)
{
    switch (ch)
    {
    case 'b':
        lineMode = 'b';
        break;

    case 'd':
        lineMode = 'd';
        break;

    case 'r':
        rasterswitch = ! rasterswitch;
        break;
    
    default:
        //prints out which key the user hit
        printf("User hit the \"%c\" key\n", ch);
        break;
    }
    //redraw the scene after keyboard input
    glutPostRedisplay();
}


//gets called when a mouse button is pressed
void mouse(int button, int state, int x, int y)
{
    //print the pixel location, and the grid location
    printf("MOUSE AT PIXEL: %d %d, GRID: %d %d\n", x, y, (int)(x / pixel_size), (int)((win_height - y) / pixel_size));
    switch (button)
    {
    case GLUT_LEFT_BUTTON: //left button
        printf("LEFT ");
        break;
    case GLUT_RIGHT_BUTTON: //right button
        printf("RIGHT ");
    default:
        printf("UNKNOWN "); //any other mouse button
        break;
    }
    if (state != GLUT_DOWN) //button released
        printf("BUTTON UP\n");
    else
        printf("BUTTON DOWN\n"); //button clicked

    //redraw the scene after mouse click
    glutPostRedisplay();
}

//gets called when the curser moves accross the scene
void motion(int x, int y)
{
    //redraw the scene after mouse movement
    glutPostRedisplay();
}

//checks for any opengl errors in the previous calls and
//outputs if they are present
void check()
{
    GLenum err = glGetError();
    if (err != GL_NO_ERROR)
    {
        printf("GLERROR: There was an error %s\n", "error");
        exit(1);
    }
}


void polyClip(Polygon &poly) {
    for(int dir = 0; dir < 4; dir++) {
        if(dir == 0) {
            polyClipLeft(poly);
        } else if(dir == 1) {
            polyClipRight(poly);
        } else  if(dir == 2) {
            polyClipBottom(poly);
        } else if(dir == 3) {
            polyClipTop(poly);
        }
    }
}

void polyClipLeft(Polygon &poly) {
    std::vector<std::vector<float>> vs;
    int n = poly.count;
    int vsLen = 0;

    for(int i = 0; i < n; i++) {
        std::vector<float> vA = poly.vertices[i];
        std::vector<float> vB = poly.vertices[(i + 1) % n];

        //outside?
        if(vA[0] < xMin && vB[0] < xMin) {
            //discard
        } else if(vA[0] >= xMin && vB[0] >= xMin) {
            //inside?
            //keep later one
            vs.push_back(vB);
            vsLen += 1;
        } else  {
            //partial            
            //find vert
            float k = (vB[1]-vA[1]) / (vB[0] - vA[0]);
            float intersectVal = vA[1] + (xMin-vA[0]) * k;
            std::vector<float> v = {xMin, intersectVal, 1.0};
            vs.push_back(v);
            vsLen += 1;  
            
            if(vB[0] > xMin) {
                vs.push_back(vB);   
                vsLen += 1;  
            }      
        }
    }

    poly.count = vsLen;
    poly.vertices = vs;
}

void polyClipRight(Polygon &poly) {
    std::vector<std::vector<float>> vs;
    int n = poly.count;
    int vsLen = 0;

    for(int i = 0; i < n; i++) {
        std::vector<float> vA = poly.vertices[i];
        std::vector<float> vB = poly.vertices[(i + 1) % n];

        //outside?
        if(vA[0] > xMax && vB[0] > xMax) {
            //discard
        } else if(vA[0] <= xMax && vB[0] <= xMax) {
            //inside?
            vs.push_back(vB);
            vsLen += 1;
        } else  {
            //partial
            //find vert
        
            float k = (vB[1]-vA[1]) / (vB[0] - vA[0]);
            float intersectVal = vA[1] + (xMax-vA[0]) * k;
            std::vector<float> v = {xMax, intersectVal, 1.0};

            vs.push_back(v);
            vsLen += 1;      
            
            if(vB[0] < xMax) {
                vs.push_back(vB);
                vsLen += 1;
            }  
        }
    }

    poly.count = vsLen;     
    poly.vertices = vs;
}

void polyClipBottom(Polygon &poly) {
    std::vector<std::vector<float>> vs;
    int n = poly.count;
    int vsLen = 0;

    for(int i = 0; i < n; i++) {
        std::vector<float> vA = poly.vertices[i];
        std::vector<float> vB = poly.vertices[(i + 1) % n];

        //outside?
        if(vA[1] < yMin && vB[1] < yMin) {
            //discard
        } else if(vA[1] >= yMin && vB[1] >= yMin) {
            //inside?
            vs.push_back(vB);
            vsLen += 1;
        } else  {
            //find vert
        
            float k = (vB[1]-vA[1]) / (vB[0] - vA[0]);
            float intersectVal = vA[0] + (yMin-vA[1]) / k;
            std::vector<float> v = {intersectVal, yMin, 1.0};

            vs.push_back(v);
            vsLen += 1;   
                 
            //partial
            if(vB[1] > yMin) {
                vs.push_back(vB);   
                 vsLen += 1; 
            }
        }
    }
    poly.count = vsLen;    
    poly.vertices = vs;
}

void polyClipTop(Polygon &poly) {
    std::vector<std::vector<float>> vs;
    int n = poly.count;
    int vsLen = 0;

    for(int i = 0; i < n; i++) {
        std::vector<float> vA = poly.vertices[i];
        std::vector<float> vB = poly.vertices[(i + 1) % n];

        //outside?
        if(vA[1] > yMax && vB[1] > yMax) {
            //discard
        } else if(vA[1] <= yMax && vB[1] <= yMax) {
            //inside?
            vs.push_back(vB);
            vsLen += 1;
        } else  { 
            //find vert
        
            float k = (vB[1]-vA[1]) / (vB[0] - vA[0]);
            float intersectVal = vA[0] + (yMax-vA[1]) / k;
            std::vector<float> v = {intersectVal, yMax, 1.0};

            vs.push_back(v);
            vsLen += 1;  
            
            //partial
            if(vB[1] < yMax) {
                vs.push_back(vB);   
                vsLen += 1;  
            }
        }
    }

    poly.count = vsLen;    
    poly.vertices = vs;
}

void exchangeV(std::vector<float> &vA, std::vector<float> &vB) {
    float tempX, tempY;
    tempX = vA[0];
    tempY = vA[1];
    vA[0] = vB[0];
    vA[1] = vB[1];
    vB[0] = tempX;
    vB[1] = tempY;
}
