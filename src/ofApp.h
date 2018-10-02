#pragma once

#include "ofMain.h"

class Particle : public ofBaseApp {
public:
    void setup();
    void update();
    void draw();
    void setDebug(bool* debug);
    void setParticles(vector<Particle>* particles);
    bool isInteracting(Particle * particle);
    void updateParticles(Particle * particle2);
	// position (centre of mass)
    ofVec2f p; 
	// velocity
    ofVec2f v;
	// mass
    float m; 
    float radius;
	//physical equations
    float K; //average kinetic energy per gas molecule
	float T; // temperature of a gas molecule
	
	
	bool updated; // has this particle been updated
	bool * debug; // pointer to apps debug flag
    vector<Particle>* particles; // pointer to other particles
};

class ofApp : public ofBaseApp{
public:
    void setup();
	void update();
    void draw();
    void keyPressed(int key);
	
    bool debug; // used for displaying debug information
    float K; // calculation of the total energy kinetic of the system by molecular gas
	float T;//calculation of the total temperature K of the system by molecular gas
	float cv;
	float CV;
	float U;
	float Q;
	float W;
    vector<Particle> particles; // map of all the particles
};

