#include "ofApp.h"
void Particle::setup() {
	radius = 8;
	m = 2 * PI * pow(radius, 2);
	// assign it's starting position
	p.x = ofRandom(radius, ofGetWidth() - radius);
	p.y = ofRandom(radius, ofGetHeight() - radius);
	// assign it's velocity
	v.x = 1;
	v.y =-1;
}
//------------------------------------------------------------------------------------------------------
void Particle::update() {
	float minX = 0;
	float minY = 0;
	float maxX = ofGetWidth();
	float maxY = ofGetHeight();
	ofVec2f reverseX(-1, 1);
	ofVec2f reverseY(1, -1);
	if (p.x > maxX || p.x < minX) {
		v *= reverseX;
	}
	if (p.y > maxY || p.y < minY) {
		v *= reverseY;
	}
	if (updated) {
		return;
	}
	for (int i = 0; i < particles->size(); i++) {
		if (p.x == particles->at(i).p.x &&
			p.y == particles->at(i).p.y) {
			continue;
		}
		Particle * particle2 = &particles->at(i);
		if (isInteracting(particle2)) {
			// Using the following model to update veloctiy and angular velocity:
			// http://www.euclideanspace.com/physics/dynamics/collision/twod/
			// get all the required information:
			// mass
			float m1 = m;
			float m2 = particle2->m;
			// radius
			float radius1 = radius;
			float radius2 = particle2->radius;
			// inertiass
			float i1 = (PI / 2)*pow(radius1, 4);
			float i2 = (PI / 2)*pow(radius2, 4);
			// velocitys
			ofVec2f v1 = v;
			ofVec2f v2 = particle2->v;
			// positions
			ofVec2f p1 = p;
			ofVec2f p2 = particle2->p;
			// relative vector of collision point to centre of mass
			ofVec2f r1 = (p1 - p2).scale(radius1 / radius2);
			ofVec2f r2 = (p2 - p1).scale(radius2 / radius1);
			// Impulse
			ofVec2f j;
			// e is the coefficient of restitution. It's fun to vary this.
			float e = 0.9;
			float k = 1 / (m1*m1) + 2 / (m1*m2) + 1 / (m2*m2) - r1.x*r1.x / (m1*i1) - r2.x*r2.x / (m1*i2) - r1.y*r1.y / (m1*i1)
				- r1.y*r1.y / (m2*i1) - r1.x*r1.x / (m2*i1) - r2.x*r2.x / (m2*i2) - r2.y*r2.y / (m1*i2)
				- r2.y*r2.y / (m2*i2) + r1.y*r1.y*r2.x*r2.x / (i1*i2) + r1.x*r1.x*r2.y*r2.y / (i1*i2) - 2 * r1.x*r1.y*r2.x*r2.y / (i1*i2);
			// set the impulse in the two dimensions
			j.x = (e + 1) / k * (v1.x - v2.x)*(1 / m1 - r1.x*r1.x / i1 + 1 / m2 - r2.x*r2.x / i2)
				- (e + 1) / k * (v1.y - v2.y)* (r1.x*r1.y / i1 + r2.x*r2.y / i2);
			j.y = -(e + 1) / k * (v1.x - v2.x) * (r1.x*r1.y / i1 + r2.x*r2.y / i2)
				+ (e + 1) / k * (v1.y - v2.y) * (1 / m1 - r1.y*r1.y / i1 + 1 / m2 - r2.y*r2.y / i2);
			// velocity and angular velocity after the collision
			ofVec2f v1f = v1 - j / m1;
			ofVec2f v2f = v2 + j / m2;
			// update this particles velocities
			v = v1f;
			particles->at(i).v = v2f;
			particles->at(i).p += particles->at(i).v;
			updated = true;
			particles->at(i).updated;
		}
	}
	// update the position using the particles velocity
	p += v;
	//kinetic energy of average translation of a gas molecule
	K = 0.5 * m * pow(v.length(), 2);// J
}
//------------------------------------------------------------------------------------------------------
void Particle::draw() {
	ofDrawCircle(p.x, p.y, radius);
	if (*debug) {
		ofDrawBitmapString(K, p.x + radius, p.y + radius);
	}
}
//------------------------------------------------------------------------------------------------------
void Particle::setDebug(bool * d) {
	debug = d;
}
//------------------------------------------------------------------------------------------------------
void Particle::setParticles(vector<Particle> *p) {
	particles = p;
}
//------------------------------------------------------------------------------------------------------
bool Particle::isInteracting(Particle * particle) {
	// calculate the distance between two particles
	float d = p.distance(particle->p);
	if (d < radius + particle->radius) {
		return true;
	}
	return false;
}
//------------------------------------------------------------------------------------------------------
void Particle::updateParticles(Particle * particle2){
}
//------------------------------------------------------------------------------------------------------
void ofApp::setup(){
    ofBackground(0);
	ofSetFrameRate(60);
    const int NUM_PARTICLES =250;
    for(int i=0; i<NUM_PARTICLES; i++){
        Particle particle;
        particle.setup();
        particles.push_back(particle);
    }
    for(int i=0; i<particles.size(); i++){
        particles[i].setDebug(&debug);
        particles[i].setParticles(&particles);
    }
}
//------------------------------------------------------------------------------------------------------
void ofApp::update() {
	for (int i = 0; i < particles.size(); i++) {
		particles[i].updated = false;
	}
	for (int i = 0; i < particles.size(); i++) {
		particles[i].update();
	}
	K = 0;
	for (int i = 0; i < particles.size(); i++) {
		K += particles[i].K;
		//temperature for one mole of ideal gas---> T= 2*K/(3nR) n=1
		float R = 8.314472; // J/mol*K
		cv = 1.5*R; //J/mol*K
		T = K/cv; // K \ K= (3/2)RT , 1/cv = 0.7/R
	    //calculation of the internal energy of the system
		U = K; // U = Cv * T-- > U = 3 / 2 * R*T & T = (2 / 3 * R)*K
		// if dU = dW + dQ & dW = -pdV but V = const implies that
		//dV = 0 then dW = 0 implies that dU = dQ
		//U = cv * T;
		//Q = cv * T;
		Q = U;
		W = U - Q;
		CV = Q/T;
	}
}
//------------------------------------------------------------------------------------------------------
void ofApp::draw(){
    for(int i=0; i<particles.size(); i++){
        particles[i].draw();
    }
    ofSetColor(255);
	
	string info = "FPS:        "     + ofToString(ofGetFrameRate(), 0) + "\n";
	info += "Timer:      "           + ofToString(ofGetElapsedTimeMillis() / 1e+3, 2) + " seconds\n";
	info += "----------------------------System data------------------------"     "\n";
	info += "Molar calorific capacity at constant volume Cv: " + ofToString(CV, 2) + " J/mol*K\n";
	info += "Total Kinetic Energy: " + ofToString(K/1e+2, 2) + " J\n";
	info += "Internal Energy: "      + ofToString(U / 1e+2, 2) + " J\n";
	info += "Heat: "                 + ofToString(Q / 1e+2, 2) + " J\n";
	info += "Work: "                 + ofToString(W, 2) + " J\n";
	info += "Total Temperature: "    + ofToString(T / 1e+2, 2) + " K\n";

	ofDrawBitmapString(info, 20, 20);
}
//------------------------------------------------------------------------------------------------------
void ofApp::keyPressed(int key){
	switch (key) {
	case 100: // d
	debug = !debug;
	}
}
//------------------------------------------------------------------------------------------------------
