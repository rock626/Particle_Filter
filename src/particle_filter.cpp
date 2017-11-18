/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#define _USE_MATH_DEFINES

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	default_random_engine gen;
	double std_x, std_y, std_theta;

	// This line creates a normal (Gaussian) distribution for x with mean as given value
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	num_particles=100;
	
	for(int i=0;i<num_particles;i++){
		Particle  p;
		p.id =i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1;

		particles.push_back(p);
	}
	// cout << "init done" << endl;
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	default_random_engine gen;

	// This line creates a normal (Gaussian) distribution for x, with mean as 0
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);

	// cout << "Prediction" << endl;
	for(int i=0;i<num_particles;i++){
		
		//add measurements to existing values
		double val1 = particles[i].theta + (yaw_rate * delta_t);
		if(yaw_rate == 0) {
			particles[i].x += (velocity * delta_t) * cos(particles[i].theta);
			particles[i].y += (velocity * delta_t) * sin(particles[i].theta);
			//particles[i].theta = particles[i].theta;
		} else {
			particles[i].x += (velocity/yaw_rate) * (sin(val1) -  sin(particles[i].theta));
			particles[i].y += (velocity/yaw_rate) * (cos(particles[i].theta) - cos(val1));
			particles[i].theta  += yaw_rate * delta_t;
		}

		//add noise
		particles[i].x +=  dist_x(gen);
		particles[i].y += dist_y(gen);
		particles[i].theta += dist_theta(gen);

		// cout << "Particle " << i << ": x:" <<particles[i].x << " y:" << particles[i].y << " theta:" << particles[i].theta << " weight:" << particles[i].weight << endl;
	}

}



void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

double euclideanDistance(double x1, double y1, double x2, double y2)
{
    double x = x1 - x2;
    double y = y1 - y2;
    double dist;

    dist = pow(x,2)+pow(y,2);           //calculating distance by euclidean formula
    dist = sqrt(dist);                  //sqrt is function in math.h

    return dist;
}

Map::single_landmark_s nearestLandMark(double sensor_range,double x_map,double y_map, const Map &map_landmarks) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	double min = std::numeric_limits<double>::max();
	// cout << "Initial min value" << min << endl;
	Map::single_landmark_s nearestLandMark;
	
	for(int i=0;i<map_landmarks.landmark_list.size();i++){
		Map::single_landmark_s landmark = map_landmarks.landmark_list[i];
		if(abs(landmark.x_f - x_map) <= sensor_range && abs(landmark.y_f - y_map) <= sensor_range){

			//eculidean distance
			double distance = euclideanDistance(x_map,y_map,landmark.x_f,landmark.y_f);
			if(distance < min){
				min = distance;
				nearestLandMark = landmark;
			}	
		}else{
			// cout << "out of range x_map:" << x_map << " y_map" << y_map << " landmark x:"<< landmark.x_f << " y:" << landmark.y_f <<endl;
		}

	}
	 return nearestLandMark;

}



void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	//clear weights vector;
	weights.clear();

	//for each particle
	double x_map;
	double y_map;
	// cout << "Update weights" << endl;
	double gaussian_norm = 1/(2*M_PI * std_landmark[0] * std_landmark[1]);
	for(int i=0;i<num_particles;i++){
		//for each observation
		// cout << "Particle " << i << ": x:" <<particles[i].x << " y:" << particles[i].y << " theta:" << particles[i].theta << endl;

		double particle_weight=1.0;
		std::vector<int> associations_new;
		std::vector<double> sense_x_new;
		std::vector<double> sense_y_new;

		std::vector<double> p_weights;
		for(int j=0;j<observations.size();j++){
			LandmarkObs l_obs = observations[j];
			// cout << "Observation id:"<< j << " x:" << l_obs.x << " y:" << l_obs.y << endl;
			x_map = particles[i].x + (cos(particles[i].theta) * l_obs.x)  - (sin(particles[i].theta)*l_obs.y);
			y_map = particles[i].y + (sin(particles[i].theta) * l_obs.x)  + (cos(particles[i].theta)*l_obs.y);
			
			// cout << "Mapping:"<< j << " x_map:" << x_map << " y_map:" << y_map << endl;
			Map::single_landmark_s mu = nearestLandMark(sensor_range,x_map,y_map,map_landmarks);
			// cout << "NearestLandMark " << mu.id_i << ": x:" <<mu.x_f << " y:" << mu.y_f << endl;

			//set associations and sense_x, sense_y
			associations_new.push_back(mu.id_i);
			sense_x_new.push_back(x_map);
			sense_y_new.push_back(y_map);

			//calculate weight
			double exp_x = pow((x_map - mu.x_f),2)/(2 * pow(std_landmark[0],2));
			double exp_y = pow((y_map - mu.y_f),2)/(2 * pow(std_landmark[1],2));
			double exponent = exp_x+exp_y;
			double weight = gaussian_norm * exp(-exponent); 
			// cout << "gaussian_norm:"<< gaussian_norm <<" exp_x:" << exp_x << " exp_y:" << exp_y << " exponent:" << exponent << endl;
			// cout << "weight:" << weight << endl;
			//final weight is by multiplying all the calculated measurement probabilities together.
			p_weights.push_back(weight);

		}

		double c_weight =1.0;
		for(int k=0;k<p_weights.size();k++){
			c_weight *= p_weights[k];
		}
		// cout << "cumulative weight:" << c_weight << endl;

		SetAssociations(particles[i],associations_new,sense_x_new,sense_y_new);
		particles[i].weight = c_weight;

		weights.push_back(c_weight);

		cout << endl;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	std::random_device rd;
    std::mt19937 gen(rd());
	std::discrete_distribution<double> distribution(weights.begin(),weights.end());

	// cout << "Before resampling: " << endl;
	// for(int i=0;i<num_particles;i++){
	// 	// cout << weights[i] << "  " << endl;
	// }

	std::vector<Particle> p_new;

	// cout << "After resampling: " << endl;

	for(int i=0;i<num_particles;i++){
		int index = distribution(gen);
		p_new.push_back(particles[index]);
		// cout << particles[index].weight << "  " << endl;
	}

	particles = p_new;
	cout << endl;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
