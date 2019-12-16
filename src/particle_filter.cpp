/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"
#include <stdio.h> 
#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <limits>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  is_initialized = true;
  std::default_random_engine gen;
  
  // This line creates a normal (Gaussian) distribution for x,y,theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

 for(int j=0; j<num_particles; j++){
     // Sample from these normal distributions like this: 
    Particle p;
    p.id = j;
    p.x = x;
    p.y = y;
    p.theta = theta;
    p.weight = 1.0;
   
   	 p.x = dist_x(gen);
     p.y = dist_y(gen);
     p.theta = dist_theta(gen);
   
    particles.push_back(p);
   }
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
      std::random_device rd;
    std::default_random_engine gen;

    for (auto &p: particles){

std::normal_distribution<double> dist_x(0, std_pos[0]); //Gaussian noise initializing based on estimates by sampling from Gaussian dist
    std::normal_distribution<double> dist_y(0, std_pos[1]); //x,y,theta are values from GPS
    std::normal_distribution<double> dist_theta(0, std_pos[2]);

        if(fabs(yaw_rate) > 0.0001){
            p.x += velocity/yaw_rate* (sin(p.theta + yaw_rate * delta_t) - sin(p.theta)) + dist_x(gen);
            p.y += velocity/yaw_rate * (cos(p.theta) - cos(p.theta + yaw_rate* delta_t)) + dist_y(gen);
        }
        else{
            p.x += velocity * delta_t *cos(p.theta) + dist_x(gen);
            p.y += velocity * delta_t * sin(p.theta)+ dist_y(gen);
        }
        p.theta = p.theta + yaw_rate*delta_t + dist_theta(gen);
    }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
      
   for(unsigned int i=0; i < observations.size(); i++){
     
     double distance = std::numeric_limits<double>::max(); 
     //observations[i].id = -1; // If no associated landmark is found assign value = -1.
     
     for (unsigned int j=0; j< predicted.size(); j++){
       double cur_dist = dist(predicted[j].x,predicted[j].y,observations[i].x,observations[i].y); 
       if ( cur_dist < distance ){
         distance = cur_dist;
         observations[i].id = predicted[j].id;
         }
      }     
   }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  // Step One: Transform obsevation into particle coordinate system 
  // Step Two: Associate new map landmarks
  // step Three: compute weights 

 
  double distance;
  double weight_normalizer = 0;
  double sig_x = std_landmark[0];
  double sig_y = std_landmark[1];
  double delta_y, delta_x;
  double gauss_normalizer = (1/(2*M_PI*sig_x*sig_y));
  
  

  for(int i=0; i < num_particles; i++){
    /**
     Find the list of landmarks that the sensor can detect for the current particle. To directly use data assiciation fucntion above.
     **/
    
    // get the particle x, y coordinates
    double p_x = particles[i].x;
    double p_y = particles[i].y;
    double p_theta = particles[i].theta;
    
    vector<LandmarkObs> valid_landmarks;

    for (unsigned int k = 0; k < map_landmarks.landmark_list.size(); k++){
      
        distance = dist(p_x,  p_y,  (double) map_landmarks.landmark_list[k].x_f, (double) map_landmarks.landmark_list[k].y_f);
        if (distance <= sensor_range){
           valid_landmarks.push_back(LandmarkObs{map_landmarks.landmark_list[k].id_i, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f });
        }
    }
    
    // for a single partcile, transform obersvations to obtain measured landmark coordinates (w.r.t particle).
    vector<LandmarkObs> new_map_coordinates; // coordinates of possible landmark w.r.t a specific particle.
    for (unsigned int j = 0; j < observations.size(); j++){
      //x_map = x_part + (cos(theta) * x_obs) - (sin(theta) * y_obs);
        double x_obv = p_x + (cos(p_theta)*observations[j].x - sin(p_theta)*observations[j].y);
      //y_map = y_part + (sin(theta) * x_obs) + (cos(theta) * y_obs);
        double y_obv = p_y + (sin(p_theta)*observations[j].x + cos(p_theta)*observations[j].y);
       new_map_coordinates.push_back(LandmarkObs{observations[j].id,x_obv, y_obv});
    }
    
    dataAssociation(valid_landmarks, new_map_coordinates);
     
    /**
    Computaiton of weights based on the mapped landmarks for a single particle
   **/
    particles[i].weight = 1.0 ;
    bool nearest_ngbr_found = false;
    for (unsigned int j = 0; j < new_map_coordinates.size(); j++){     
      for(unsigned int k  = 0; k < valid_landmarks.size(); k++){
        if (valid_landmarks[k].id == new_map_coordinates[j].id){
          delta_x = (valid_landmarks[k].x - new_map_coordinates[j].x);
          delta_y = valid_landmarks[k].y - new_map_coordinates[j].y;
          //std::cout << "weight "<< particles[i].weight << std::endl << std::flush; 
          //std::cout << "delta_x "<< delta_x << std::endl << std::flush;
          //std::cout << "sig_x "<< exp(-(pow((delta_x/sig_x),2))) << std::endl << std::flush;
          //std::cout << "gauss_normalizer "<< gauss_normalizer << std::endl << std::flush;
          particles[i].weight *= (gauss_normalizer * exp(-( pow((delta_x/sig_x),2) + pow((delta_y/sig_y),2) ) ) );
          // particles[i].weight *= 1 ;
          nearest_ngbr_found = true;
          //std::cout << "new landmark found" << std::endl << std::flush;
          //std::cout << "new landmark "<< new_map_coordinates[j].id << std::endl << std::flush;
          //std::cout << "new landmark "<< valid_landmarks[k].id<< std::endl << std::flush;
          //std::cout << "modified_weight "<< particles[i].weight << std::endl << std::flush;
          break;}
        }
      
       if (!nearest_ngbr_found){
         //std::cout << "no ngbr found "<< std::endl << std::flush;
            particles[i].weight = 0;
          }
      
      }
 
    //std::cout << "particle weight " << particles[i].weight << std::endl << std::flush;
    //std::cout << "number of valid landmark " << valid_landmarks.size() << std::endl << std::flush;
    //std::cout << "number of observations  " << observations.size() << std::endl << std::flush;
    
     weight_normalizer += particles[i].weight;
   }
  
  for (int i = 0; i < num_particles; i++) {
    particles[i].weight /= weight_normalizer;
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
   // Get weights and max weight.
  vector<double> weights;
  double Weight_max = std::numeric_limits<double>::min();
  for(int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
    if ( particles[i].weight > Weight_max ) {
      Weight_max = particles[i].weight;
    }
  }

  
  // Obtain distributions and pick index value.
  std::default_random_engine gen;
  std::uniform_real_distribution<double> distDouble(0.0, 2*Weight_max);
  std::uniform_int_distribution<int> distInt(0, num_particles - 1); 
  
  int index = distInt(gen);

  // From class implement the algorithm for resampling.
  double beta = 0.0;
  vector<Particle> Particles_resampled;
  for(int i = 0; i < num_particles; i++) {
    beta += distDouble(gen);
    while( beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    //std::cout << "particle index " <<index << std::endl << std::flush;
    Particles_resampled.push_back(particles[index]);
  }
  particles = Particles_resampled;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}