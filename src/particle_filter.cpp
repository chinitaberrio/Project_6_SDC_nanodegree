/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

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
  num_particles = 200;  // Set the number of particles

  // set standard deviations for x, y and therta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  std::default_random_engine gen;

  // initializing all the particles
  for (int i = 0; i < num_particles; ++i)
  {
    particles.push_back(Particle());
    particles[i].id = i;
    particles[i].x  = dist_x(gen);
    particles[i].y  = dist_y(gen);
    particles[i].theta  = dist_theta(gen);
    particles[i].weight = 1.0;
    weights.push_back(1.0);

   }

 is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  normal_distribution<double> dist_x(0.0, std_pos[0]);
  normal_distribution<double> dist_y(0.0, std_pos[1]);
  normal_distribution<double> dist_theta(0.0, std_pos[2]);

  std::default_random_engine gen;

  for (int i = 0; i < particles.size(); ++i)
  {
    // checking for yawrate (I can't divide by 0.0)
    if (fabs(yaw_rate) > 0.001){
      particles[i].x  = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta)) + dist_x(gen);
      particles[i].y  = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta)-cos(particles[i].theta+yaw_rate*delta_t)) + dist_y(gen);
      particles[i].theta  = particles[i].theta + yaw_rate*delta_t + dist_theta(gen);
    }
    else { // keeping the same direction
      particles[i].x  = particles[i].x + (velocity*delta_t*cos(particles[i].theta)) + dist_x(gen);
      particles[i].y  = particles[i].y + (velocity*delta_t*sin(particles[i].theta)) + dist_y(gen);
      particles[i].theta  = particles[i].theta + dist_theta(gen);
    }
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
  // iterating throught the observations
  for(int i = 0; i < observations.size(); i++)
  {
    // initialize the min distance with a very large number
    double min_dist = 987563241.0;

    // iterating through the predicted landmarks (map)
    for (int j = 0; j < predicted.size(); ++j) {
      // calculating the distance between the two
      double d = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
      // checking for the min distance to the observation
      if (d <= min_dist){
        min_dist = d;
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

  // for each particle
  for (int i = 0; i < particles.size(); ++i) {

    // Transforming the observations to the map frame
    vector<LandmarkObs> observations_in_map;
    for (int j = 0; j < observations.size(); ++j) {
      observations_in_map.push_back(LandmarkObs());
      observations_in_map[j].x = particles[i].x + (cos(particles[i].theta) * observations[j].x) - (sin(particles[i].theta) * observations[j].y);
      observations_in_map[j].y = particles[i].y + (sin(particles[i].theta) * observations[j].x) + (cos(particles[i].theta) * observations[j].y);
    }

    // Obtaining the landmarks within the sensor range
    vector<LandmarkObs> maplm_within_range;
    for (int j = 0; j < map_landmarks.landmark_list.size(); ++j) {
      if (dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f) <= sensor_range) {
        LandmarkObs map_obs = {};
        map_obs.x  = double(map_landmarks.landmark_list[j].x_f);
        map_obs.y  = double(map_landmarks.landmark_list[j].y_f);
        map_obs.id = map_landmarks.landmark_list[j].id_i;
        maplm_within_range.push_back(map_obs);
      }
    }

    // Data association
    dataAssociation(maplm_within_range, observations_in_map);


    particles[i].weight = 1.0;

    // Probabilities
    for (int j = 0; j < observations_in_map.size(); j++){
      // Current observation
      LandmarkObs current_obs = observations_in_map[j];
      // Nearest map landmark
      // NOTE: Map id starts with 1
      LandmarkObs nearest_map = {};
      nearest_map.x = map_landmarks.landmark_list[(current_obs.id-1)].x_f;
      nearest_map.y = map_landmarks.landmark_list[(current_obs.id-1)].y_f;

      // Calculating multivariate gaussian
      // Normalization term
      double gauss_norm;
      gauss_norm = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
      // Exponent term
      double exponent;
      exponent = (pow(current_obs.x - nearest_map.x, 2) / (2 * pow(std_landmark[0], 2)))
                  + (pow(current_obs.y - nearest_map.y, 2) / (2 * pow(std_landmark[1], 2)));
      double weight = gauss_norm * exp(-exponent);

      // Combining probabilities
      particles[i].weight *= weight;

    }
    weights[i]=particles[i].weight;

  }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  // Using randomg generator to resample the particles
  std::default_random_engine gen;
  std::discrete_distribution<int> index_dist(weights.begin(), weights.end());

  // vector of new particles
  vector<Particle> new_particles(particles.size());

  // populating the new particles vector
  for (int i = 0; i < particles.size(); ++i) {
    new_particles[i] = particles[index_dist(gen)];
  }

  //replacing the particles for the new resampled particles
  particles = new_particles;

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
