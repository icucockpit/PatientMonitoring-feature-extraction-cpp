/*
 * interestpoint.h
 *
 *  Created on: Aug 3, 2017
 *      Author: aed
 */

#ifndef INTERESTPOINT_H_
#define INTERESTPOINT_H_

//Includes ONLY needed for interestpoint.h
//OpenCV
#include "opencv2/core.hpp"
//C++
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <iomanip>
//Boost
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

//#include "boost/date_time/posix_time/posix_time_types.hpp" //no i/o just types
#include "boost/date_time/posix_time/posix_time.hpp" //include all types plus i/o

namespace ip {

struct interestpoint {
  // Interest point type
  // TIP = 0
  // BRISK = 1
  int iptype;

  // Position
  int x;
  int y;
  int frameno;
  boost::posix_time::ptime video_timestamp;

  // TIP related
  float sigma_l_sq; // corresponds to spatial size for scale-space representation
  float tau_l_sq;
  float s; // TIP scaling factor sigma_i_sq = s * sigma_l_sq; tau_i_sq = s * tau_l_sq
  float H_value; // H value at local maximum

  //BRISK related
  //based on float cv::KeyPoint::size
  //diameter of the meaningful keypoint neighborhood ??
  float brisk_size;

  // HOG descriptors
  int hog_patch_size;
  std::vector<float> HOG_descriptors;
  std::vector<cv::Point> HOG_desc_locations; // e.g. HOG descriptor locations

  // BRISK descriptors
  std::vector<float> BRISK_descriptors;

  // Other
  std::string roilabel;
};


//void drawIP(cv::Mat& mt_dst, int x, int y, float size_r);
//void drawIP(cv::Mat& mt_dst, int x, int y, float size_r, cv::Scalar color);
//void drawIP(cv::Mat& mt_dst, int x, int y, float size_r, cv::Scalar color, const float& scale);
//void drawAllTIPs(cv::Mat& mt_dst, std::vector<interestpoint>& vc_ip, const float& scale);
//void drawAllBRISKs(cv::Mat& mt_dst, std::vector<interestpoint>& vc_ip, const float& scale);

void updateTimestampToAllIPs(std::vector<interestpoint>& vc_ip, boost::posix_time::ptime timestamp);
void appendIPsAndDescrToCsvFile(std::ofstream& ofs_myfile, const std::vector<interestpoint>& vc_ip, std::string& attr_a, std::string& attr_b, std::string& attr_c, std::string& attr_d, std::string& attr_e, std::string& attr_f);

void appendIPsAndHOGDescrToCsvFile(std::ofstream& ofs_myfile, const std::vector<interestpoint>& vc_ip, const std::vector<std::string>& vc_labels);
void appendIPsAndBRISKDescrToCsvFile(std::ofstream& ofs_myfile, const std::vector<interestpoint>& vc_ip, const std::vector<std::string>& vc_labels);
void appendIPsAndHOGBRISKDescrToCsvFile(std::ofstream& ofs_myfile, const std::vector<interestpoint>& vc_ip, const std::vector<std::string>& vc_labels);

void convertIPVectorToBRISKKeyPoints(const std::vector<interestpoint>& vc_ip, std::vector<cv::KeyPoint>& vc_kp_in);
void convertBRISKKeyPointsAndBRISKDescriptorsToIPVector(const std::vector<cv::KeyPoint>& vc_kp_in, const cv::Mat& decr_in, const int& brisk_descr_size, const int& frameno, std::vector<interestpoint>& vc_ip);
void convertBRISKKeyPointsToIPVector(const std::vector<cv::KeyPoint>& vc_kp_in, const int& frameno, std::vector<interestpoint>& vc_ip_out);
void addBRISKDescrToIPVector(const std::vector<cv::KeyPoint>& vc_kp_in, const cv::Mat& decr_in, const int& brisk_descr_size, const int& frameno, std::vector<interestpoint>& vc_ip);

void checkIfIPsAreInsideRectAndGiveLabel(std::vector<interestpoint>& vc_ip, const cv::Rect& roirect, const std::string& ifinsidelabel, const std::string& ifoutsidelabel);
void checkIfIPsBRISKEmptyAndAddDummyIP(std::vector<interestpoint>& vc_ip, const int& frameno, const int& descr_size);
void checkIfIPsHOGEmptyAndAddDummyIP(std::vector<interestpoint>& vc_ip, const int& frameno, const int& descr_size);

} // namespace ip

//non-intrusive serialization for struct interestpoint and cv::Point
namespace boost {
	namespace serialization {
		template<class Archive>
		inline void serialize(Archive & ar, ip::interestpoint & t, const unsigned int file_version){
			ar & t.x;
			ar & t.y;
			ar & t.frameno;
			ar & t.sigma_l_sq;
			ar & t.tau_l_sq;
			ar & t.s;
			ar & t.H_value;
			ar & t.hog_patch_size;
			ar & t.HOG_descriptors;
			ar & t.HOG_desc_locations; //TODO UPDATE WITH ALL struct elements
		}
		template<class Archive>
		inline void serialize(Archive & ar, cv::Point & t, const unsigned int file_version){
			ar & t.x;
			ar & t.y;
		}
	} // namespace serialization
} // namespace boost


#endif /* INTERESTPOINT_H_ */
