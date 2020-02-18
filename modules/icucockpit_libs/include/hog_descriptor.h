/*
 * hogdescriptor.h
 *
 *  Created on: Jul 3, 2017
 *      Author: aed
 */

#ifndef HOGDESCRIPTOR_H_
#define HOGDESCRIPTOR_H_

//Includes ONLY needed for hogdescriptor.h
// C++
#include <vector>
// OpenCV
#include "opencv2/core.hpp"
// Mat
#include "interest_point.h" //for interest point structure stip::interestpoint


namespace hog {

void testHOG();
int testHOGOnVideo(cv::CommandLineParser& parser);
int testHOGOnVideo2(cv::CommandLineParser& parser);

void computeDense3x3HOG(const cv::Mat& mt_src, const int& framenumber, std::vector<ip::interestpoint>& ls_locations, const int& step, const int& patchsize, const int& numbins);
void compute3x3HOGForAllTIPs(const cv::Mat& mt_src, std::vector<ip::interestpoint>& vc_ip_io, const int& numbins);
void compute3x3HOGForAllBRISKIPs(const cv::Mat& mt_src, std::vector<ip::interestpoint>& vc_ip_io, const int& numbins);
int compute3x3HOGAtPoint(const cv::Mat& mt_src, std::vector<float>& vc_descriptors, std::vector<cv::Point>& vc_locations, const cv::Point& position, const int& patchsize, const int& numbins);
void compute3x3HOG( const cv::Mat& mt_src, std::vector<float>& vc_descriptors, std::vector<cv::Point>& vc_locations, const int& numbins);

void draw3x3HOGVisualizationForAllIPs(cv::Mat& mt_ioimage, std::vector<ip::interestpoint>& ls_STIPs, const int& numbins, const float& scale);
void draw3x3HOGVisualizationAtPoint(cv::Mat& mt_ioimage, const std::vector<float>& vc_descriptors, const std::vector<cv::Point>& vc_locations, const cv::Point& position, const int& patchsize, const int& numbins, const float& scale);
cv::Mat get3x3HOGVisualizationAtPoint(const cv::Mat& mt_src, const std::vector<float>& vc_descriptors, const std::vector<cv::Point>& vc_locations, const cv::Point& position, const int& patchsize, const int& numbins);
cv::Mat get3x3HOGVisualization(const cv::Mat& color_origImg, const std::vector<float>& descriptorValues, const int& numbins);

class hogdescriptor {
public:
	hogdescriptor();
	virtual ~hogdescriptor();
};

} /* namespace hog */

#endif /* HOGDESCRIPTOR_H_ */
