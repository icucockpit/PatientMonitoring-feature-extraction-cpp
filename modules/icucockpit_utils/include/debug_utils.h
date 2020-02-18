#ifndef DEBUG_UTILS_H_
#define DEBUG_UTILS_H_

//Includes only needed for debug_utils.h
//C++
#include <string>

//OpenCV
#include "opencv2/core.hpp"

//helpful utilities
void dispMatProp(const cv::Mat& mt_src, const std::string& name);
void dispMatDim(const cv::Mat& mt_src, const std::string& name);
void dispMat(const cv::Mat& mt_src, const std::string& name);
void dispMinMax(const cv::Mat& mt_src, const std::string& name);
void showFloatMat(const cv::Mat& mt_src, const std::string& name);
void drawCross(cv::Mat& mt_dst, int x, int y, int length = 4);


#endif /* DEBUG_UTILS_H_ */
