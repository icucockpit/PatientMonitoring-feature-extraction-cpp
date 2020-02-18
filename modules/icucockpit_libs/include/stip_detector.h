/*
 * stipdetector.h
 *
 *  Created on: May 19, 2017
 *      Author: aed
 */

#ifndef STIPDETECTOR_H_
#define STIPDETECTOR_H_

//#define SHOWVIDEO
//#define TIME_MEASUREMENTS

// Includes ONLY needed for stipdetector.h
// C++
#include <vector>
#include <list>
#include <fstream>
// OpenCV
#include "opencv2/core.hpp"
// Boost
#include "boost/circular_buffer.hpp"
// Mat interest point
#include "interest_point.h"

namespace tip {

struct params {
	int frame_height, frame_width;
	double startingpostition;
	float sigma_sq, tau_sq, kappa, thresh;
};

int testSTIPDetection_V1(cv::CommandLineParser& parser);

class stipdetector {
public:
	//stipdetector(const int& frame_height, const int& frame_width, const int& startingpostition, const float& sigma_sq, const float& tau_sq, const float& kappa, const float& thresh);
	stipdetector(const tip::params& stip_parameters);
	virtual ~stipdetector();
	std::vector<ip::interestpoint> processFrame(const cv::Mat& mt_source, std::vector<ip::interestpoint>& ls_IPs);
	int getDelay(void){return delay; };

	//DEV:
	int processFrameTemporal(const cv::Mat& mt_source, std::vector<ip::interestpoint>& ls_ip, double& corresponfing_frameno);

private:
	//UNDER CONSTRUCTION
	std::vector<ip::interestpoint> processFrameFast(const cv::Mat& mt_source, std::vector<ip::interestpoint>& ls_IPs);

	cv::Mat returnGaussianKernel(int ksize, float sigma_sq);
	void temporalConvolution(boost::circular_buffer<cv::Mat>& cb_src, cv::Mat& mt_dst, const cv::Mat& mt_kernel);
	void temporalBoxFilter(boost::circular_buffer<cv::Mat>& cb_src, cv::Mat& mt_dst, const int& boxkernelsize, const int& n_iter);
	void rowBoxFilter(cv::Mat& mt_src, cv::Mat& mt_dst, const int& boxkernelsize, const int& n_iter);
	void colBoxFilter(cv::Mat& mt_src, cv::Mat& mt_dst, const int& boxkernelsize, const int& n_iter);
	void findLocalMaxima3D(boost::circular_buffer<cv::Mat>& cb_src, std::vector<ip::interestpoint>& ls_local_maxima, int xy_neighbourhood, int t_neighbourhood, int framenumber);
	void findLocalMaxima3DThreshold(boost::circular_buffer<cv::Mat>& cb_src, std::vector<ip::interestpoint>& ls_local_maxima, int& xy_neighbourhood, int& t_neighbourhood, double& frameno_being_processed, float& thresh);

	unsigned int SPATIAL_GAUSSIAN_KERNEL_SIZE, SPATIAL_RESIZED_GAUSSIAN_KERNEL_SIZE, TEMPORAL_GAUSSIAN_KERNEL_SIZE, TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE, DERIV_KERNEL_SIZE, HARRIS_MAT_BUFFER_SIZE, center_frame_index;

	double frameno_being_processed, frameno_starting_position;
	int delay, xy_localmaxneighbourhood, t_localmaxneighbourhood;

	//Necessary circular buffers holding cv::Mat.
	boost::circular_buffer<cv::Mat> cb_input_frame, cb_input_frame_flt, cb_LxLx, cb_LyLy, cb_LtLt, cb_LxLy, cb_LxLt, cb_LyLt, cb_H;

	float resize_factor, sigma_l, tau_l;
	int height_resized, width_resized;

	//Gaussian kernel for temporal convolution, normalized to sum(all elements) = 1
	float tau_l_sq;
	cv::Mat mt_1D_Gauss_kernel_tau_l_sq;

	//Gaussian kernel for spatial convolution, normalized to sum(all elements) = 1
	float sigma_l_sq;
	cv::Mat mt_1D_Gauss_kernel_sigma_l_sq;

	//Scaled Gaussian kernel for temporal convolution, normalized to sum(all elements) = 1
	float s, tau_i_sq;
	cv::Mat mt_1D_Gauss_kernel_tau_i_sq;

	//Scaled Gaussian kernel for spatial convolution, normalized to sum(all elements) = 1
	float sigma_i_sq;
	cv::Mat mt_1D_Gauss_kernel_sigma_i_sq;

	//Partial normalized derivative kernel (Sobel...?)
	cv::Mat mt_sob_deriv_kernel_x, mt_sob_deriv_kernel_y, mt_deriv_kernel_x, mt_deriv_kernel_y;

	cv::Mat mt_Lt, mt_Lx, mt_Ly;

	cv::Mat mt_LxLx_filtered_temporally, mt_LxLx_filtered_temporally_and_spacially, mt_LyLy_filtered_temporally, mt_LyLy_filtered_temporally_and_spacially, mt_LtLt_filtered_temporally, mt_LtLt_filtered_temporally_and_spacially, mt_LxLy_filtered_temporally,t mt_LxLy_filtered_temporally_and_spacially, mt_LxLt_filtered_temporally, mt_LxLt_filtered_temporally_and_spacially, mt_LyLt_filtered_temporally, mt_LyLt_filtered_temporally_and_spacially;

	cv::Mat mt_tmp1, mt_tmp2;
	float k; // Harris kappa
	float threshold; //threshold for omitting weak points

	//std::vector<ip::interestpoint> ls_IPs;

	//For ProcessFrameFast in addition
	cv::Mat mt_current_frame_flt_temp;
	int filtersizewatch;
};

} /* namespace tip */

#endif /* STIPDETECTOR_H_ */
