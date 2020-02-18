/*
 * stipdetector.cpp
 *
 *  Created on: May 19, 2017
 *      Author: aed
 */

#include "stip_detector.h"

// ADDITIONAL includes needed by stip_detector.cpp
// C++
#include <math.h>
#include <iostream>
#include <time.h> //for calculating elapsed time
// OpenCV
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"
// Mat debug utilities
#include "debug_utils.h"

namespace tip {

int testSTIPDetection_V1(cv::CommandLineParser& parser){
	cv::VideoCapture cap;

	std::string in_file_name = parser.get<std::string>("input_video_file");

	cap.open(in_file_name);

	if( !cap.isOpened() )
	{
		std::cout << "***Could not initialize capturing...***\n"<<"ERROR: Unable to open video file."<<std::endl;
		return -1;
	}

	std::string out_file_name = parser.get<std::string>("output_file");

	std::ofstream myfile;
	myfile.open ((char*)out_file_name.c_str(), std::ios::trunc); //Any current content is discarded
	if (!myfile.is_open())
	{
		std::cerr<<"ERROR: Unable to open output file for STIP results."<<std::endl;
		return -1;

	}
	myfile<<"# "<<in_file_name<<std::endl;

	//IMPORTANT resize factor of input video to make things faster
	float resize = 0.5; //0.5 = the half
	//IMPORTANT

	int capheight = cv::saturate_cast<int>(cap.get(CV_CAP_PROP_FRAME_HEIGHT)*resize);
	int capwidth = cv::saturate_cast<int>(cap.get(CV_CAP_PROP_FRAME_WIDTH)*resize);

	tip::params stip_params;
	stip_params.kappa = 0.0005;
	stip_params.thresh = 1e-14;
	stip_params.frame_height = capheight;
	stip_params.frame_width = capwidth;
	stip_params.startingpostition = 0;

	//Instantiate all STIP detectors
	stip_params.sigma_sq = 2;
	stip_params.tau_sq = 2;
	tip::stipdetector mystipdetector2_2(stip_params);
	int delay = mystipdetector2_2.getDelay();
	stip_params.sigma_sq = 4;
	stip_params.tau_sq = 2;
	tip::stipdetector mystipdetector4_2(stip_params);
	stip_params.sigma_sq = 8;
	stip_params.tau_sq = 2;
	tip::stipdetector mystipdetector8_2(stip_params);
	stip_params.sigma_sq = 16;
	stip_params.tau_sq = 2;
	tip::stipdetector mystipdetector16_2(stip_params);
	stip_params.sigma_sq = 32;
	stip_params.tau_sq = 2;
	tip::stipdetector mystipdetector32_2(stip_params);
	stip_params.sigma_sq = 64;
	stip_params.tau_sq = 2;
	tip::stipdetector mystipdetector64_2(stip_params);
	stip_params.sigma_sq = 128;
	stip_params.tau_sq = 2;
	tip::stipdetector mystipdetector128_2(stip_params);
	stip_params.sigma_sq = 256;
	stip_params.tau_sq = 2;
	tip::stipdetector mystipdetector256_2(stip_params);
	stip_params.sigma_sq = 512;
	stip_params.tau_sq = 2;
	tip::stipdetector mystipdetector512_2(stip_params);


	//--------------------------------------------------
	//thresh = 1e-8;
	stip_params.sigma_sq = 2;
	stip_params.tau_sq = 4;
	tip::stipdetector mystipdetector2_4(stip_params);
	stip_params.sigma_sq = 4;
	stip_params.tau_sq = 4;
	tip::stipdetector mystipdetector4_4(stip_params);
	stip_params.sigma_sq = 8;
	stip_params.tau_sq = 4;
	tip::stipdetector mystipdetector8_4(stip_params);
	stip_params.sigma_sq = 16;
	stip_params.tau_sq = 4;
	tip::stipdetector mystipdetector16_4(stip_params);
	stip_params.sigma_sq = 32;
	stip_params.tau_sq = 4;
	tip::stipdetector mystipdetector32_4(stip_params);
	stip_params.sigma_sq = 64;
	stip_params.tau_sq = 4;
	tip::stipdetector mystipdetector64_4(stip_params);
	stip_params.sigma_sq = 128;
	stip_params.tau_sq = 4;
	tip::stipdetector mystipdetector128_4(stip_params);
	stip_params.sigma_sq = 256;
	stip_params.tau_sq = 4;
	tip::stipdetector mystipdetector256_4(stip_params);
	stip_params.sigma_sq = 512;
	stip_params.tau_sq = 4;
	tip::stipdetector mystipdetector512_4(stip_params);

/*#ifdef VISUALIZE_MAIN
	cv::VideoCapture cap2;
	cap2.open(in_file_name);
	if( !cap2.isOpened() )
	{
		std::cout << "***Could not initialize capturing...***\n";
		std::cerr<<"ERROR: Unable to open video file."<<std::endl;
		return -1;
	}
	cv::Mat captured_frame_vis0;// = cv::Mat::zeros(cap2.get(CV_CAP_PROP_FRAME_HEIGHT), cap2.get(CV_CAP_PROP_FRAME_WIDTH), cap2.get(CV_CAP_PROP_FORMAT));
	cv::Mat captured_frame_vis = cv::Mat::zeros(capheight, capwidth, cap2.get(CV_CAP_PROP_FORMAT));
	cv::namedWindow("Input_STIPS", cv::WINDOW_NORMAL);
	cv::resizeWindow("Input_STIPS", 1200, 400);
	bool paused = false;
	bool firstrun = true;
	int frameno_being_processed = 0;
#endif
*/
    
#ifdef TIME_MEASUREMENTS_MAIN
	//for calculating elapsed time
	clock_t ETSstart, ETSend;
	unsigned int fps_counter = 0;
	double duration = 0.0;
	double duration_sum = 0.0;
#endif


	std::vector<ip::interestpoint> ls_local_maxima2_2, ls_local_maxima4_2,ls_local_maxima8_2, ls_local_maxima16_2, ls_local_maxima32_2, ls_local_maxima64_2, ls_local_maxima128_2, ls_local_maxima256_2, ls_local_maxima512_2, ls_local_maxima2_4, ls_local_maxima4_4, ls_local_maxima8_4, ls_local_maxima16_4, ls_local_maxima32_4, ls_local_maxima64_4, ls_local_maxima128_4, ls_local_maxima256_4, ls_local_maxima512_4;

	cv::Mat mt_captured_frame0, mt_captured_frame;
	//***********************************************
	//Start the loop
	//***********************************************
	for(;;)
	{

#ifdef VISUALIZE_MAIN
		if(!paused)
		{
#endif
#ifdef TIME_MEASUREMENTS_MAIN
			ETSstart = clock();
#endif
			// Capture next frame
			cap >> mt_captured_frame0;
			if(mt_captured_frame0.empty())
				break;
			cv::resize(mt_captured_frame0, mt_captured_frame, cv::Size(), resize, resize, CV_INTER_AREA);
			// Calculate STIP points
			mystipdetector2_2.processFrame(mt_captured_frame, ls_local_maxima2_2);
			mystipdetector4_2.processFrame(mt_captured_frame, ls_local_maxima4_2);
			mystipdetector8_2.processFrame(mt_captured_frame, ls_local_maxima8_2);
			mystipdetector16_2.processFrame(mt_captured_frame, ls_local_maxima16_2);
			mystipdetector32_2.processFrame(mt_captured_frame, ls_local_maxima32_2);
			mystipdetector64_2.processFrame(mt_captured_frame, ls_local_maxima64_2);
			mystipdetector128_2.processFrame(mt_captured_frame, ls_local_maxima128_2);
			mystipdetector256_2.processFrame(mt_captured_frame, ls_local_maxima256_2);
			mystipdetector512_2.processFrame(mt_captured_frame, ls_local_maxima512_2);

			mystipdetector2_4.processFrame(mt_captured_frame, ls_local_maxima2_4);
			mystipdetector4_4.processFrame(mt_captured_frame, ls_local_maxima4_4);
			mystipdetector8_4.processFrame(mt_captured_frame, ls_local_maxima8_4);
			mystipdetector16_4.processFrame(mt_captured_frame, ls_local_maxima16_4);
			mystipdetector32_4.processFrame(mt_captured_frame, ls_local_maxima32_4);
			mystipdetector64_4.processFrame(mt_captured_frame, ls_local_maxima64_4);
			mystipdetector128_4.processFrame(mt_captured_frame, ls_local_maxima128_4);
			mystipdetector256_4.processFrame(mt_captured_frame, ls_local_maxima256_4);
			mystipdetector512_4.processFrame(mt_captured_frame, ls_local_maxima512_4);
			//calculate descriptors

#ifdef TIME_MEASUREMENTS_MAIN
			ETSend = clock();
			duration = (double(ETSend - ETSstart) / CLOCKS_PER_SEC);
			duration_sum += duration;
			fps_counter++;
			if (fps_counter >= 10)
			{
				duration_sum = 0.0;
				fps_counter = 0.0;
			}
#endif


#ifdef VISUALIZE_MAIN
			if (frameno_being_processed - delay >= 0)// start visualization
			{
				if (frameno_being_processed - delay == 0)//this needs to be only set once, hence "=="
					cap2.set(CV_CAP_PROP_POS_FRAMES, frameno_being_processed - delay);

				cap2 >> captured_frame_vis0;
				if(captured_frame_vis0.empty())
					break;

				captured_frame_vis = captured_frame_vis0;

				/*drawAllSTIPs(captured_frame_vis, ls_local_maxima2_2);
				drawAllSTIPs(captured_frame_vis, ls_local_maxima4_2);
				drawAllSTIPs(captured_frame_vis, ls_local_maxima8_2);
				drawAllSTIPs(captured_frame_vis, ls_local_maxima16_2);
				drawAllSTIPs(captured_frame_vis, ls_local_maxima32_2);
				drawAllSTIPs(captured_frame_vis, ls_local_maxima64_2);
				drawAllSTIPs(captured_frame_vis, ls_local_maxima128_2);
				drawAllSTIPs(captured_frame_vis, ls_local_maxima256_2);
				drawAllSTIPs(captured_frame_vis, ls_local_maxima512_2);
				drawAllSTIPs(captured_frame_vis, ls_local_maxima2_4);
				drawAllSTIPs(captured_frame_vis, ls_local_maxima4_4);
				drawAllSTIPs(captured_frame_vis, ls_local_maxima8_4);
				drawAllSTIPs(captured_frame_vis, ls_local_maxima16_4);
				drawAllSTIPs(captured_frame_vis, ls_local_maxima32_4);
				drawAllSTIPs(captured_frame_vis, ls_local_maxima64_4);
				drawAllSTIPs(captured_frame_vis, ls_local_maxima128_4);
				drawAllSTIPs(captured_frame_vis, ls_local_maxima256_4);
				drawAllSTIPs(captured_frame_vis, ls_local_maxima512_4);
*/

			}
			cv::imshow("Input_STIPS", captured_frame_vis);
			frameno_being_processed++;
		}//if(!paused) end
		if (firstrun == true)
		{
			firstrun = false;
			paused = true;
			std::cout<<"Execution is paused. Press p to toggle pause. Press Esc to exit at any time."<<std::endl;
		}

		char c = (char)cv::waitKey(20);
		if(c == 27) //Esc
			break;
		switch(c)
		{
			case 'p':
				paused = !paused;
				std::cout<<"Paused = "<<paused<<std::endl;
				break;
			default:
				break;
		}
#endif

	}
	myfile.close();
	return 0; //success
}


stipdetector::stipdetector(const tip::params& stip_parameters){
	// TODO Auto-generated constructor stub

#ifdef TIME_MEASUREMENTS
    ETSstart = clock();
    ETSend = clock();
#endif

	k = stip_parameters.kappa;
	s = 2; //scaling factor sigma_i_sq = s * sigma_l_sq; tau_i_sq = s * tau_l_sq
	threshold = stip_parameters.thresh;
	frameno_starting_position = stip_parameters.startingpostition;
	frameno_being_processed = stip_parameters.startingpostition;
	xy_localmaxneighbourhood = 2;
	t_localmaxneighbourhood = 1;

	//internal resize factor for gaussian filtering
	//resize_factor = sqrt(sigma_sq); //resize_factor = sigma !!!
	resize_factor = 1.0;

	sigma_l = sqrt(stip_parameters.sigma_sq) / resize_factor; //sigma_r should be 1 if rs factor = sigma
	sigma_l_sq = round(sigma_l*sigma_l);
	sigma_i_sq = s * sigma_l_sq;

	tau_l_sq = stip_parameters.tau_sq;
	tau_l = sqrt(tau_l_sq);
	tau_i_sq = s * tau_l_sq;

	height_resized = cv::saturate_cast<int>(stip_parameters.frame_height/resize_factor);
	width_resized = cv::saturate_cast<int>(stip_parameters.frame_width/resize_factor);
    SPATIAL_GAUSSIAN_KERNEL_SIZE = cv::saturate_cast<int>(3*sigma_l)+1;
    TEMPORAL_GAUSSIAN_KERNEL_SIZE = 7;//TODO kept constant to keep delay constant for different tau

    SPATIAL_RESIZED_GAUSSIAN_KERNEL_SIZE = SPATIAL_GAUSSIAN_KERNEL_SIZE * s + 1;
   
    TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE = TEMPORAL_GAUSSIAN_KERNEL_SIZE * s + 1;

   	DERIV_KERNEL_SIZE = 3;//must be odd && >=3
   	center_frame_index = (DERIV_KERNEL_SIZE-1)/2; //temporal

	HARRIS_MAT_BUFFER_SIZE = 2 * t_localmaxneighbourhood + 1;
	delay = (TEMPORAL_GAUSSIAN_KERNEL_SIZE-1)/2  + (DERIV_KERNEL_SIZE-1)/2 + (TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE-1)/2 + (HARRIS_MAT_BUFFER_SIZE-1)/2;
	
	//Gaussian kernel for spatial convolution, normalized to sum(all elements) = 1
	mt_1D_Gauss_kernel_sigma_l_sq = cv::getGaussianKernel(SPATIAL_GAUSSIAN_KERNEL_SIZE, sigma_l, CV_32F).t();
	//dispMat(mt_1D_Gauss_kernel_sigma_l_sq, "mt_1D_Gauss_kernel_sigma_l_sq");

	//Gaussian kernel for temporal convolution, normalized to sum(all elements) = 1
	mt_1D_Gauss_kernel_tau_l_sq = cv::getGaussianKernel(TEMPORAL_GAUSSIAN_KERNEL_SIZE, tau_l, CV_32F).t();
	//dispMat(mt_1D_Gauss_kernel_tau_l_sq, "mt_1D_Gauss_kernel_tau_l_sq");

	//Scaled Gaussian kernel for spatial convolution, normalized to sum(all elements) = 1
	mt_1D_Gauss_kernel_sigma_i_sq = cv::getGaussianKernel(SPATIAL_RESIZED_GAUSSIAN_KERNEL_SIZE, sqrt(sigma_i_sq), CV_32F).t();
	//dispMat(mt_1D_Gauss_kernel_sigma_i_sq, "mt_1D_Gauss_kernel_sigma_i_sq");

	//Scaled Gaussian kernel for temporal convolution, normalized to sum(all elements) = 1
	mt_1D_Gauss_kernel_tau_i_sq = cv::getGaussianKernel(TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE, sqrt(tau_i_sq), CV_32F).t();
	//dispMat(mt_1D_Gauss_kernel_tau_i_sq, "mt_1D_Gauss_kernel_tau_i_sq");

	//Partial normalized derivative kernel (Sobel...?)
	getDerivKernels(mt_sob_deriv_kernel_x, mt_sob_deriv_kernel_y, 1, 1, DERIV_KERNEL_SIZE, true, CV_32F);
	//printMat(mt_sob_deriv_kernel_x, "mt_sob_deriv_kernel_x");
	//printMat(mt_sob_deriv_kernel_y, "mt_sob_deriv_kernel_y");
	mt_deriv_kernel_x = mt_sob_deriv_kernel_x.t();
	mt_deriv_kernel_y = mt_deriv_kernel_x.t();

	//Set capacity of all circular buffers
	cb_input_frame.set_capacity(TEMPORAL_GAUSSIAN_KERNEL_SIZE);
	cb_input_frame_flt.set_capacity(DERIV_KERNEL_SIZE);
	cb_LxLx.set_capacity(TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE);
	cb_LyLy.set_capacity(TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE);
	cb_LtLt.set_capacity(TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE);
	cb_LxLy.set_capacity(TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE);
	cb_LxLt.set_capacity(TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE);
	cb_LyLt.set_capacity(TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE);
	cb_H.set_capacity(HARRIS_MAT_BUFFER_SIZE);


	//Fill all buffers with empty matrixes
	for (unsigned int i = 0; i < TEMPORAL_GAUSSIAN_KERNEL_SIZE; i++)
	{
		cv::Mat mt_empty_frame = cv::Mat::zeros(height_resized, width_resized, CV_32F);
		cb_input_frame.push_back(mt_empty_frame);
	}

	for (unsigned int i = 0; i < DERIV_KERNEL_SIZE; i++)
	{
		cv::Mat mt_empty_frame = cv::Mat::zeros(height_resized, width_resized, CV_32F);
		cb_input_frame_flt.push_back(mt_empty_frame);
	}

	for (unsigned int i = 0; i < TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE; i++)
	{
		cv::Mat mt_empty_frame = cv::Mat::zeros(height_resized, width_resized, CV_32F);
		cb_LxLx.push_back(mt_empty_frame);
		cb_LyLy.push_back(mt_empty_frame);
		cb_LtLt.push_back(mt_empty_frame);
		cb_LxLy.push_back(mt_empty_frame);
		cb_LxLt.push_back(mt_empty_frame);
		cb_LyLt.push_back(mt_empty_frame);
	}

	for (unsigned int i = 0; i < HARRIS_MAT_BUFFER_SIZE; i++)
	{
		cv::Mat mt_empty_frame = cv::Mat::zeros(stip_parameters.frame_height, stip_parameters.frame_width, CV_32F); //not resized frame size because it is being resized back
		cb_H.push_back(mt_empty_frame);
	}


	//FOR PROCESS FRAME FAST
	mt_current_frame_flt_temp = cv::Mat::zeros( height_resized, width_resized, CV_32F );
	mt_LxLx_filtered_temporally = cv::Mat::zeros( height_resized, width_resized, CV_32F );
	mt_LyLy_filtered_temporally = cv::Mat::zeros( height_resized, width_resized, CV_32F );
	mt_LtLt_filtered_temporally = cv::Mat::zeros( height_resized, width_resized, CV_32F );
	mt_LxLy_filtered_temporally = cv::Mat::zeros( height_resized, width_resized, CV_32F );
	mt_LxLt_filtered_temporally = cv::Mat::zeros( height_resized, width_resized, CV_32F );
	mt_LyLt_filtered_temporally = cv::Mat::zeros( height_resized, width_resized, CV_32F );
	filtersizewatch = 0;

}

stipdetector::~stipdetector() {
	// TODO Auto-generated destructor stub
}


std::vector<ip::interestpoint> stipdetector::processFrame(const cv::Mat& mt_source, std::vector<ip::interestpoint>& ls_IPs){

	//Clear list of local maxima from previous frame
	ls_IPs.clear();

	//Preprocessing
	cv::Mat captured_frame_gray;
	cvtColor(mt_source, captured_frame_gray, CV_BGR2GRAY);

	//CAUTION !! MUST BE HERE BECAUSE OTHERWISE ONLY THE SAME FRAME
	//GETS INTO THE BUFFER !!!
	cv::Mat mt_input_frame;
	captured_frame_gray.convertTo(mt_input_frame, CV_32F, 1.0/255.0);
	//dispMatProp(mt_input_frame, "mt_input_frame");
	//dispMinMax(mt_input_frame, "mt_input_frame");

	if (resize_factor != 1.0)
	{
		//decimate the image for faster convolution with gaussian kernel
		cv::Mat mt_captured_frame_resized;
		cv::resize(mt_input_frame, mt_captured_frame_resized, cv::Size(), (1.0/resize_factor), (1.0/resize_factor), CV_INTER_AREA);
		//dispMatProp(mt_captured_frame_resized, "mt_captured_frame_resized");
		assert(mt_captured_frame_resized.cols == width_resized);
		assert(mt_captured_frame_resized.rows == height_resized);
		//add to input buffer
		cb_input_frame.push_back(mt_captured_frame_resized);
//		#ifdef SHOWVIDEO
//		showFloatMat(mt_captured_frame_resized,"mt_captured_frame_resized");
//		#endif
	}
	else
	{
		//add to input buffer
		cb_input_frame.push_back(mt_input_frame);
		//#ifdef SHOWVIDEO
		//showFloatMat(mt_input_frame,"mt_captured_frame");
		//#endif
	}

	#ifdef TIME_MEASUREMENTS
		ETSstart = clock();
	#endif

	//STEP 1: Apply scale-space representation
	//a) Temporal convolution with the Gaussian
	cv::Mat mt_current_frame_flt_temp;
	temporalConvolution(cb_input_frame, mt_current_frame_flt_temp, mt_1D_Gauss_kernel_tau_l_sq);
//	#ifdef SHOWVIDEO
//		showFloatMat(mt_current_frame_flt_temp,"flt_temp");
//	#endif

	//b) Spatial convolution with Gaussian to the above result
	cv::Mat mt_current_frame_flt_temp_spac;
	cv::sepFilter2D(mt_current_frame_flt_temp, mt_current_frame_flt_temp_spac, -1, mt_1D_Gauss_kernel_sigma_l_sq, mt_1D_Gauss_kernel_sigma_l_sq);
	//#ifdef SHOWVIDEO
	//	showFloatMat(mt_current_frame_flt_temp_spac, "flt_temp_spac");
	//#endif

	//c) Add to filtered input buffer
	cb_input_frame_flt.push_back(mt_current_frame_flt_temp_spac);

	//STEP 2: Calculate partial derivatives
	//a) Partial derivative temporal dimension
	temporalConvolution(cb_input_frame_flt, mt_Lt, mt_deriv_kernel_x);

	//#ifdef TIME_MEASUREMENTS
	//	ETSend = clock();
	//	std::cout << "STEP 1, scale-space representation takes:" << (double(ETSend //- ETSstart) / CLOCKS_PER_SEC) << "seconds." << '\n';
	//#endif

//	#ifdef TIME_MEASUREMENTS
//		ETSstart = clock();
//	#endif
	//b) Partial derivative spatial dimension x (uses horizontal kernel)
	cv::filter2D(cb_input_frame_flt[(DERIV_KERNEL_SIZE-1)/2], mt_Lx, CV_32F, mt_deriv_kernel_x);

	//b) Partial Derivative spatial dimension y
	cv::filter2D(cb_input_frame_flt[(DERIV_KERNEL_SIZE-1)/2], mt_Ly, CV_32F, mt_deriv_kernel_y);

	//c) Create spatio-temporal second-moment matrix
	cv::Mat mt_LxLx = mt_Lx.mul(mt_Lx);
	cv::Mat mt_LyLy = mt_Ly.mul(mt_Ly);
	cv::Mat mt_LtLt = mt_Lt.mul(mt_Lt);
	cv::Mat mt_LxLy = mt_Lx.mul(mt_Ly);
	cv::Mat mt_LxLt = mt_Lx.mul(mt_Lt);
	cv::Mat mt_LyLt = mt_Ly.mul(mt_Lt);

//	#ifdef TIME_MEASUREMENTS
//		ETSend = clock();
//		std::cout << "STEP 2, partial derivative takes: " << (double(ETSend - ETSstart) / CLOCKS_PER_SEC) << "seconds." << '\n';
//	#endif

	//d) push to buffer
	cb_LxLx.push_back(mt_LxLx);
	cb_LyLy.push_back(mt_LyLy);
	cb_LtLt.push_back(mt_LtLt);
	cb_LxLy.push_back(mt_LxLy);
	cb_LxLt.push_back(mt_LxLt);
	cb_LyLt.push_back(mt_LyLt);

//	#ifdef TIME_MEASUREMENTS
//		ETSstart = clock();
//	#endif

	//a) Temporal convolution with the scaled Gaussian
	temporalConvolution(cb_LxLx, mt_LxLx_filtered_temporally, mt_1D_Gauss_kernel_tau_i_sq);
	//b) Spatial convolution with the scaled Gaussian to the above result
	cv::sepFilter2D(mt_LxLx_filtered_temporally, mt_LxLx_filtered_temporally_and_spacially, -1, mt_1D_Gauss_kernel_sigma_i_sq, mt_1D_Gauss_kernel_sigma_i_sq);

	//a) Temporal convolution with the scaled Gaussian
	temporalConvolution(cb_LyLy, mt_LyLy_filtered_temporally, mt_1D_Gauss_kernel_tau_i_sq);
	//b) Spatial convolution with the scaled Gaussian to the above result
	cv::sepFilter2D(mt_LyLy_filtered_temporally, mt_LyLy_filtered_temporally_and_spacially, -1, mt_1D_Gauss_kernel_sigma_i_sq, mt_1D_Gauss_kernel_sigma_i_sq);

	//a) Temporal convolution with the scaled Gaussian
	temporalConvolution(cb_LtLt, mt_LtLt_filtered_temporally, mt_1D_Gauss_kernel_tau_i_sq);
	//b) Spatial convolution with the scaled Gaussian to the above result
	cv::sepFilter2D(mt_LtLt_filtered_temporally, mt_LtLt_filtered_temporally_and_spacially, -1, mt_1D_Gauss_kernel_sigma_i_sq, mt_1D_Gauss_kernel_sigma_i_sq);

	//a) Temporal convolution with the scaled Gaussian
	temporalConvolution(cb_LxLy, mt_LxLy_filtered_temporally, mt_1D_Gauss_kernel_tau_i_sq);
	//b) Spatial convolution with the scaled Gaussian to the above result
	cv::sepFilter2D(mt_LxLy_filtered_temporally, mt_LxLy_filtered_temporally_and_spacially, -1, mt_1D_Gauss_kernel_sigma_i_sq, mt_1D_Gauss_kernel_sigma_i_sq);

	//a) Temporal convolution with the scaled Gaussian
	temporalConvolution(cb_LxLt, mt_LxLt_filtered_temporally, mt_1D_Gauss_kernel_tau_i_sq);
	//b) Spatial convolution with the scaled Gaussian to the above result
	cv::sepFilter2D(mt_LxLt_filtered_temporally, mt_LxLt_filtered_temporally_and_spacially, -1, mt_1D_Gauss_kernel_sigma_i_sq, mt_1D_Gauss_kernel_sigma_i_sq);

	//a) Temporal convolution with the scaled Gaussian
	temporalConvolution(cb_LyLt, mt_LyLt_filtered_temporally, mt_1D_Gauss_kernel_tau_i_sq);
	//b) Spatial convolution with the scaled Gaussian to the above result
	cv::sepFilter2D(mt_LyLt_filtered_temporally, mt_LyLt_filtered_temporally_and_spacially, -1, mt_1D_Gauss_kernel_sigma_i_sq, mt_1D_Gauss_kernel_sigma_i_sq);

//	#ifdef TIME_MEASUREMENTS
//		ETSend = clock();
//		std::cout << "STEP 3, Second time filtering takes: " << (double(ETSend - ETSstart) / CLOCKS_PER_SEC) << "seconds." << '\n';
//	#endif

//	#ifdef TIME_MEASUREMENTS
//		ETSstart = clock();
//	#endif

	//Step 3 calculate Harris matrix
	//H=det(mu)- k*trace^3(mu);
	//a) det(mu) :
	cv::multiply(mt_LxLx_filtered_temporally_and_spacially, mt_LyLy_filtered_temporally_and_spacially, mt_tmp1);
	cv::multiply(mt_LtLt_filtered_temporally_and_spacially, mt_tmp1, mt_tmp1);
	cv::multiply(mt_LxLy_filtered_temporally_and_spacially, mt_LyLt_filtered_temporally_and_spacially, mt_tmp2);
	cv::multiply(mt_LxLt_filtered_temporally_and_spacially, mt_tmp2, mt_tmp2, 2);
	cv::add(mt_tmp1, mt_tmp2, mt_tmp1);
	cv::multiply(mt_LxLt_filtered_temporally_and_spacially, mt_LxLt_filtered_temporally_and_spacially, mt_tmp2);
	cv::multiply(mt_LyLy_filtered_temporally_and_spacially, mt_tmp2, mt_tmp2);
	cv::subtract(mt_tmp1, mt_tmp2, mt_tmp1);
	cv::multiply(mt_LyLt_filtered_temporally_and_spacially, mt_LyLt_filtered_temporally_and_spacially, mt_tmp2);
	cv::multiply(mt_LxLx_filtered_temporally_and_spacially, mt_tmp2, mt_tmp2);
	cv::subtract(mt_tmp1, mt_tmp2, mt_tmp1);
	cv::multiply(mt_LxLy_filtered_temporally_and_spacially, mt_LxLy_filtered_temporally_and_spacially, mt_tmp2);
	cv::multiply(mt_LtLt_filtered_temporally_and_spacially, mt_tmp2, mt_tmp2);
	cv::subtract(mt_tmp1, mt_tmp2, mt_tmp1); //mt_tmp 1 = det(mu)

//	#ifdef TIME_MEASUREMENTS
//		ETSend = clock();
//		std::cout << "det(mu) calculation takes: " << (double(ETSend - ETSstart) / CLOCKS_PER_SEC) << "seconds." << '\n';
	#endif

//	#ifdef SHOWVIDEO
//		showFloatMat(mt_tmp1, "det(mu)");
		//dispMinMax(mt_tmp1,"det(mu)");
//	#endif

//	#ifdef TIME_MEASUREMENTS
//		ETSstart = clock();//
//	#endif

	//b) trace(mu)^3
	cv::add(mt_LxLx_filtered_temporally_and_spacially, mt_LyLy_filtered_temporally_and_spacially, mt_tmp2);
	cv::add(mt_LtLt_filtered_temporally_and_spacially, mt_tmp2, mt_tmp2);
	cv::pow(mt_tmp2, 3, mt_tmp2); //mt_tmp2 = trace(mu)^3

//	#ifdef TIME_MEASUREMENTS
//		ETSend = clock();
//		std::cout << "trace^2(mu) calculation takes: " << (double(ETSend - ETSstart) / CLOCKS_PER_SEC) << "seconds." << '\n';
//	#endif

//	#ifdef SHOWVIDEO
//		showFloatMat(mt_tmp2, "trace^3(mu)");
		//dispMinMax(mt_tmp2,"trace^3(mu)");
//	#endif

	//c) H=det(mu)- k*trace^3(mu);
	cv::Mat mt_H;
	cv::scaleAdd(mt_tmp2, -k, mt_tmp1, mt_H);

	if (resize_factor != 1.0)
	{
		//scale back up to original resolution
		cv::Mat mt_H_origres;
		cv::resize(mt_H, mt_H_origres, cv::Size(), resize_factor, resize_factor, CV_INTER_AREA);
		//add to buffer
		cb_H.push_back(mt_H_origres);
	//	#ifdef SHOWVIDEO
	//		showFloatMat(mt_H_origres, "mt_H_origres");
	//		dispMinMax(mt_H_origres, "mt_H_origres");
	//	#endif
	}
	else
	{
		//add to buffer
		cb_H.push_back(mt_H);
	//	#ifdef SHOWVIDEO
	//		showFloatMat(mt_H, "mt_H");
	//		dispMinMax(mt_H, "mt_H");
	//	#endif
	}


	//Detect local positive maxima
	double corresponfing_frameno = frameno_being_processed - delay;
	findLocalMaxima3DThreshold(cb_H, ls_IPs, xy_localmaxneighbourhood, t_localmaxneighbourhood, corresponfing_frameno, threshold);

	//delete STIPS in first frames because of initialization
	if ((frameno_being_processed - delay) < 15)
	{
		ls_IPs.clear();
	}
 
/*
#ifdef SHOWVIDEO
		//draw cross
		cv::Mat mt_src = cb_H[t_localmaxneighbourhood];
		cv::Mat mt_tmp;
		cv::normalize(mt_src, mt_tmp, 0, 255, cv::NORM_MINMAX);
		mt_tmp.convertTo(mt_tmp, CV_8U);
		cv::cvtColor(mt_tmp, mt_tmp, CV_GRAY2BGR);
		for(std::vector<ip::interestpoint>::iterator it = ls_IPs.begin(); it != ls_IPs.end();)
		{
			drawCross(mt_tmp, (*it).x, (*it).y, 4);
			it++;
		}

		cv::namedWindow( "pointinH", cv::WINDOW_NORMAL );
		cv::resizeWindow("pointinH", 250, 200);
		cv::imshow( "pointinH", mt_tmp );
#endif
*/

	frameno_being_processed++;
	return ls_IPs;
}

int stipdetector::processFrameTemporal(const cv::Mat& mt_source, std::vector<ip::interestpoint>& ls_IPs, double& corresponfing_frameno){

	//Clear list of local maxima from previous frame
	ls_IPs.clear();

	//Preprocessing
	cv::Mat captured_frame_gray;
	cvtColor(mt_source, captured_frame_gray, CV_BGR2GRAY);

	//Convert to floating point from 0 to 1
	//CAUTION !! MUST BE HERE BECAUSE OTHERWISE ONLY THE SAME FRAME
	//GETS INTO THE BUFFER !!!
	cv::Mat mt_input_frame;
	captured_frame_gray.convertTo(mt_input_frame, CV_32F, 1.0/255.0);
	//dispMatProp(mt_input_frame, "mt_input_frame");
	//dispMinMax(mt_input_frame, "mt_input_frame");

	if (resize_factor != 1.0)
	{
		//decimate the image for faster convolution with gaussian kernel
		cv::Mat mt_captured_frame_resized;
		cv::resize(mt_input_frame, mt_captured_frame_resized, cv::Size(), (1.0/resize_factor), (1.0/resize_factor), CV_INTER_AREA);
		//dispMatProp(mt_captured_frame_resized, "mt_captured_frame_resized");
		assert(mt_captured_frame_resized.cols == width_resized);
		assert(mt_captured_frame_resized.rows == height_resized);
		//add to input buffer
		cb_input_frame.push_back(mt_captured_frame_resized);
//		#ifdef SHOWVIDEO
//		showFloatMat(mt_captured_frame_resized,"mt_captured_frame_resized");
//		#endif
	}
	else
	{
		//add to input buffer
		cb_input_frame.push_back(mt_input_frame);
//		#ifdef SHOWVIDEO
//		showFloatMat(mt_input_frame,"mt_captured_frame");
//		#endif
	}

	//#ifdef TIME_MEASUREMENTS
	//	ETSstart = clock();
	//#endif

	//STEP 1: Apply scale-space representation
	//a) Temporal convolution with the Gaussian
	cv::Mat mt_current_frame_flt_temp;
	temporalConvolution(cb_input_frame, mt_current_frame_flt_temp, mt_1D_Gauss_kernel_tau_l_sq);
//	#ifdef SHOWVIDEO
//		showFloatMat(mt_current_frame_flt_temp,"flt_temp");
//	#endif


	//c) Add to filtered input buffer
	cb_input_frame_flt.push_back(mt_current_frame_flt_temp);

	//STEP 2: Calculate partial derivatives
	//a) Partial derivative temporal dimension
	temporalConvolution(cb_input_frame_flt, mt_Lt, mt_deriv_kernel_x);

	//#ifdef TIME_MEASUREMENTS
	//	ETSend = clock();
	//	std::cout << "STEP 1, scale-space representation takes:" << (double(ETSend - ETSstart) / CLOCKS_PER_SEC) << "seconds." << '\n';
	//#endif


	//c) Create spatio-temporal second-moment matrix
	cv::Mat mt_LtLt = mt_Lt.mul(mt_Lt); //t

	//d) push to buffer
	cb_LtLt.push_back(mt_LtLt);

	//a) Temporal convolution with the scaled Gaussian
	temporalConvolution(cb_LtLt, mt_LtLt_filtered_temporally, mt_1D_Gauss_kernel_tau_i_sq);

	//c) H=det(mu)- k*trace^3(mu);
	cv::Mat mt_H;
	//cv::scaleAdd(mt_tmp2, -0.05, mt_tmp1, mt_H);

	mt_H = mt_LtLt_filtered_temporally;

	if (resize_factor != 1.0)
	{
		//scale back up to original resolution
		cv::Mat mt_H_origres;
		cv::resize(mt_H, mt_H_origres, cv::Size(), resize_factor, resize_factor, CV_INTER_AREA);
		//add to buffer
		cb_H.push_back(mt_H_origres);
//		#ifdef SHOWVIDEO
//			showFloatMat(mt_H_origres, "mt_H_origres");
//			dispMinMax(mt_H_origres, "mt_H_origres");
//		#endif
	}
	else
	{
		//add to buffer
		cb_H.push_back(mt_H);
//		#ifdef SHOWVIDEO
//			showFloatMat(mt_H, "mt_H");
//			dispMinMax(mt_H, "mt_H");
		//#endif
	}


	//Detect local positive maxima
	corresponfing_frameno = frameno_being_processed - delay;
	findLocalMaxima3DThreshold(cb_H, ls_IPs, xy_localmaxneighbourhood, t_localmaxneighbourhood, corresponfing_frameno , threshold);

	//delete STIPS in first frames because of initialization etc.
	if ((corresponfing_frameno) < (frameno_starting_position + delay)) //TODO is delay ok ?????
	{
		ls_IPs.clear();
	}
/*
#ifdef SHOWVIDEO
		//draw cross
		cv::Mat mt_src = cb_H[t_localmaxneighbourhood];
		cv::Mat mt_tmp;
		cv::normalize(mt_src, mt_tmp, 0, 255, cv::NORM_MINMAX);
		mt_tmp.convertTo(mt_tmp, CV_8U);
		cv::cvtColor(mt_tmp, mt_tmp, CV_GRAY2BGR);
		for(std::vector<ip::interestpoint>::iterator it = ls_IPs.begin(); it != ls_IPs.end();)
		{
			drawCross(mt_tmp, (*it).x, (*it).y, 4);
			//drawSTIP(mt_tmp, round((*it).x/resize_factor), ((*it).y/resize_factor), sigma_l_sq);
			it++;
		}

		cv::namedWindow( "pointinH", cv::WINDOW_NORMAL );
		cv::resizeWindow("pointinH", 250, 200);
		cv::imshow( "pointinH", mt_tmp );
#endif
*/

	frameno_being_processed++;
	return 0;

}

std::vector<ip::interestpoint> stipdetector::processFrameFast(const cv::Mat& mt_source, std::vector<ip::interestpoint>& ls_IPs){
	//processFrameFast is using BOX filter instead of Gaussian filter

	//Clear list of local maxima from previous frame
	ls_IPs.clear();

	//Preprocessing
	cv::Mat captured_frame_gray;
	cvtColor(mt_source, captured_frame_gray, CV_BGR2GRAY);

	//CAUTION !! MUST BE HERE BECAUSE OTHERWISE ONLY THE SAME FRAME
	//GETS INTO THE BUFFER !!!
	cv::Mat mt_input_frame;
	captured_frame_gray.convertTo(mt_input_frame, CV_32F, 1.0/255.0);


	//add to input buffer
	cb_input_frame.push_back(mt_input_frame);

	//STEP 1: Apply scale-space representation
	//a) Temporal convolution with the Gaussian
	//cv::Mat mt_current_frame_flt_temp;
	//temporalConvolution(cb_input_frame, mt_current_frame_flt_temp, mt_1D_Gauss_kernel_tau_l_sq);
	temporalBoxFilter(cb_input_frame, mt_current_frame_flt_temp, TEMPORAL_GAUSSIAN_KERNEL_SIZE, 1);
	//mt_current_frame_flt_temp += cb_input_frame.front()/(float)TEMPORAL_GAUSSIAN_KERNEL_SIZE; //add the latest
	//mt_current_frame_flt_temp -= cb_input_frame.back()/(float)TEMPORAL_GAUSSIAN_KERNEL_SIZE; //remove the earliest

//#ifdef SHOWVIDEO
//	showFloatMat(mt_current_frame_flt_temp,"flt_temp");
//#endif

	//b) Spatial convolution with Gaussian to the above result
	//Box filter !!
	cv::Mat mt_current_frame_flt_temp_spac(mt_current_frame_flt_temp.rows, mt_current_frame_flt_temp.cols, mt_current_frame_flt_temp.depth());
	cv::Mat mt_current_frame_flt_temp_spac_tmp(mt_current_frame_flt_temp.rows, mt_current_frame_flt_temp.cols, mt_current_frame_flt_temp.depth());
	rowBoxFilter(mt_current_frame_flt_temp, mt_current_frame_flt_temp_spac_tmp, SPATIAL_GAUSSIAN_KERNEL_SIZE, 1);
	colBoxFilter(mt_current_frame_flt_temp_spac_tmp, mt_current_frame_flt_temp_spac, SPATIAL_GAUSSIAN_KERNEL_SIZE, 1);

//#ifdef SHOWVIDEO
//	showFloatMat(mt_current_frame_flt_temp_spac, "flt_temp_spac");
//#endif

	//c) Add to filtered input buffer
	cb_input_frame_flt.push_back(mt_current_frame_flt_temp_spac);

	//STEP 2: Calculate partial derivatives
	//a) Partial derivative temporal dimension
	temporalConvolution(cb_input_frame_flt, mt_Lt, mt_deriv_kernel_x);

	//b) Partial derivative spatial dimension x (uses horizontal kernel)
	cv::filter2D(cb_input_frame_flt[(DERIV_KERNEL_SIZE-1)/2], mt_Lx, CV_32F, mt_deriv_kernel_x);

	//TODO The function filter2D does actually compute correlation, not the convolution:
	//That is, the kernel is not mirrored around the anchor point. If you need a real
	//convolution, flip the kernel using flip() and set the new anchor to
	//(kernel.cols - anchor.x - 1, kernel.rows - anchor.y - 1).
	//
	//The function uses the DFT-based algorithm in case of sufficiently large kernels
	//(~``11 x 11`` or larger) and the direct algorithm (that uses the engine retrieved
	//by createLinearFilter() ) for small kernels.

	//b) Partial Derivative spatial dimension y
	cv::filter2D(cb_input_frame_flt[(DERIV_KERNEL_SIZE-1)/2], mt_Ly, CV_32F, mt_deriv_kernel_y);

	//c) Create spatio-temporal second-moment matrix and push to buffer
	cv::Mat mt_LxLx = mt_Lx.mul(mt_Lx);
	cb_LxLx.push_back(mt_LxLx);
	cv::Mat mt_LyLy = mt_Ly.mul(mt_Ly);
	cb_LyLy.push_back(mt_LyLy);
	cv::Mat mt_LtLt = mt_Lt.mul(mt_Lt);
	cb_LtLt.push_back(mt_LtLt);
	cv::Mat mt_LxLy = mt_Lx.mul(mt_Ly);
	cb_LxLy.push_back(mt_LxLy);
	cv::Mat mt_LxLt = mt_Lx.mul(mt_Lt);
	cb_LxLt.push_back(mt_LxLt);
	cv::Mat mt_LyLt = mt_Ly.mul(mt_Lt);
	cb_LyLt.push_back(mt_LyLt);

	//a) Temporal convolution with the scaled Gaussian
	temporalBoxFilter(cb_LxLx, mt_LxLx_filtered_temporally, TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE, 1);
	//mt_LxLx_filtered_temporally += cb_LxLx.front()/(float)TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE; //add the latest
	//mt_LxLx_filtered_temporally -= cb_LxLx.back()/(float)TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE; //remove the earliest


	//b) Spatial convolution with the scaled Gaussian to the above result
	cv::Mat mt_LxLx_filtered_temporally_and_spacially(mt_LxLx_filtered_temporally.rows, mt_LxLx_filtered_temporally.cols, mt_LxLx_filtered_temporally.depth());
	cv::Mat mt_LxLx_filtered_temporally_and_spacially_tmp(mt_LxLx_filtered_temporally.rows, mt_LxLx_filtered_temporally.cols, mt_LxLx_filtered_temporally.depth());
	rowBoxFilter(mt_LxLx_filtered_temporally, mt_LxLx_filtered_temporally_and_spacially_tmp, SPATIAL_RESIZED_GAUSSIAN_KERNEL_SIZE, 1);
	colBoxFilter(mt_LxLx_filtered_temporally_and_spacially_tmp, mt_LxLx_filtered_temporally_and_spacially, SPATIAL_RESIZED_GAUSSIAN_KERNEL_SIZE, 1);

	//a) Temporal convolution with the scaled Gaussian
	temporalBoxFilter(cb_LyLy, mt_LyLy_filtered_temporally, TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE, 1);

	//mt_LyLy_filtered_temporally += cb_LyLy.front()/(float)TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE; //add the latest
	//mt_LyLy_filtered_temporally -= cb_LyLy.back()/(float)TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE; //remove the earliest

	//b) Spatial convolution with the scaled Gaussian to the above result
	cv::Mat mt_LyLy_filtered_temporally_and_spacially(mt_LyLy_filtered_temporally.rows, mt_LyLy_filtered_temporally.cols, mt_LyLy_filtered_temporally.depth());
	cv::Mat mt_LyLy_filtered_temporally_and_spacially_tmp(mt_LyLy_filtered_temporally.rows, mt_LyLy_filtered_temporally.cols, mt_LyLy_filtered_temporally.depth());
	rowBoxFilter(mt_LxLx_filtered_temporally, mt_LyLy_filtered_temporally_and_spacially_tmp, SPATIAL_RESIZED_GAUSSIAN_KERNEL_SIZE, 1);
	colBoxFilter(mt_LyLy_filtered_temporally_and_spacially_tmp, mt_LyLy_filtered_temporally_and_spacially, SPATIAL_RESIZED_GAUSSIAN_KERNEL_SIZE, 1);

	//a) Temporal convolution with the scaled Gaussian
	temporalBoxFilter(cb_LtLt, mt_LtLt_filtered_temporally, TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE, 1);
	//mt_LtLt_filtered_temporally += cb_LtLt.front()/(float)TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE; //add the latest
	//mt_LtLt_filtered_temporally -= cb_LtLt.back()/(float)TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE; //remove the earliest

	//b) Spatial convolution with the scaled Gaussian to the above result
	cv::Mat mt_LtLt_filtered_temporally_and_spacially(mt_LtLt_filtered_temporally.rows, mt_LtLt_filtered_temporally.cols, mt_LtLt_filtered_temporally.depth());
	cv::Mat mt_LtLt_filtered_temporally_and_spacially_tmp(mt_LtLt_filtered_temporally.rows, mt_LtLt_filtered_temporally.cols, mt_LtLt_filtered_temporally.depth());
	rowBoxFilter(mt_LtLt_filtered_temporally, mt_LtLt_filtered_temporally_and_spacially_tmp, SPATIAL_RESIZED_GAUSSIAN_KERNEL_SIZE, 1);
	colBoxFilter(mt_LtLt_filtered_temporally_and_spacially_tmp, mt_LtLt_filtered_temporally_and_spacially, SPATIAL_RESIZED_GAUSSIAN_KERNEL_SIZE, 1);

	//a) Temporal convolution with the scaled Gaussian
	temporalBoxFilter(cb_LxLy, mt_LxLy_filtered_temporally, TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE, 1);
	//mt_LxLy_filtered_temporally += cb_LxLy.front()/(float)TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE; //add the latest
	//mt_LxLy_filtered_temporally -= cb_LxLy.back()/(float)TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE; //remove the earliest
	//b) Spatial convolution with the scaled Gaussian to the above result
	cv::Mat mt_LxLy_filtered_temporally_and_spacially(mt_LxLy_filtered_temporally.rows, mt_LxLy_filtered_temporally.cols, mt_LxLy_filtered_temporally.depth());
	cv::Mat mt_LxLy_filtered_temporally_and_spacially_tmp(mt_LxLy_filtered_temporally.rows, mt_LxLy_filtered_temporally.cols, mt_LxLy_filtered_temporally.depth());
	rowBoxFilter(mt_LxLy_filtered_temporally, mt_LxLy_filtered_temporally_and_spacially_tmp, SPATIAL_RESIZED_GAUSSIAN_KERNEL_SIZE, 1);
	colBoxFilter(mt_LxLy_filtered_temporally_and_spacially_tmp, mt_LxLy_filtered_temporally_and_spacially, SPATIAL_RESIZED_GAUSSIAN_KERNEL_SIZE, 1);

	//a) Temporal convolution with the scaled Gaussian
	temporalBoxFilter(cb_LxLt, mt_LxLt_filtered_temporally, TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE, 1);
	//mt_LxLt_filtered_temporally += cb_LxLt.front()/(float)TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE; //add the latest
	//mt_LxLt_filtered_temporally -= cb_LxLt.back()/(float)TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE; //remove the earliest
	//b) Spatial convolution with the scaled Gaussian to the above result
	cv::Mat mt_LxLt_filtered_temporally_and_spacially(mt_LxLt_filtered_temporally.rows, mt_LxLt_filtered_temporally.cols, mt_LxLt_filtered_temporally.depth());
	cv::Mat mt_LxLt_filtered_temporally_and_spacially_tmp(mt_LxLt_filtered_temporally.rows, mt_LxLt_filtered_temporally.cols, mt_LxLt_filtered_temporally.depth());
	rowBoxFilter(mt_LxLt_filtered_temporally, mt_LxLt_filtered_temporally_and_spacially_tmp, SPATIAL_RESIZED_GAUSSIAN_KERNEL_SIZE, 1);
	colBoxFilter(mt_LxLt_filtered_temporally_and_spacially_tmp, mt_LxLt_filtered_temporally_and_spacially, SPATIAL_RESIZED_GAUSSIAN_KERNEL_SIZE, 1);

	//a) Temporal convolution with the scaled Gaussian
	temporalBoxFilter(cb_LyLt, mt_LyLt_filtered_temporally, TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE, 1);
	//mt_LyLt_filtered_temporally += cb_LyLt.front()/(float)TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE; //add the latest
	//mt_LyLt_filtered_temporally -= cb_LyLt.back()/(float)TEMPORAL_RESIZED_GAUSSIAN_KERNEL_SIZE; //remove the earliest
	//b) Spatial convolution with the scaled Gaussian to the above result
	cv::Mat mt_LyLt_filtered_temporally_and_spacially(mt_LyLt_filtered_temporally.rows, mt_LyLt_filtered_temporally.cols, mt_LyLt_filtered_temporally.depth());
	cv::Mat mt_LyLt_filtered_temporally_and_spacially_tmp(mt_LyLt_filtered_temporally.rows, mt_LyLt_filtered_temporally.cols, mt_LyLt_filtered_temporally.depth());
	rowBoxFilter(mt_LyLt_filtered_temporally, mt_LyLt_filtered_temporally_and_spacially_tmp, SPATIAL_RESIZED_GAUSSIAN_KERNEL_SIZE, 1);
	colBoxFilter(mt_LyLt_filtered_temporally_and_spacially_tmp, mt_LyLt_filtered_temporally_and_spacially, SPATIAL_RESIZED_GAUSSIAN_KERNEL_SIZE, 1);

	//Step 3 calc Harris matrix
	//H=det(mu)- k*trace^3(mu);
	//a) det(mu) :
	cv::multiply(mt_LxLx_filtered_temporally_and_spacially, mt_LyLy_filtered_temporally_and_spacially, mt_tmp1);
	cv::multiply(mt_LtLt_filtered_temporally_and_spacially, mt_tmp1, mt_tmp1);
	cv::multiply(mt_LxLy_filtered_temporally_and_spacially, mt_LyLt_filtered_temporally_and_spacially, mt_tmp2);
	cv::multiply(mt_LxLt_filtered_temporally_and_spacially, mt_tmp2, mt_tmp2, 2);
	cv::add(mt_tmp1, mt_tmp2, mt_tmp1);
	cv::multiply(mt_LxLt_filtered_temporally_and_spacially, mt_LxLt_filtered_temporally_and_spacially, mt_tmp2);
	cv::multiply(mt_LyLy_filtered_temporally_and_spacially, mt_tmp2, mt_tmp2);
	cv::subtract(mt_tmp1, mt_tmp2, mt_tmp1);
	cv::multiply(mt_LyLt_filtered_temporally_and_spacially, mt_LyLt_filtered_temporally_and_spacially, mt_tmp2);
	cv::multiply(mt_LxLx_filtered_temporally_and_spacially, mt_tmp2, mt_tmp2);
	cv::subtract(mt_tmp1, mt_tmp2, mt_tmp1);
	cv::multiply(mt_LxLy_filtered_temporally_and_spacially, mt_LxLy_filtered_temporally_and_spacially, mt_tmp2);
	cv::multiply(mt_LtLt_filtered_temporally_and_spacially, mt_tmp2, mt_tmp2);
	cv::subtract(mt_tmp1, mt_tmp2, mt_tmp1); //mt_tmp 1 = det(mu)


//#ifdef SHOWVIDEO
//	showFloatMat(mt_tmp1, "det(mu)");
	//dispMinMax(mt_tmp1,"det(mu)");
//#endif

	//b) trace(mu)^3
	cv::add(mt_LxLx_filtered_temporally_and_spacially, mt_LyLy_filtered_temporally_and_spacially, mt_tmp2);
	cv::add(mt_LtLt_filtered_temporally_and_spacially, mt_tmp2, mt_tmp2);
	cv::pow(mt_tmp2, 3, mt_tmp2); //mt_tmp2 = trace(mu)^3

//#ifdef SHOWVIDEO
//	showFloatMat(mt_tmp2, "trace^3(mu)");
	//dispMinMax(mt_tmp2,"trace^3(mu)");
//#endif

	//c) H=det(mu)- k*trace^3(mu);
	cv::Mat mt_H;
	cv::scaleAdd(mt_tmp2, -k, mt_tmp1, mt_H);
	cb_H.push_back(mt_H);

//#ifdef SHOWVIDEO
//	showFloatMat(mt_H, "mt_H");
//	dispMinMax(mt_H, "mt_H");
//#endif

	//Detect local positive maxima
	double corresponfing_frameno = frameno_being_processed - delay;
	findLocalMaxima3DThreshold(cb_H, ls_IPs, xy_localmaxneighbourhood, t_localmaxneighbourhood, corresponfing_frameno , threshold);

/*
#ifdef SHOWVIDEO
		//draw cross
		cv::Mat mt_src = cb_H[t_localmaxneighbourhood];
		cv::Mat mt_tmp;
		cv::normalize(mt_src, mt_tmp, 0, 255, cv::NORM_MINMAX);
		mt_tmp.convertTo(mt_tmp, CV_8U);
		cv::cvtColor(mt_tmp, mt_tmp, CV_GRAY2BGR);
		for(std::vector<ip::interestpoint>::iterator it = ls_IPs.begin(); it != ls_IPs.end();)
		{
			drawCross(mt_tmp, (*it).x, (*it).y, 4);
			//drawSTIP(mt_tmp, round((*it).x/resize_factor), ((*it).y/resize_factor), sigma_l_sq);
			it++;
		}

		cv::namedWindow( "pointinH", cv::WINDOW_NORMAL );
		cv::resizeWindow("pointinH", 250, 200);
		cv::imshow( "pointinH", mt_tmp );
#endif
*/

	frameno_being_processed++;
	return ls_IPs;
}
cv::Mat stipdetector::returnGaussianKernel(int ksize, float sigma_sq)
{
	//anchor point is center
	//check size is odd
	assert(ksize % 2 == 1);
	assert(ksize > 3);
	cv::Mat mt_kernel(1, ksize, CV_32F);

	int jmin = - (ksize - 1)/2;
	int jmax = (ksize - 1)/2;
	for (int j = jmin; j <= jmax; j++)
	{
		mt_kernel.at<float>(0, j + jmax) = 1/(2*CV_PI*sigma_sq)*std::exp(-j*j/2*sigma_sq);
	}
	return mt_kernel;
}

void stipdetector::temporalConvolution(boost::circular_buffer<cv::Mat>& cb_src, cv::Mat& mt_dst, const cv::Mat& mt_kernel)
{
	//Central frame is frame in the middle e.g. frame 1 in 3 sized buffer
	//mt_kernel must be a row vector
	//dispMatdim(mt_kernel, "mt_kernel");
	//Assert size of mask equals size of buffer
	assert(cb_src.capacity() == (unsigned int)mt_kernel.cols);
	// Check if input_buffer_raw is empty.
	assert(!cb_src.empty());
	assert(cb_src.size() == (unsigned int)mt_kernel.cols);
	//check kernel size
	assert(mt_kernel.cols >=3);
	//TODO if kernel is smaller than buffer pad kernel with zeros ???

	//Restructure the frames in the buffer into one matrix where
	//each frame is one row
	cv::Mat mt_frames_in_rows(cb_src.size(), (cb_src[0].rows * cb_src[0].cols), cb_src[0].type());
	//dispMatProp(mt_frames_in_rows, "mt_frames_in_rows");
	for(unsigned int i=0; i < cb_src.size(); i++)
	{
	     (cb_src[i].reshape(0,1)).copyTo(mt_frames_in_rows.row(i));
	}

//	dispMatDim(mt_frames_in_rows, "mt_frames_in_rows");
//	dispMatProp(mt_frames_in_rows, "mt_frames_in_rows");

	//Matrix multiplication replaces convolution:
	//https://en.wikipedia.org/wiki/Convolution
	//Thus when g has finite support in the set { − M , − M + 1 , … , M − 1 , M }
	//(representing, for instance, a finite impulse response),
	//a finite summation may be used
	//also commutativity is used f * g = g * f
	cv::Mat mt_rowmat_B;
	mt_rowmat_B = mt_kernel * mt_frames_in_rows;

//	dispMatDim(mt_rowmat_B, "mt_rowmat_B");
//	dispMatProp(mt_rowmat_B, "mt_rowmat_B");
	mt_dst = mt_rowmat_B.reshape(0, cb_src[0].rows);
}


void stipdetector::temporalBoxFilter(boost::circular_buffer<cv::Mat>& cb_src, cv::Mat& mt_dst, const int& boxkernelsize, const int& n_iter)
{
	//Central frame is frame in the middle e.g. frame 1 in 3 sized buffer
	mt_dst = cv::Mat::zeros( cb_src[0].rows, cb_src[0].cols, CV_32F );
	for(unsigned int i=0; i < cb_src.size(); i++)
	{
		mt_dst += cb_src[i];
	}

	mt_dst = mt_dst/(float)boxkernelsize;

}


void stipdetector::rowBoxFilter(cv::Mat& mt_src, cv::Mat& mt_dst, const int& boxkernelsize, const int& n_iter)
{
	//this is a little bit faster than rowBoxFilter
	assert(boxkernelsize % 2 == 1); // is odd
	assert(mt_src.size() == mt_dst.size() && mt_src.depth() == mt_dst.depth());

	int border = (boxkernelsize-1)/2;
	float boxkernelsum;

	// Constructs a larger image to fit both the image and the border
	cv::Mat mt_buf(mt_src.rows, mt_src.cols + border*2, mt_src.depth());
	cv::Mat mt_src_temp;
	mt_src.copyTo(mt_src_temp);
	for (int i = 0; i < n_iter; ++i)// iterations
	{
		cv::copyMakeBorder(mt_src_temp, mt_buf, 0, 0, border, border, cv::BORDER_REPLICATE); //Border left right only!!
		for (int row = 0; row < mt_buf.rows; ++row)
		{
			boxkernelsum = 0;
			for (int col = 0; col < boxkernelsize; ++col)
				boxkernelsum += mt_buf.at<float>(row, col);
			mt_dst.at<float>(row, 0) = boxkernelsum/boxkernelsize;

			for (int col = 1; col < mt_dst.cols; ++col)
			{
				boxkernelsum = boxkernelsum - mt_buf.at<float>(row, col-1) + mt_buf.at<float>(row, col+boxkernelsize-1);
				mt_dst.at<float>(row, col) = boxkernelsum/boxkernelsize;
			}
		}
		if(i < n_iter - 1)
			mt_dst.copyTo(mt_src_temp);
	}
}

void stipdetector::colBoxFilter(cv::Mat& mt_src, cv::Mat& mt_dst, const int& boxkernelsize, const int& n_iter)
{
	//this is faster than colBoxFilter !!
	assert(boxkernelsize % 2 == 1); //is odd
	assert(mt_src.size() == mt_dst.size() && mt_src.depth() == mt_dst.depth());

	int border = (boxkernelsize-1)/2;
	float boxkernelsum;

	////BORDER top bottom only!!
	//// constructs a larger image to fit both the image and the border
	cv::Mat mt_buf(mt_src.rows + border*2, mt_src.cols, mt_src.depth());
	cv::Mat mt_src_temp;
	mt_src.copyTo(mt_src_temp);
	for (int i = 0; i < n_iter; ++i)// iterations
	{
		cv::copyMakeBorder(mt_src_temp, mt_buf, border, border, 0, 0, cv::BORDER_REPLICATE);
		for (int col = 0; col < mt_buf.cols; ++col)
		{
			boxkernelsum = 0;
			for (int row = 0; row < boxkernelsize; ++row)
				boxkernelsum += mt_buf.at<float>(row, col);
			mt_dst.at<float>(0, col) = boxkernelsum/boxkernelsize;

			for (int row = 1; row < mt_dst.rows; ++row)
			{
				boxkernelsum = boxkernelsum - mt_buf.at<float>(row-1, col) + mt_buf.at<float>(row+boxkernelsize-1, col);
				mt_dst.at<float>(row, col) = boxkernelsum/boxkernelsize;
			}
		}
		if(i < n_iter - 1)
			mt_dst.copyTo(mt_src_temp);
	}
}

void stipdetector::findLocalMaxima3D(boost::circular_buffer<cv::Mat>& cb_src, std::vector<ip::interestpoint>& ls_local_maxima, int xy_neighbourhood, int t_neighbourhood, int frameno_being_processed)
{
	//TODO CAUTION !! it does NOT search on the border of the 3D area !!
	//only in the region minus the neighbourhood

	//TODO all checks
	// Check if buffer is empty.
	assert(!cb_src.empty());
	assert(cb_src.size() == (unsigned int)(2*t_neighbourhood+1));

	//TODO threshold ?? based on histogram
	//should find most common value range and threshold by making
	//everything outside that range = 0 ??? what if negative valuse exist?
	//then checking in the loop below if center is 0 then it is for sure not a maximum....

	int min_t = t_neighbourhood;
	int max_t = cb_src.size() - t_neighbourhood - 1;
	int min_x = xy_neighbourhood;
	int min_y = xy_neighbourhood;
	int max_x = cb_src[0].cols - xy_neighbourhood - 1;
	int max_y = cb_src[0].rows - xy_neighbourhood - 1;

	int min_t_nb, max_t_nb, min_x_nb, min_y_nb, max_x_nb, max_y_nb;

	//float value_at_center = cb_src[i].at<float>(y,x);
	float value_at_center;
	float value_at_neighbour;
	bool value_at_center_is_max;
	for (int t = min_t; t <= max_t; t++)
	{
		for (int y = min_y; y <= max_y; y++)
		{
			for (int x = min_x; x <= max_x; x++)
			{

				value_at_center = cb_src[t].at<float>(y,x);
				if(value_at_center <= 0)//TODO is this assumption ok ???
					continue;

				//here check if already thresholded based on histogram
				//otherwise continue
				value_at_center_is_max = true; //we don't know jet but assume
				//create neigbourhood:
				min_t_nb = t - t_neighbourhood;
				max_t_nb = t + t_neighbourhood;
				min_x_nb = x - xy_neighbourhood;
				max_x_nb = x + xy_neighbourhood;
				min_y_nb = y - xy_neighbourhood;
				max_y_nb = y + xy_neighbourhood;

				for (int t_nb = min_t_nb; t_nb <= max_t_nb; t_nb++)
				{
					for (int y_nb = min_y_nb; y_nb <= max_y_nb; y_nb++)
					{
						for (int x_nb = min_x_nb; x_nb <= max_x_nb; x_nb++)
						{
							//skip center point
							if (t_nb == t && x_nb == x && y_nb == y)
								continue;

							value_at_neighbour = cb_src[t_nb].at<float>(y_nb, x_nb);
							
							if (value_at_center <= value_at_neighbour)
							{
								value_at_center_is_max = false;
								break;
							}
						}
						if (value_at_center_is_max == false)
							break;
					}
					if (value_at_center_is_max == false)
						break;
				}
				//if this is true it means that the neighbourhood was checked and nothing was
				//greater than or equal to the value at center
				if (value_at_center_is_max == true)
				{
					//found a local maximum...
					ip::interestpoint tmp_localmax;
					tmp_localmax.x = x;
					tmp_localmax.y = y;
					//tmp_localmax.t = t;
					tmp_localmax.H_value = value_at_center;
					tmp_localmax.frameno = frameno_being_processed;
					tmp_localmax.sigma_l_sq = sigma_l_sq;
					tmp_localmax.tau_l_sq = tau_l_sq;
					tmp_localmax.s = s;
					ls_local_maxima.push_back(tmp_localmax);
					//skip the remaining neigbourhood, at least in one positive direction
					//this should save time
					x += xy_neighbourhood;
				}

			}
		}

	}


}

void stipdetector::findLocalMaxima3DThreshold(boost::circular_buffer<cv::Mat>& cb_src, std::vector<ip::interestpoint>& ls_local_maxima, int& xy_neighbourhood, int& t_neighbourhood, double& frameno_being_processed, float& thresh)
{
	//TODO CAUTION !! it does NOT search on the border of the 3D area !!
	//only in the region minus the neighbourhood

	//TODO all checks
	// Check if buffer is empty.
	assert(!cb_src.empty());
	assert(cb_src.size() == (unsigned int)(2*t_neighbourhood+1));

	//TODO threshold ?? based on histogram
	//should find most common value range and threshold by making
	//everything outside that range = 0 ??? what if negative valuse exist?
	//then checking in the loop below if center is 0 then it is for sure not a maximum....

	int min_t = t_neighbourhood;
	int max_t = cb_src.size() - t_neighbourhood - 1;
	int min_x = xy_neighbourhood;
	int min_y = xy_neighbourhood;
	int max_x = cb_src[0].cols - xy_neighbourhood - 1;
	int max_y = cb_src[0].rows - xy_neighbourhood - 1;

	int min_t_nb, max_t_nb, min_x_nb, min_y_nb, max_x_nb, max_y_nb;

	//float value_at_center = cb_src[i].at<float>(y,x);
	float value_at_center;
	float value_at_neighbour;
	bool value_at_center_is_max;
	for (int t = min_t; t <= max_t; t++)
	{
		for (int y = min_y; y <= max_y; y++)
		{
			for (int x = min_x; x <= max_x; x++)
			{

				value_at_center = cb_src[t].at<float>(y,x);
				if(value_at_center <= 0)//TODO is this assumption ok ???
					continue;

				//here check if already thresholded based on histogram
				//otherwise continue
				value_at_center_is_max = true; //we don't know jet but assume
				//create neigbourhood:
				min_t_nb = t - t_neighbourhood;
				max_t_nb = t + t_neighbourhood;
				min_x_nb = x - xy_neighbourhood;
				max_x_nb = x + xy_neighbourhood;
				min_y_nb = y - xy_neighbourhood;
				max_y_nb = y + xy_neighbourhood;

				for (int t_nb = min_t_nb; t_nb <= max_t_nb; t_nb++)
				{
					for (int y_nb = min_y_nb; y_nb <= max_y_nb; y_nb++)
					{
						for (int x_nb = min_x_nb; x_nb <= max_x_nb; x_nb++)
						{
							//skip center point
							if (t_nb == t && x_nb == x && y_nb == y)
								continue;

							value_at_neighbour = cb_src[t_nb].at<float>(y_nb, x_nb);
							
							if (value_at_center <= value_at_neighbour)
							{
								value_at_center_is_max = false;
								break;
							}
						}
						if (value_at_center_is_max == false)
							break;
					}
					if (value_at_center_is_max == false)
						break;
				}
				//if this is true it means that the neighbourhood was checked and nothing was
				//greater than or equal to the value at center
				if (value_at_center_is_max == true)
				{
					if(value_at_center > thresh)
					{
						ip::interestpoint tmp_localmax;
						tmp_localmax.x = x;
						tmp_localmax.y = y;
						tmp_localmax.H_value = value_at_center;
						tmp_localmax.frameno = frameno_being_processed;
						tmp_localmax.sigma_l_sq = sigma_l_sq;
						tmp_localmax.tau_l_sq = tau_l_sq;
						tmp_localmax.s = s;
						ls_local_maxima.push_back(tmp_localmax);
					}
					//skip the remaining neigbourhood, at least in one positive direction
					//this should save time
					x += xy_neighbourhood;
				}
			}
		}
	}

}


} /* namespace stip */
