/**
    main.cpp
    Purpose: Runs the various video feature extraction methods

    @author M. Pediaditis IBM
    @version 0.1 3/11/17
*/

//TODO 20170522
//normalize output against various sigma_sq and tau_sq so that threshold is applicable => different thresholds defined
//multithreading ?? with pthreads or http://www.boost.org/doc/libs/1_63_0/doc/html/thread.html
//OR: do convolution using FFT...

//#define USE_DLIB_FACEDETECTION
//#define DEBUG_VIDEO_WRITER

// C++
#include <iostream>
// OpenCV
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/features2d.hpp" //for BRISK etc.
// Debug utilities
#include "debug_utils.h"
// Interest points class
#include "interest_point.h"
// STIP detector
#include "stip_detector.h"
// HOG descriptors
#include "hog_descriptor.h"
// Annotations
#include "roi_annotation.h" //class for handling yml annotations
#include "xml_annotation.h" //class for handling xml annotation file

#ifdef USE_DLIB_FACEDETECTION
// dlib face detector
// set dialect flags -mavx for fastest face detection
#include "dlib/opencv.h"
#include "dlib/image_processing/frontal_face_detector.h"
//#include "dlib/image_processing/render_face_detections.h"
#include "dlib/image_processing.h"
//#include "dlib/gui_widgets.h"

static cv::Rect dlibRectangleToCVRect(dlib::rectangle r)
{
	return cv::Rect(cv::Point2i(r.left(), r.top()), cv::Point2i(r.right() + 1, r.bottom() + 1));
}
#endif

int run_fe_for_all_features(cv::CommandLineParser& parser);

//UNUSED, OLD OR UNDER DEVELOPMENT
/*
#define SHOWVIDEO_MAIN
int runTIP_HOG_OR_BRISK_AnnXml(cv::CommandLineParser& parser);
int runBRISK_BRISK_AnnXml(cv::CommandLineParser& parser);
int runDenseHOGExtration(cv::CommandLineParser& parser);
int runBRISKExtration(cv::CommandLineParser& parser, const float& resizefactor);

int runSTIPDetection_V1(cv::CommandLineParser& parser, const float& resizefactor);
*/

static void help()
{
    std::cout
    << "\nThis program does...\n"
    << std::endl;
}

//for CommandLineParser in OpenCV 3:
const cv::String keys =
	"{help h usage ?                  |      |Print help message}"
	"{ivf input_video_file            |      |Input video file for feature extraction}"
	"{ifl input_video_file_list       |      |List (.txt) of all input video files for feature extraction}"
	"{ann annotation_file             |      |Input file (.xml) with annotations for input video file}"
	"{anl annotation_file_list        |      |List (.txt) of all input files (.xml) with annotations}"
	"{lbl video_label                 |      |Label for input video file or whole input file list if no annotation file exists}"
	"{otp output_path                 |<none>|Path to write output file(s) with interest points and descriptors (feature_type.csv will be appended to video file name)}"
	"{tph extract_tip_hog_features    |      |Feature type: tip_hog}"
	"{tpb extract_tip_brisk_features  |      |Feature type: tip_brisk}"
	"{brh extract_brisk_hog_features  |      |Feature type: brisk_hog}"
	"{brb extract_brisk_brisk_features|      |Feature type: brisk_brisk}"
	"{vst video_start_date_time       |      |Video start date time. Format: 2017-10-20 01:36:34.655}"
	"{pst processing_start_date_time  |      |Processing start date time. Format: 2017-10-20 01:36:34.655}"
	"{pet processing_end_date_time    |      |Processing end date time. Format: 2017-10-20 01:36:34.655}"
	"{rix region_of_interest_x        |      |Top left x position of ROI to calculate features from}"
	"{riy region_of_interest_y        |      |Top left y position of ROI to calculate features from}"
	"{riw region_of_interest_w        |      |Width of ROI to calculate features from}"
	"{rih region_of_interest_h        |      |Height of ROI to calculate features from}"
	"{ofs annotation_data_offset      |0     |Annotation data offset in ms with respect to video start date time (t(a) = t(v) + ofs), default = 0}"
	"{vrb verbose                     |      |Print debug information}"
	"{shv show_video                  |      |Show video}"
	"{fdt face_detection_enabled      |      |Enable face detection}"
	"{rsf resize_factor               |1.0   |Resize factor for resizing input video prior to IP detection. Default = 1.0}"
	"{thr tip_threshold               |0.0001|TIP threshold. Default = 0.0001}"
	"{ssq tip_sigma_squared           |2     |TIP sigma_sq. Default = 2}"
	"{tsq tip_tau_squared             |2     |TIP tau_sq. Default = 2}"
    ;

int main(int argc, const char** argv)
{
    cv::CommandLineParser parser(argc, argv, keys);
    if (parser.has("help"))
    {
        help();
        return 0;
    }

    int return_value = 0;
    if (parser.has("extract_tip_hog_features") ||
    	parser.has("extract_tip_brisk_features") ||
		parser.has("extract_brisk_hog_features") ||
		parser.has("extract_brisk_brisk_features"))
    {
	    return_value = run_fe_for_all_features(parser);
    }
    return return_value;
}


/**
    Runs feature extraction using both TIP and BRISK interest points. HOG and BRISK descriptors are calculated for both interest point types

    @param parser, the parser object for getting input data and parameters
    @return success: 0 = ok
*/
int run_fe_for_all_features(cv::CommandLineParser& parser){

	// Verbose
	bool verbose = false;
	if (parser.has("verbose"))
		verbose = true;

	// Show video
	bool show_video = false;
	if (parser.has("show_video"))
		show_video = true;

	// Face detection
	bool face_detection = false;
	if (parser.has("face_detection_enabled"))
	{
		face_detection = true;
		#ifndef USE_DLIB_FACEDETECTION
		std::cerr<<"ERROR: Face detection is not possible. Define USE_DLIB_FACEDETECTION in the preprocessor directives."<<std::endl;
		return -1;
		#endif
	}


	// Video start date time
	boost::posix_time::ptime video_start_posix_date_time;
	std::string str_video_start_date_time = parser.get<std::string>("video_start_date_time");
	// Format check
	assert(str_video_start_date_time[4] == '-');
	assert(str_video_start_date_time[7] == '-');
	assert(str_video_start_date_time[10] == ' ');
	assert(str_video_start_date_time[13] == ':');
	assert(str_video_start_date_time[16] == ':');
	assert(str_video_start_date_time[19] == '.');
	video_start_posix_date_time = boost::posix_time::time_from_string(str_video_start_date_time);
	//std::cout<<"Video start date time (ISO extended): "<< boost::posix_time::to_iso_extended_string(video_start_posix_date_time)<<std::endl;

	// Processing start date time
	std::string str_processing_start_date_time;
	boost::posix_time::ptime processing_start_posix_date_time;
	if (parser.has("processing_start_date_time"))
	{
		// Processing start date time
		str_processing_start_date_time = parser.get<std::string>("processing_start_date_time");
		// Format check
		assert(str_processing_start_date_time[4] == '-');
		assert(str_processing_start_date_time[7] == '-');
		assert(str_processing_start_date_time[10] == ' ');
		assert(str_processing_start_date_time[13] == ':');
		assert(str_processing_start_date_time[16] == ':');
		assert(str_processing_start_date_time[19] == '.');
		processing_start_posix_date_time = boost::posix_time::time_from_string(str_processing_start_date_time);
	}
	else
	{
		processing_start_posix_date_time = video_start_posix_date_time;
	}
	//std::cout<<"Processing start date time (ISO extended): "<< boost::posix_time::to_iso_extended_string(processing_start_posix_date_time)<<std::endl;

	// Processing end date time
	std::string processing_end_date_time;
	boost::posix_time::ptime processing_end_posix_date_time;
	if (parser.has("processing_end_date_time"))
	{
		// Processing start date time
		processing_end_date_time = parser.get<std::string>("processing_end_date_time");
		// Format check
		assert(processing_end_date_time[4] == '-');
		assert(processing_end_date_time[7] == '-');
		assert(processing_end_date_time[10] == ' ');
		assert(processing_end_date_time[13] == ':');
		assert(processing_end_date_time[16] == ':');
		assert(processing_end_date_time[19] == '.');
		processing_end_posix_date_time = boost::posix_time::time_from_string(processing_end_date_time);
	}
	else
	{
		processing_end_posix_date_time = boost::posix_time::pos_infin;
	}

	// Apply ROI for selecting a sub-region inside the video
	cv::Rect subregion_roi;
	bool subregion_roi_needed = false;
	if (parser.has("region_of_interest_x") && parser.has("region_of_interest_y") &&parser.has("region_of_interest_w") &&parser.has("region_of_interest_h"))
	{
		subregion_roi.x = parser.get<int>("region_of_interest_x");
		subregion_roi.y = parser.get<int>("region_of_interest_y");
		subregion_roi.width = parser.get<int>("region_of_interest_w");
		subregion_roi.height = parser.get<int>("region_of_interest_h");
		subregion_roi_needed = true;
	}

	// Annotation data offset (t(a) = t(v) + ofs)
	double annotation_data_offset_msec = parser.get<double>("annotation_data_offset");
	double current_pos_msec;

	// IMPORTANT resize factor of input video to make things faster
	float resize_factor = parser.get<float>("resize_factor");

	// Input video for TIP detection
	std::string in_file_name = parser.get<std::string>("input_video_file");
	cv::VideoCapture cap;
	cap.open(in_file_name, cv::CAP_FFMPEG);
	if( !cap.isOpened() )
	{
		std::cout<< "***Could not initialize capturing...***"<<std::endl<<"ERROR: Unable to open video file."<<std::endl;
		return -1;
	}


	// Allow the provision of a label for the whole video
	std::string video_label;
	bool provided_video_label = false;
	if (parser.has("video_label"))
	{
		provided_video_label = true;
		video_label = parser.get<std::string>("video_label");
	}

	// Set video position for processing start
	boost::posix_time::time_duration video_position_start_time_from_zero;
	video_position_start_time_from_zero = processing_start_posix_date_time - video_start_posix_date_time;
	cap.set(CV_CAP_PROP_POS_MSEC, (double)video_position_start_time_from_zero.total_milliseconds());
	double video_position_start_frame = cap.get(CV_CAP_PROP_POS_FRAMES);

	// For descriptor extraction and visualization
	cv::VideoCapture cap2;
	cap2.open(in_file_name, cv::CAP_FFMPEG);
	if( !cap2.isOpened() )
	{
		std::cout << "***Could not initialize capturing...***"<<std::endl<<"ERROR: Unable to open video file."<<std::endl;
		return -1;
	}
	bool cap2_position_set = false;

	// For writing output data (TIPs, HOGs and BRISKs) to csv
	std::string out_file_path = "./";
	if(parser.has("output_path"))
		out_file_path = parser.get<std::string>("output_path");
	else
	{
		std::cout << "***Parameter output_path must be defined***"<<std::endl;
		return -1;
	}
	std::string name, ext;
	size_t dot = in_file_name.find_last_of(".");
	size_t path_end = in_file_name.find_last_of("/\\");
	if (dot != std::string::npos)
	{
		name = in_file_name.substr(path_end+1, dot);
		ext  = in_file_name.substr(dot, in_file_name.size() - dot);
	}

	// Extract tip_hog features
	bool extract_tip_hog_features = false;
	if (parser.has("extract_tip_hog_features"))
		extract_tip_hog_features = true;

	// Extract tip_brisk features
	bool extract_tip_brisk_features = false;
	if (parser.has("extract_tip_brisk_features"))
		extract_tip_brisk_features = true;

	// Extract brisk_hog features
	bool extract_brisk_hog_features = false;
	if (parser.has("extract_brisk_hog_features"))
		extract_brisk_hog_features = true;

	// Extract brisk_brisk features
	bool extract_brisk_brisk_features = false;
	if (parser.has("extract_brisk_brisk_features"))
		extract_brisk_brisk_features = true;

	std::ofstream myfile_tip_hog;
	if (extract_tip_hog_features)
	{
		std::string out_file_name_tip_hog = out_file_path + name + "_tip_hog.csv";
		myfile_tip_hog.open ((char*)out_file_name_tip_hog.c_str(), std::ios::trunc); //Any current content is discarded
		if (!myfile_tip_hog.is_open())
		{
			std::cerr<<"ERROR: Unable to open output file for writing tip_hog features."<<std::endl;
			return -1;
		}
	}

	std::ofstream myfile_tip_brisk;
	if (extract_tip_brisk_features)
	{
		std::string out_file_name_tip_brisk = out_file_path + name + "_tip_brisk.csv";
		myfile_tip_brisk.open ((char*)out_file_name_tip_brisk.c_str(), std::ios::trunc); //Any current content is discarded
		if (!myfile_tip_brisk.is_open())
		{
			std::cerr<<"ERROR: Unable to open output file for writing tip_brisk features."<<std::endl;
			return -1;
		}
	}

	std::ofstream myfile_brisk_hog;
	if (extract_brisk_hog_features)
	{
		std::string out_file_name_brisk_hog = out_file_path + name + "_brisk_hog.csv";
		myfile_brisk_hog.open ((char*)out_file_name_brisk_hog.c_str(), std::ios::trunc); //Any current content is discarded
		if (!myfile_brisk_hog.is_open())
		{
			std::cerr<<"ERROR: Unable to open output file for writing brisk_hog features."<<std::endl;
			return -1;
		}
	}

	std::ofstream myfile_brisk_brisk;
	if (extract_brisk_brisk_features)
	{
		std::string out_file_name_brisk_brisk = out_file_path + name + "_brisk_brisk.csv";
		myfile_brisk_brisk.open ((char*)out_file_name_brisk_brisk.c_str(), std::ios::trunc); //Any current content is discarded
		if (!myfile_brisk_brisk.is_open())
		{
			std::cerr<<"ERROR: Unable to open output file for writing brisk_brisk features."<<std::endl;
			return -1;
		}
	}

	// Open USZ xml annotation file
	std::string annotation_file_name;
	std::vector<anxml::xml_event> vc_xml_events;
	std::vector<std::string> vc_event_labels;
	int labevc_found;
	bool read_xml = false;
	if (parser.has("annotation_file"))
	{
		read_xml = true;
		annotation_file_name = parser.get<std::string>("annotation_file");
		if(!anxml::loadXmlFile(annotation_file_name, vc_xml_events))
		{
			std::cerr<<"ERROR: Problem in loading xml annotation file."<<std::endl;
			return -1;
		}
	}

	// ***********************************************
	// Initialize TIP detector
	// ***********************************************
	// TODO cap.get(CV_CAP_PROP_FRAME_HEIGHT) doesn't work when gstreamer is used for reading videos
	//int cheight = 1080;
	//int cwidth = 1920;

	int cwidth = cap.get(CV_CAP_PROP_FRAME_WIDTH);
	int cheight = cap.get(CV_CAP_PROP_FRAME_HEIGHT);

	// If a smaller area of the video is being analyzed
	if(subregion_roi_needed)
	{
		assert(subregion_roi.x >= 0);
		assert(subregion_roi.y >= 0);
		assert(subregion_roi.x + subregion_roi.width <= cwidth);
		assert(subregion_roi.y + subregion_roi.height <= cheight);
		cwidth = subregion_roi.width;
		cheight = subregion_roi.height;
	}
	int resized_width = cv::saturate_cast<int>(cwidth*resize_factor);
	int resized_height = cv::saturate_cast<int>(cheight*resize_factor);

	float tip_threshold = parser.get<float>("tip_threshold");
	float tip_sigma_squared = parser.get<float>("tip_sigma_squared");
	float tip_tau_squared = parser.get<float>("tip_tau_squared");

	//TIP parameters
	tip::params tip_params;
	//keep IP if(value_at_center > thresh)
	tip_params.kappa = 0.0005;
	tip_params.thresh = tip_threshold; //KTH
	tip_params.frame_height = resized_height;
	tip_params.frame_width = resized_width;
	tip_params.startingpostition = video_position_start_frame;

	//Instantiate all TIP detectors
	tip_params.sigma_sq = tip_sigma_squared;
	tip_params.tau_sq = tip_tau_squared;
	//tip::stipdetector mytipdetector_00(tip_params);
	//tip_params.sigma_sq = 64;
	//tip_params.tau_sq = 32;
	//tip::stipdetector mytipdetector_01(tip_params);
	//tip_params.sigma_sq = 32;
	//tip_params.tau_sq = 64;
	//tip::stipdetector mytipdetector_02(tip_params);

	//FOR ICU = 4
	//tip_params.sigma_sq = 4;
	//tip_params.tau_sq = 4;

	//FOR ICU = 2
	//tip_params.sigma_sq = 2;
	//tip_params.tau_sq = 2;

	//FOR UCF = 2
	//tip_params.thresh = 0.0001; //KTH
	//tip_params.sigma_sq = 8;
	//tip_params.tau_sq = 4;

	tip::stipdetector mytipdetector_03(tip_params);
	//tip_params.sigma_sq = 252;
	//tip_params.tau_sq = 258;
	//tip::stipdetector mytipdetector_04(tip_params);
	//tip_params.sigma_sq = 254;
	//tip_params.tau_sq = 258;
	//tip::stipdetector mytipdetector_05(tip_params);

	double corresponding_frame;

	//std::vector<ip::interestpoint> vc_tips_01;
	//std::vector<ip::interestpoint> vc_tips_02;
	std::vector<ip::interestpoint> vc_ip_tip_hog_03;
	//std::vector<ip::interestpoint> vc_tips_04;
	//std::vector<ip::interestpoint> vc_tips_05;
	//std::vector<ip::interestpoint> vc_tips_06;


	#ifdef DEBUG_VIDEO_WRITER
	// ***********************************************
	// Video Writer
	// ***********************************************
    cv::Size S = cv::Size(resized_width, resized_height);
	//cv::Size S = cv::Size(resized_height, resized_width); //because it is being rotated
    cv::VideoWriter outputVideo;                                        // Open the output
    std::string OUT_VIDEO_NAME = "out_icu_anfall_4.avi";
    outputVideo.open(OUT_VIDEO_NAME, CV_FOURCC('D', 'I', 'V', 'X'), cap2.get(CV_CAP_PROP_FPS), S, true);

    if (!outputVideo.isOpened())
    {
        std::cout  << "Could not open the output video for write: " << OUT_VIDEO_NAME << std::endl;
        return -1;
    }
	#endif

	// ***********************************************
	// ROI rectangle for special label e.g. face region
	// ***********************************************

	// For resized image (adjust to resized images):
	cv::Rect myROI(0,0,0,0);
	cv::Rect myROI_resized(0,0,0,0);
	// Labels to give if interest point is inside or outside the ROI
	std::string inlabel = "INROI";
	std::string outlabel = "OUTROI";

	#ifdef USE_DLIB_FACEDETECTION
	// ***********************************************
	// dlib face detector
	// ***********************************************
	// Load face detection and pose estimation models.
	dlib::frontal_face_detector f_detector = dlib::get_frontal_face_detector();
	//dlib::shape_predictor pose_model;
	//dlib::deserialize("shape_predictor_68_face_landmarks.dat") >> pose_model;
	#endif

	// ***********************************************
	// HOG descriptors
	// ***********************************************
	int HOG_number_of_bins = 5;
	int hog_descriptor_size = HOG_number_of_bins * 9;

	// ***********************************************
	// Initialize BRISK detector
	// ***********************************************
	// Vector of detected interest points and their descriptors
	std::vector<cv::KeyPoint> vc_keypoints_tip_brisk;
	std::vector<cv::KeyPoint> vc_keypoints_brisk_brisk;
	cv::Mat mt_descriptors_tip_brisk;
	cv::Mat mt_descriptors_brisk_brisk;
	// BRISK descriptors
	cv::Ptr<cv::Feature2D> b;
	b = cv::BRISK::create(50, 2, 1.0f); // cv::BRISK::create (int thresh = 30, int octaves = 3, float patternScale = 1.0f)
	int brisk_descriptor_size = 64;

	std::vector<ip::interestpoint> vc_ip_tip_brisk_03;
	std::vector<ip::interestpoint> vc_ip_brisk_hog;
	std::vector<ip::interestpoint> vc_ip_brisk_brisk;

	// Create windows for visualization
	std::string windowname_tip = "TIP_HOG";
	std::string windowname_brisk = "BRISK_HOG";
	std::string text_for_video = "";
	cv::Scalar color_for_text;

	cv::Mat mt_captured_frame_vis_tip = cv::Mat::zeros(resized_height, resized_width, CV_8UC3);
	cv::Mat mt_captured_frame_vis_brisk = cv::Mat::zeros(resized_height, resized_width, CV_8UC3);

	cv::Mat mt_captured_frame;
	cv::Mat mt_captured_frame_resized;
	cv::Mat mt_captured_frame_descr;
	cv::Mat mt_captured_frame_descr_resized;
	bool paused = false;

	// ***********************************************
	// Start the loop
	// ***********************************************
	for(;;)
	//for(int i = 1; i<=20; i++)
	{
		if(!paused)
		{
			cap >> mt_captured_frame;
			if(mt_captured_frame.empty())
				break;

			//Apply subregion ROI
			if(subregion_roi_needed)
				mt_captured_frame = mt_captured_frame(subregion_roi);

			// Reduce the resolution of the input frame to make things faster
			if (resize_factor != 1.0)
				cv::resize(mt_captured_frame, mt_captured_frame_resized, cv::Size(), resize_factor, resize_factor, CV_INTER_AREA);
			else
				mt_captured_frame_resized = mt_captured_frame;

			//***********************************************
			//TIP DETECTION
			//***********************************************
			// FOR VIZ !!@!@> remove
			mytipdetector_03.processFrameTemporal(mt_captured_frame_resized, vc_ip_tip_hog_03, corresponding_frame);
			//corresponding_frame = cap.get(CV_CAP_PROP_POS_FRAMES)-1;
			if (corresponding_frame >= video_position_start_frame)// Start video capture for descriptor extraction
			{
				if (!cap2_position_set)// This needs to be set only once
				{
					cap2.set(CV_CAP_PROP_POS_FRAMES, corresponding_frame);
					cap2_position_set = true;
				}
				cap2 >> mt_captured_frame_descr;
				if(mt_captured_frame_descr.empty())
					break;

				//Apply subregion ROI
				if(subregion_roi_needed)
					mt_captured_frame_descr = mt_captured_frame_descr(subregion_roi);

				// Reduce the resolution of the input frame to make things faster
				if (resize_factor != 1.0)
					cv::resize(mt_captured_frame_descr, mt_captured_frame_descr_resized, cv::Size(), resize_factor, resize_factor, CV_INTER_AREA);
				else
					mt_captured_frame_descr_resized = mt_captured_frame_descr;

				// Frame position in date time format
				current_pos_msec = cap2.get(CV_CAP_PROP_POS_MSEC);
				boost::posix_time::ptime video_current_position_time(video_start_posix_date_time + boost::posix_time::milliseconds(current_pos_msec));
				boost::posix_time::ptime video_current_position_time_corrected(video_current_position_time + boost::posix_time::milliseconds(annotation_data_offset_msec));
            
				if (video_current_position_time > processing_end_posix_date_time)
					break;

				// Read label data from xml file
				if (read_xml)
				{
					//Search for labels inside the xml annotation data
					labevc_found = anxml::findEvent(video_current_position_time_corrected, vc_xml_events, vc_event_labels);
				}
				else if (provided_video_label)
				{
					vc_event_labels.clear();
					vc_event_labels.push_back(video_label);
				}

				#ifdef USE_DLIB_FACEDETECTION
				//***********************************************
				//FACE DETECTION
				//***********************************************
				if (face_detection)
				{
					dlib::cv_image<dlib::bgr_pixel> cimg(mt_captured_frame_descr_resized);
					// Detect faces
					std::vector<dlib::rectangle> faces = f_detector(cimg);
					myROI = dlibRectangleToCVRect(faces[0]);
				}
				#endif


				myROI_resized.x = cv::saturate_cast<int>(myROI.x*resize_factor);
				myROI_resized.y = cv::saturate_cast<int>(myROI.y*resize_factor);
				myROI_resized.height = cv::saturate_cast<int>(myROI.height*resize_factor);
				myROI_resized.width = cv::saturate_cast<int>(myROI.width*resize_factor);


				//***********************************************
				//DESCRIPTOR EXTRACTION FOR TIP POINTS
				//***********************************************
				if(verbose)
				{
					std::cout << "Descriptor frame: "<<corresponding_frame<<std::endl<<"Size of vc_tips: "<<vc_ip_tip_hog_03.size()<<std::endl;
				}

				// Compute HOG descriptors for each TIP interest point
				// Patchsize is defined by sigma_l_sq from TIP
				hog::compute3x3HOGForAllTIPs(mt_captured_frame_descr_resized, vc_ip_tip_hog_03, HOG_number_of_bins);

				// Compute BRISK descriptors for each TIP interest point
				ip::convertIPVectorToBRISKKeyPoints(vc_ip_tip_hog_03, vc_keypoints_tip_brisk);
				b->compute(mt_captured_frame_descr_resized, vc_keypoints_tip_brisk, mt_descriptors_tip_brisk);
        
				//***********************************************
				//BRISK KEYPOINT DETECTION
				//***********************************************
				b->detect(mt_captured_frame_descr_resized, vc_keypoints_brisk_brisk, cv::Mat());
				// or detect and compute descriptors in one step
				//b->detectAndCompute(mt_captured_frame, cv::Mat(), vc_keypoints_tip_brisk, mt_descriptors_tip_brisk, false);

				//***********************************************
				//DESCRIPTOR EXTRACTION FOR BRISK KEYPOINT
				//***********************************************

				// Compute BRISK descriptors for each BRISK interest point
				b->compute(mt_captured_frame_descr_resized, vc_keypoints_brisk_brisk, mt_descriptors_brisk_brisk);

				// Compute HOG descriptors for each BRISK key point
				ip::convertBRISKKeyPointsToIPVector(vc_keypoints_brisk_brisk, corresponding_frame, vc_ip_brisk_hog);
				hog::compute3x3HOGForAllBRISKIPs(mt_captured_frame_descr_resized, vc_ip_brisk_hog, HOG_number_of_bins);
				if (verbose)
					std::cout << "vc_ip_brisk_hog.size() = "<<vc_ip_brisk_hog.size()<<std::endl;


				// Check if empty and add one dummy IP
				checkIfIPsBRISKEmptyAndAddDummyIP(vc_ip_tip_brisk_03, corresponding_frame, brisk_descriptor_size);
				checkIfIPsHOGEmptyAndAddDummyIP(vc_ip_tip_hog_03, corresponding_frame, hog_descriptor_size);
				checkIfIPsBRISKEmptyAndAddDummyIP(vc_ip_brisk_brisk, corresponding_frame, brisk_descriptor_size);
				checkIfIPsHOGEmptyAndAddDummyIP(vc_ip_brisk_hog, corresponding_frame, hog_descriptor_size);

				// Check if ip is inside or outside the ROI
				ip::checkIfIPsAreInsideRectAndGiveLabel(vc_ip_tip_brisk_03, myROI_resized, inlabel, outlabel);
				ip::checkIfIPsAreInsideRectAndGiveLabel(vc_ip_tip_hog_03, myROI_resized, inlabel, outlabel);
				ip::checkIfIPsAreInsideRectAndGiveLabel(vc_ip_brisk_brisk, myROI_resized, inlabel, outlabel);
				ip::checkIfIPsAreInsideRectAndGiveLabel(vc_ip_brisk_hog, myROI_resized, inlabel, outlabel);

				// Add timestamp information to all IPs
				updateTimestampToAllIPs(vc_ip_tip_brisk_03, video_current_position_time);
				updateTimestampToAllIPs(vc_ip_tip_hog_03, video_current_position_time);
				updateTimestampToAllIPs(vc_ip_brisk_brisk, video_current_position_time);
				updateTimestampToAllIPs(vc_ip_brisk_hog, video_current_position_time);

				// Append interest points, descriptors and labels to .csv file
				if (extract_tip_brisk_features)
					ip::appendIPsAndBRISKDescrToCsvFile(myfile_tip_brisk, vc_ip_tip_brisk_03, vc_event_labels);

				if (extract_tip_hog_features)
					ip::appendIPsAndHOGDescrToCsvFile(myfile_tip_hog, vc_ip_tip_hog_03, vc_event_labels);

				if (extract_brisk_brisk_features)
					ip::appendIPsAndBRISKDescrToCsvFile(myfile_brisk_brisk, vc_ip_brisk_brisk, vc_event_labels);

				if (extract_brisk_hog_features)
					ip::appendIPsAndHOGDescrToCsvFile(myfile_brisk_hog, vc_ip_brisk_hog, vc_event_labels);

			}
		}

	}//for(;;) end
	if (extract_tip_brisk_features)
		myfile_tip_brisk.close();

	if (extract_tip_hog_features)
		myfile_tip_hog.close();

	if (extract_brisk_brisk_features)
		myfile_brisk_brisk.close();

	if (extract_brisk_hog_features)
		myfile_brisk_hog.close();
	return 0;
}

