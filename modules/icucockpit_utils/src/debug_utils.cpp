/*
 * debug_utils.cpp
 *
 *  Created on: Jul 26, 2017
 *      Author: aed
 */

#include "debug_utils.h"

//Additional includes needed by debug_utils.cpp
//C++
#include <math.h>
#include <iostream>

//OpenCV
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"

void dispMatProp(const cv::Mat& mt_src, const std::string& name)
{
	std::cout << "Matrix name: " << name << std::endl;
	switch (mt_src.depth())
	{
	case 0:
		std::cout << name << " depth: CV_8U" << std::endl;
		break;
	case 1:
		std::cout << name << " depth: CV_8S" << std::endl;
		break;
	case 2:
		std::cout << name << " depth: CV_16U" << std::endl;
		break;
	case 3:
		std::cout << name << " depth: CV_16S" << std::endl;
		break;
	case 4:
		std::cout << name << " depth: CV_32S" << std::endl;
		break;
	case 5:
		std::cout << name << " depth: CV_32F" << std::endl;
		break;
	case 6:
		std::cout << name << " depth: CV_64F" << std::endl;
		break;
	case 7:
		std::cout << name << "Matrix depth: CV_USRTYPE1" << std::endl;
		break;
	default:
		std::cout << name << "Matrix depth: UNKNOWN" << std::endl;
		break;
	}
	std::cout << name << " channels: " << mt_src.channels()  << name << " width (columns): " << mt_src.size().width << std::endl<< name << " height (rows): " << mt_src.size().height << std::endl;
}

void dispMatDim(const cv::Mat& mt_src, const std::string& name)
{
	std::cout<<name<<" rows: "<<mt_src.rows<<name<<" cols: "<<mt_src.cols<<std::endl;
}

void dispMat(const cv::Mat& mt_src, const std::string& name)
{
	std::cout<<name<<" = "<< std::endl <<" "<<mt_src<<std::endl<<std::endl;
}

void dispMinMax(const cv::Mat& mt_src, const std::string& name)
{
	//For finding global min max
	double min, max;
	cv::minMaxLoc(mt_src, &min, &max);
    std::cout<<name<<" min = "<<min<<std::endl<<name<<" max = "<<max<<std::endl;
}

void showFloatMat(const cv::Mat& mt_src, const std::string& name)
{
	cv::Mat mt_tmp;
	cv::normalize(mt_src, mt_tmp, 0, 255, cv::NORM_MINMAX);
	//dispMinMax(mt_tmp, "mt_tmp");
	mt_tmp.convertTo(mt_tmp, CV_8U);
	cv::namedWindow( name, cv::WINDOW_NORMAL );
	cv::resizeWindow(name, 250, 200);
	cv::imshow( name, mt_tmp );
}

void drawCross(cv::Mat& mt_dst, int x, int y, int length)
{
	cv::Scalar color = CV_RGB(255, 0, 0);
	int linetype = 8;// CV_AA;
	line(mt_dst, cv::Point(x-length,y), cv::Point(x+length,y), color , 1,linetype, 0 );
	line(mt_dst, cv::Point(x,y-length), cv::Point(x,y+length), color , 1,linetype, 0 );
}

