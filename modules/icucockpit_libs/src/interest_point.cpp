/*
 * interestpoint.cpp
 *
 *  Created on: Aug 3, 2017
 *      Author: aed
 */

#include "interest_point.h"

//ADDITIONAL includes needed by interest_point.cpp
//OpenCV
#include "opencv2/imgproc.hpp"
//#include <opencv2/highgui.hpp>


namespace ip {

/*
void drawIP(cv::Mat& mt_dst, int x, int y, float size_r)
{
	//cv::Scalar color = CV_RGB(255, 0, 0);
	cv::Scalar color = CV_RGB(0, 255, 0);
	int linetype = 8;// CV_AA;
	int radius = round(size_r);
	//line(Mat& img, Point pt1, Point pt2, const Scalar& color, int thickness=1, int lineType=8, int shift=0)
	line(mt_dst, cv::Point(x-radius,y), cv::Point(x+radius,y), color , 1,linetype, 0 );
	line(mt_dst, cv::Point(x,y-radius), cv::Point(x,y+radius), color , 1,linetype, 0 );
	//circle(Mat& img, Point center, int radius, const Scalar& color, int thickness=1, int lineType=8, int shift=0)
	circle(mt_dst, cv::Point(x, y), radius, color, 1, linetype, 0);
}

void drawIP(cv::Mat& mt_dst, int x, int y, float size_r, cv::Scalar color)
{
	//cv::Scalar color = CV_RGB(255, 0, 0);
	//cv::Scalar color = CV_RGB(0, 255, 0);
	int linetype = 8;// CV_AA;
	int radius = round(size_r);
	//line(Mat& img, Point pt1, Point pt2, const Scalar& color, int thickness=1, int lineType=8, int shift=0)
	line(mt_dst, cv::Point(x-radius,y), cv::Point(x+radius,y), color , 1,linetype, 0 );
	line(mt_dst, cv::Point(x,y-radius), cv::Point(x,y+radius), color , 1,linetype, 0 );
	//circle(Mat& img, Point center, int radius, const Scalar& color, int thickness=1, int lineType=8, int shift=0)
	circle(mt_dst, cv::Point(x, y), radius, color, 1, linetype, 0);
}

void drawIP(cv::Mat& mt_dst, int x, int y, float size_r, cv::Scalar color, const float& scale)
{
	int linetype = 8;// CV_AA;
	//int radius = cv::saturate_cast<int>(2*sqrt(sigma_l_sq)*scale);
	int radius = round(size_r * scale);
	x = cv::saturate_cast<int>(x * scale);
	y = cv::saturate_cast<int>(y * scale);
	//line(Mat& img, Point pt1, Point pt2, const Scalar& color, int thickness=1, int lineType=8, int shift=0)
	line(mt_dst, cv::Point(x-radius,y), cv::Point(x+radius,y), color , 1, linetype, 0 );
	line(mt_dst, cv::Point(x,y-radius), cv::Point(x,y+radius), color , 1, linetype, 0 );
	//circle(Mat& img, Point center, int radius, const Scalar& color, int thickness=1, int lineType=8, int shift=0)
	circle(mt_dst, cv::Point(x, y), radius, color, 1, linetype, 0);
}


void drawAllTIPs(cv::Mat& mt_dst, std::vector<interestpoint>& vc_ip, const float& scale)
{
	cv::Scalar blue = CV_RGB(0, 0, 255);
	cv::Scalar green = CV_RGB(0, 255, 0);
	cv::Scalar red = CV_RGB(255, 0, 0);
	//if vc_ip is empty, for loop exits immediately, no assertion necessary
	for(std::vector<interestpoint>::iterator it = vc_ip.begin(); it != vc_ip.end();)
	{
		if((*it).x >=0  && (*it).y >= 0)//... no IP found
		{
			if ((*it).roilabel == "INROI")
				drawIP(mt_dst, (*it).x, (*it).y, (*it).hog_patch_size/2, red, scale);
			else if ((*it).roilabel == "OUTROI")
				drawIP(mt_dst, (*it).x, (*it).y, (*it).hog_patch_size/2, blue, scale);
		}
		it++;
	}
}

void drawAllBRISKs(cv::Mat& mt_dst, std::vector<interestpoint>& vc_ip, const float& scale)
{
	cv::Scalar blue = CV_RGB(0, 0, 255);
	cv::Scalar green = CV_RGB(0, 255, 0);
	cv::Scalar red = CV_RGB(255, 0, 0);
	//if vc_ip is empty, for loop exits immediately, no assertion necessary
	for(std::vector<interestpoint>::iterator it = vc_ip.begin(); it != vc_ip.end();)
	{
		if((*it).x >=0  && (*it).y >= 0)//... no IP found
		{
			if ((*it).roilabel == "INROI")
				drawIP(mt_dst, (*it).x, (*it).y, (*it).brisk_size, red, scale);
			else if ((*it).roilabel == "OUTROI")
				drawIP(mt_dst, (*it).x, (*it).y, (*it).brisk_size, blue, scale);
		}
		it++;
	}
}
*/

void updateTimestampToAllIPs(std::vector<interestpoint>& vc_ip, boost::posix_time::ptime timestamp)
{
	for(std::vector<interestpoint>::iterator it = vc_ip.begin(); it != vc_ip.end(); it++)
		{
			(*it).video_timestamp = timestamp;
		}
}

void appendIPsAndDescrToCsvFile(std::ofstream& ofs_myfile, const std::vector<interestpoint>& vc_ip, std::string& attr_a, std::string& attr_b, std::string& attr_c, std::string& attr_d, std::string& attr_e, std::string& attr_f)
{
	assert(ofs_myfile.is_open());
	//if vc_ip is empty for loop exits immediately, no assertion necessary
	for(std::vector<interestpoint>::const_iterator it = vc_ip.begin(); it != vc_ip.end(); it++)
	{
		if ((*it).HOG_descriptors.empty())
			continue;

		ofs_myfile	<<(*it).y
					<<","<<(*it).x
					<<","<<(*it).frameno
					<<","<<boost::posix_time::to_iso_extended_string((*it).video_timestamp);
		for(std::vector<float>::const_iterator it2 = (*it).HOG_descriptors.begin(); it2 != (*it).HOG_descriptors.end(); it2++)
		{
			ofs_myfile <<","<<(*it2);
		}
		ofs_myfile <<","<<attr_a<<","<<attr_b<<","<<attr_c<<","<<attr_d<<","<<attr_e<<","<<attr_f;
		ofs_myfile <<std::endl;
	}
}

void appendIPsAndHOGDescrToCsvFile(std::ofstream& ofs_myfile, const std::vector<interestpoint>& vc_ip, const std::vector<std::string>& vc_labels)
{
	assert(ofs_myfile.is_open());
	std::vector<std::string> vc_labevc_internal=vc_labels;
	int labevc_max_no = 20;
	std::string no_label = "no_label";

	//deal with missing labels
	int vc_size = vc_labevc_internal.size();
	if (vc_size > labevc_max_no)
		std::cerr<<"WARNING: number of labels in this frame are more than 20, only the first 20 will be saved!"<<std::endl;
	if (vc_size < labevc_max_no)
	{
		for (int i = 0; i < (labevc_max_no - vc_size); i++)
			vc_labevc_internal.push_back(no_label);
	}

	//if vc_ip is empty for loop exits immediately, no assertion necessary
	for(std::vector<interestpoint>::const_iterator it = vc_ip.begin(); it != vc_ip.end(); it++)
	{
//		if ((*it).HOG_descriptors.empty())
//			continue;

		ofs_myfile	<<(*it).y
					<<","<<(*it).x
					<<","<<(*it).frameno
					<<","<<boost::posix_time::to_iso_extended_string((*it).video_timestamp);
		for(std::vector<float>::const_iterator it2 = (*it).HOG_descriptors.begin(); it2 != (*it).HOG_descriptors.end(); it2++)
		{
			ofs_myfile <<","<<(*it2);
		}
		ofs_myfile	<<","<<(*it).roilabel;
		for(std::vector<std::string>::iterator vc_it = vc_labevc_internal.begin(); vc_it != vc_labevc_internal.begin() + labevc_max_no ; vc_it++)
		{
			ofs_myfile <<","<<(*vc_it);
		}
		ofs_myfile<<std::endl;
	}
}

void appendIPsAndBRISKDescrToCsvFile(std::ofstream& ofs_myfile, const std::vector<interestpoint>& vc_ip, const std::vector<std::string>& vc_labels)
{
	assert(ofs_myfile.is_open());
	std::vector<std::string> vc_labevc_internal=vc_labels;
	int labevc_max_no = 20;
	std::string no_label = "no_label";

	//deal with missing labels
	int vc_size = vc_labevc_internal.size();
	if (vc_size > labevc_max_no)
		std::cerr<<"WARNING: number of labels in this frame are more than 20, only the first 20 will be saved!"<<std::endl;
	if (vc_size < labevc_max_no)
	{
		for (int i = 0; i < (labevc_max_no - vc_size); i++)
			vc_labevc_internal.push_back(no_label);
	}

	//if vc_ip is empty for loop exits immediately, no assertion necessary
	for(std::vector<interestpoint>::const_iterator it = vc_ip.begin(); it != vc_ip.end(); it++)
	{
//		if ((*it).BRISK_descriptors.empty())
//			continue;

		ofs_myfile	<<(*it).y
					<<","<<(*it).x
					<<","<<(*it).frameno
					<<","<<boost::posix_time::to_iso_extended_string((*it).video_timestamp);
		for(std::vector<float>::const_iterator it2 = (*it).BRISK_descriptors.begin(); it2 != (*it).BRISK_descriptors.end(); it2++)
		{
			ofs_myfile <<","<<(*it2);
		}
		ofs_myfile	<<","<<(*it).roilabel;
		for(std::vector<std::string>::iterator vc_it = vc_labevc_internal.begin(); vc_it != vc_labevc_internal.begin() + labevc_max_no ; vc_it++)
		{
			ofs_myfile <<","<<(*vc_it);
		}
		ofs_myfile<<std::endl;
	}
}

void appendIPsAndHOGBRISKDescrToCsvFile(std::ofstream& ofs_myfile, const std::vector<interestpoint>& vc_ip, const std::vector<std::string>& vc_labels)
{
	assert(ofs_myfile.is_open());
	std::vector<std::string> vc_labevc_internal=vc_labels;
	int labevc_max_no = 20;
	std::string no_label = "no_label";

	//deal with missing labels
	int vc_size = vc_labevc_internal.size();
	if (vc_size > labevc_max_no)
		std::cerr<<"WARNING: number of labels in this frame are more than 20, only the first 20 will be saved!"<<std::endl;
	if (vc_size < labevc_max_no)
	{
		for (int i = 0; i < (labevc_max_no - vc_size); i++)
			vc_labevc_internal.push_back(no_label);
	}

	//if vc_ip is empty for loop exits immediately, no assertion necessary
	for(std::vector<interestpoint>::const_iterator it = vc_ip.begin(); it != vc_ip.end(); it++)
	{
		ofs_myfile	<<(*it).y
					<<","<<(*it).x
					<<","<<(*it).frameno
					<<","<<boost::posix_time::to_iso_extended_string((*it).video_timestamp);
		for(std::vector<float>::const_iterator it2 = (*it).HOG_descriptors.begin(); it2 != (*it).HOG_descriptors.end(); it2++)
		{
			ofs_myfile <<","<<(*it2);
		}
		for(std::vector<float>::const_iterator it2 = (*it).BRISK_descriptors.begin(); it2 != (*it).BRISK_descriptors.end(); it2++)
		{
			ofs_myfile <<","<<(*it2);
		}
		ofs_myfile	<<","<<(*it).roilabel;
		for(std::vector<std::string>::iterator vc_it = vc_labevc_internal.begin(); vc_it != vc_labevc_internal.begin() + labevc_max_no ; vc_it++)
		{
			ofs_myfile <<","<<(*vc_it);
		}
		ofs_myfile<<std::endl;
	}
}

void convertIPVectorToBRISKKeyPoints(const std::vector<interestpoint>& vc_ip, std::vector<cv::KeyPoint>& vc_kp_out)
{
	vc_kp_out.clear();
	cv::KeyPoint temp_keypoint;
	for(std::vector<interestpoint>::const_iterator it = vc_ip.begin(); it != vc_ip.end(); it++)
	{
		temp_keypoint.pt.x = (float)(*it).x;
		temp_keypoint.pt.y = (float)(*it).y;
		// Size is defined by sigma_l_sq from TIP ??

		temp_keypoint.size = (*it).hog_patch_size/2;
		vc_kp_out.push_back(temp_keypoint);
	}
}

void convertBRISKKeyPointsAndBRISKDescriptorsToIPVector(const std::vector<cv::KeyPoint>& vc_kp_in, const cv::Mat& decr_in, const int& brisk_descr_size, const int& frameno, std::vector<interestpoint>& vc_ip)
{
	if (!vc_ip.empty())
		vc_ip.clear();
	interestpoint ip_temp;
	int idx = 0;
	std::vector<cv::KeyPoint>::const_iterator it = vc_kp_in.begin(), it_end = vc_kp_in.end();
	for( ; it != it_end; ++it, ++idx )
	{
		ip_temp.x = cvRound((*it).pt.x);
		ip_temp.y = cvRound((*it).pt.y);
		ip_temp.frameno = frameno;
		ip_temp.brisk_size = (*it).size;
		decr_in.row(idx).copyTo(ip_temp.BRISK_descriptors); //TODO !!! copies uint8 to float
		vc_ip.push_back(ip_temp);
	}
}

void convertBRISKKeyPointsToIPVector(const std::vector<cv::KeyPoint>& vc_kp_in, const int& frameno, std::vector<interestpoint>& vc_ip_out)
{
	if (!vc_ip_out.empty())
		vc_ip_out.clear();
	interestpoint ip_temp;
	std::vector<cv::KeyPoint>::const_iterator it = vc_kp_in.begin(), it_end = vc_kp_in.end();
	for(; it != it_end; ++it)
	{
		ip_temp.x = cvRound((*it).pt.x);
		ip_temp.y = cvRound((*it).pt.y);
		ip_temp.frameno = frameno;
		ip_temp.brisk_size = (*it).size;
		vc_ip_out.push_back(ip_temp);
	}
}

void addBRISKDescrToIPVector(const std::vector<cv::KeyPoint>& vc_kp_in, const cv::Mat& decr_in, const int& brisk_descr_size, const int& frameno, std::vector<interestpoint>& vc_ip)
{
	std::vector<float> temp_vector;
	std::vector<float> brisk_zero_vector(brisk_descr_size, 0.0); //vector with zeros
	bool match_found;
	for(std::vector<interestpoint>::iterator it_ip = vc_ip.begin() ; it_ip != vc_ip.end(); ++it_ip)
	{
		int idx = 0;
		match_found = false;
		for(std::vector<cv::KeyPoint>::const_iterator it_kp = vc_kp_in.begin(); it_kp != vc_kp_in.end(); ++it_kp, ++idx)
		{
			if (((*it_ip).x == (*it_kp).pt.x) && ((*it_ip).y == (*it_kp).pt.y))
			{
				// match found
				match_found = true;
				decr_in.row(idx).copyTo(temp_vector); //TODO !!! copies uint8 to float
				(*it_ip).BRISK_descriptors = temp_vector;
				break;
			}
		}
		if (match_found == false)
			(*it_ip).BRISK_descriptors = brisk_zero_vector;
	}
}


void checkIfIPsAreInsideRectAndGiveLabel(std::vector<interestpoint>& vc_ip, const cv::Rect& roirect, const std::string& ifinsidelabel, const std::string& ifoutsidelabel)
{
	int xmin = roirect.x;
	int xmax = roirect.x + roirect.width;
	int ymin = roirect.y;
	int ymax = roirect.y + roirect.height;
	for(std::vector<interestpoint>::iterator it = vc_ip.begin(); it != vc_ip.end(); it++)
	{
		if( ((*it).x >= xmin) && ((*it).x < xmax) && ((*it).y >= ymin) && ((*it).y < ymax) )
			(*it).roilabel = ifinsidelabel;
		else
			(*it).roilabel = ifoutsidelabel;
	}
}


void checkIfIPsBRISKEmptyAndAddDummyIP(std::vector<interestpoint>& vc_ip, const int& frameno, const int& descr_size)
{
	if (vc_ip.empty())
	{
		interestpoint ip_temp;
		std::vector<float> temp_vector(descr_size, 0); //vector with zeros
		ip_temp.x = -1;
		ip_temp.y = -1;
		ip_temp.sigma_l_sq = -1;
		ip_temp.frameno = frameno;
		ip_temp.BRISK_descriptors = temp_vector; //TODO !!! copies uint8 to float
		vc_ip.push_back(ip_temp);
	}
}

void checkIfIPsHOGEmptyAndAddDummyIP(std::vector<interestpoint>& vc_ip, const int& frameno, const int& descr_size)
{
	if (vc_ip.empty())
	{
		interestpoint ip_temp;
		std::vector<float> temp_vector(descr_size, 0.0); //vector with zeros
		ip_temp.x = -1;
		ip_temp.y = -1;
		ip_temp.sigma_l_sq = -1;
		ip_temp.frameno = frameno;
		ip_temp.HOG_descriptors = temp_vector; //TODO !!! copies uint8 to float
		vc_ip.push_back(ip_temp);
	}
}

} /* namespace ip */
