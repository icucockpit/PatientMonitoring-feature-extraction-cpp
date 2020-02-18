/*
 * hogdescriptor.cpp
 *
 *  Created on: Jul 3, 2017
 *      Author: aed
 */


//#define SHOWVIDEO_HOGDESCRIPTOR

#include "hog_descriptor.h"

//ADDITIONAL includes needed by hogdescriptor.cpp
//C++
#include <iostream>
#include <time.h> //for calculating elapsed time
#include <math.h>

//OpenCV
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/videoio.hpp"
#include "opencv2/objdetect.hpp"

namespace hog {

void testHOG(){

	//create test image
	cv::Mat mt_patch = cv::Mat(cv::Size(54,54), CV_8UC3);
	mt_patch = cv::Scalar(255,255,255);

	cv::circle(mt_patch, cv::Point(27, 27), 20, cv::Scalar(0,0,0), 1);

	std::vector<float> vc_descriptors;
	std::vector<cv::Point> vc_locations;

	hog::compute3x3HOG(mt_patch, vc_descriptors, vc_locations, 20);

}

int testHOGOnVideo(cv::CommandLineParser& parser){

	std::string in_file_name = parser.get<std::string>("input_video_file");
	cv::VideoCapture cap;
	cap.open(in_file_name);

	if( !cap.isOpened() )
	{
		std::cout << "***Could not initialize capturing...***\n";
		std::cerr<<"ERROR: Unable to open video file."<<std::endl;
		return -1;
	}

	//for moving the ROI
	int patchsize = 54;
	int capheight = cap.get(CV_CAP_PROP_FRAME_HEIGHT);
	int capwidth = cap.get(CV_CAP_PROP_FRAME_WIDTH);
	int Dx = 3;
	int Dy = 3;
	int Dx_add = Dx;
	int Dy_add = Dy;
	int border = 40;
	int minx = border;
	int miny = border;
	int maxx = capwidth - patchsize - border;
	int maxy = capheight - patchsize - border;
	int x = minx;
	int y = miny;

	cv::Mat mt_captured_frame0;
	cv::Mat mt_captured_frame;
	bool paused = false;
	//***********************************************
	//Start the loop
	//***********************************************
	for(;;)
	{
		if(!paused)
		{
			cap >> mt_captured_frame0;
			if(mt_captured_frame0.empty())
				break;
			//std::cout<<"captured frame number: "<<cap.get(CV_CAP_PROP_POS_FRAMES)-1<<std::endl;

			//moving patch...
			cv::Rect STIPROI(x, y, patchsize, patchsize);
			if (x <= minx)
				Dx_add = Dx;
			if (y <= minx)
				Dy_add = Dy;
			if (x >= maxx)
				Dx_add = -1 * Dx;
			if (y >= maxy)
				Dy_add = -1 * Dy;
			x += Dx_add;
			y += Dy_add;

			cv::Mat mt_patch = mt_captured_frame0(STIPROI);
			std::vector<float> vc_descriptors;
			std::vector<cv::Point> vc_locations;
			hog::compute3x3HOG(mt_patch, vc_descriptors, vc_locations, 20);

			//imshow( "patch", hog::get3x3HOGVisualization(mt_patch, vc_descriptors, 20));

			//show patch on input video...
			//cv::rectangle(mt_captured_frame0, STIPROI, cv::Scalar(0,0,255), 2);
			//cv::namedWindow("mt_input_frame", cv::WINDOW_NORMAL);
			//cv::resizeWindow("mt_input_frame", 1000, 500);
			//cv::imshow( "mt_input_frame", mt_captured_frame0);

		}

		char c = (char)cv::waitKey(0);
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

	}//for(;;) end
	return 0;
}

int testHOGOnVideo2(cv::CommandLineParser& parser){

	std::string in_file_name = parser.get<std::string>("input_video_file");
	cv::VideoCapture cap;
	cap.open(in_file_name);

	if( !cap.isOpened() )
	{
		std::cerr<<" ***Could not initialize capturing...***\nERROR: Unable to open video file."<<std::endl;
		return -1;
	}

	//for moving the ROI
	int patchsize = 54;
	int capheight = cap.get(CV_CAP_PROP_FRAME_HEIGHT);
	int capwidth = cap.get(CV_CAP_PROP_FRAME_WIDTH);
	int Dx = 3;
	int Dy = 3;
	int Dx_add = Dx;
	int Dy_add = Dy;
	int border = 40;
	int minx = border;
	int miny = border;
	int maxx = capwidth - patchsize - border;
	int maxy = capheight - patchsize - border;
	int x = minx;
	int y = miny;

	cv::Mat mt_captured_frame0;
	cv::Mat mt_captured_frame;
	bool paused = false;
	//***********************************************
	//Start the loop
	//***********************************************
	for(;;)
	{
		if(!paused)
		{
			cap >> mt_captured_frame0;
			if(mt_captured_frame0.empty())
				break;
    
			if (x <= minx)
				Dx_add = Dx;
			if (y <= minx)
				Dy_add = Dy;
			if (x >= maxx)
				Dx_add = -1 * Dx;
			if (y >= maxy)
				Dy_add = -1 * Dy;
			x += Dx_add;
			y += Dy_add;

			std::vector<float> vc_descriptors;
			std::vector<cv::Point> vc_locations;
			hog::compute3x3HOGAtPoint(mt_captured_frame0, vc_descriptors, vc_locations, cv::Point(x,y), patchsize, 20);
			
		}

		char c = (char)cv::waitKey(0);
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

	}//for(;;) end
	return 0;
}

void computeDense3x3HOG(const cv::Mat& mt_src, const int& framenumber, std::vector<ip::interestpoint>& ls_locations, const int& step, const int& patchsize, const int& numbins)
{
	ls_locations.clear();
	ip::interestpoint temp_interestpoint;
	temp_interestpoint.H_value= 0.0;
	temp_interestpoint.frameno = framenumber;
	temp_interestpoint.sigma_l_sq = 4.0; //TODO dummy values for drawing
	temp_interestpoint.tau_l_sq = 4.0; //TODO dummy values for drawing

	int half_patchsize = (patchsize - 1)/2;
	for (int x = half_patchsize; x < mt_src.cols - half_patchsize; x+= step)
	{
		for (int y = half_patchsize; y < mt_src.rows - half_patchsize; y+= step)
			{
				temp_interestpoint.x = x;
				temp_interestpoint.y = y;
				temp_interestpoint.hog_patch_size = patchsize;
				hog::compute3x3HOGAtPoint(mt_src, temp_interestpoint.HOG_descriptors, temp_interestpoint.HOG_desc_locations, cv::Point(x, y), patchsize, numbins);
				ls_locations.push_back(temp_interestpoint);
			}
	}

}

void compute3x3HOGForAllTIPs(const cv::Mat& mt_src, std::vector<ip::interestpoint>& vc_ip_io, const int& numbins)
{
	//If vc_ip_io is empty, for loop exits immediately, no assertion necessary
	for(std::vector<ip::interestpoint>::iterator it_vc_ip = vc_ip_io.begin(); it_vc_ip != vc_ip_io.end();)
	{
		// Patchsize is defined by sigma_l_sq from TIP
		(*it_vc_ip).hog_patch_size = (cv::saturate_cast<int>(18 * sqrt((*it_vc_ip).sigma_l_sq)/3))*3;
		if (hog::compute3x3HOGAtPoint(mt_src, (*it_vc_ip).HOG_descriptors, (*it_vc_ip).HOG_desc_locations, cv::Point((*it_vc_ip).x, (*it_vc_ip).y), (*it_vc_ip).hog_patch_size, numbins) != 0)
		{
			it_vc_ip = vc_ip_io.erase(it_vc_ip);
		}
		else
		{
			it_vc_ip++;
		}
	}
}

void compute3x3HOGForAllBRISKIPs(const cv::Mat& mt_src, std::vector<ip::interestpoint>& vc_ip_io, const int& numbins)
{
	int patchsize;
	//If vc_ip_io is empty, for loop exits immediately, no assertion necessary
	for(std::vector<ip::interestpoint>::iterator it_vc_ip = vc_ip_io.begin(); it_vc_ip != vc_ip_io.end();)
	{
		// Patchsize is defined brisk_size
		patchsize = round((*it_vc_ip).brisk_size * 2 / 3) * 3; //must be divisible through 3 for 3x3 hog
		(*it_vc_ip).hog_patch_size = patchsize;
		if (hog::compute3x3HOGAtPoint(mt_src, (*it_vc_ip).HOG_descriptors, (*it_vc_ip).HOG_desc_locations, cv::Point((*it_vc_ip).x, (*it_vc_ip).y), (*it_vc_ip).hog_patch_size, numbins) != 0)
		{
			it_vc_ip = vc_ip_io.erase(it_vc_ip);
		}
		else
		{
			it_vc_ip++;
		}
	}
}

int compute3x3HOGAtPoint(const cv::Mat& mt_src, std::vector<float>& vc_descriptors, std::vector<cv::Point>& vc_locations, const cv::Point& position, const int& patchsize, const int& numbins)
{
	std::vector<float> zero_descr_vector((numbins * 9), 0.0); //vector with zeros
	int x_upperleft = position.x - cv::saturate_cast<int>(patchsize/2);
	int y_upperleft = position.y - cv::saturate_cast<int>(patchsize/2);
	cv::Rect STIPROI(x_upperleft, y_upperleft, patchsize, patchsize);
	//If whole patch is inside the image..
	if( 0 <= STIPROI.x && 0 <= STIPROI.width && STIPROI.x + STIPROI.width <= mt_src.cols && 0 <= STIPROI.y && 0 <= STIPROI.height && STIPROI.y + STIPROI.height <= mt_src.rows)
	{
		cv::Mat mt_patch = mt_src(STIPROI);
		hog::compute3x3HOG(mt_patch, vc_descriptors, vc_locations, numbins);
	}
	else // If patch is outside the image, return zero HOG
	{
		vc_descriptors = zero_descr_vector;
		return -1;
	}
	return 0;
}

//Calculates HOG on whole mt_src dividing it in 3x3 cells
//mt_src must have same width and height
//mt_src width and height must be whole multiples of 3
void compute3x3HOG(const cv::Mat& mt_src, std::vector<float>& vc_descriptors, std::vector<cv::Point>& vc_locations, const int& numbins)
{
	assert(mt_src.cols == mt_src.rows);
	assert((mt_src.cols % 3 == 0) && (mt_src.rows % 3 == 0));

	cv::HOGDescriptor hog;
    hog.winSize = mt_src.size();
    hog.blockSize = mt_src.size();
    hog.cellSize = cv::Size(mt_src.rows/3, mt_src.cols/3); // = 1/3 (blocksize)
    hog.blockStride = cv::Size(1, 1); // can be whatever >0; num of blocks is 1 if winSize = blockSize
    hog.nbins = numbins;

    cv::Mat mt_src_gray;
    cvtColor( mt_src, mt_src_gray, CV_BGR2GRAY );
	hog.compute(mt_src_gray, vc_descriptors, cv::Size(0,0), cv::Size(0,0), vc_locations);
}

void draw3x3HOGVisualizationForAllIPs(cv::Mat& mt_ioimage, std::vector<ip::interestpoint>& ls_IPs, const int& numbins, const float& scale)
{
	std::vector<cv::Point> vc_locations;
	//if ls_STIPs is empty for loop exits immediately, no assertion necessary
	for(std::vector<ip::interestpoint>::iterator it = ls_IPs.begin(); it != ls_IPs.end();)
	{
		if(!(*it).HOG_descriptors.empty() && (*it).x >= 0 && (*it).y >=0)
		{
			draw3x3HOGVisualizationAtPoint(mt_ioimage, (*it).HOG_descriptors, vc_locations, cv::Point((*it).x, (*it).y), (*it).hog_patch_size, numbins, scale);
		}
		it++;
	}
}

// adapted from http://www.juergenwiki.de/old_wiki/doku.php?id=public:hog_descriptor_computation_and_visualization
void draw3x3HOGVisualizationAtPoint(cv::Mat& mt_ioimage, const std::vector<float>& vc_descriptors, const std::vector<cv::Point>& vc_locations, const cv::Point& position, const int& patchsize, const int& numbins, const float& scale)
{

	assert(patchsize % 3 == 0);
	assert(!vc_descriptors.empty());

	//patch size
	const int DIMX = patchsize;
	const int DIMY = patchsize;
	int cellSize = patchsize / 3; // CHANGED BY MAT
	int gradientBinSize = numbins; // CHANGED BY MAT
	float radRangeForOneBin = (float)(CV_PI/(float)gradientBinSize); // dividing 180 into no of bins, how large (in rad) is one bin?

	//visualization scale, to see the lines better
	float linescale = 2.5;
	int linewidth = 1;

	//zooming on output image
	float zoomFac = scale;
    //patch position
    int x_upperleft = (position.x - cv::saturate_cast<int>(patchsize/2)) * zoomFac;
    int y_upperleft = (position.y - cv::saturate_cast<int>(patchsize/2)) * zoomFac;

	// prepare data structure: 9 orientation / gradient strenghts for each cell
	int cells_in_x_dir = DIMX / cellSize;
	int cells_in_y_dir = DIMY / cellSize;
	float*** gradientStrengths = new float**[cells_in_y_dir];
	int** cellUpdateCounter   = new int*[cells_in_y_dir];
	for (int y=0; y<cells_in_y_dir; y++)
	{
		gradientStrengths[y] = new float*[cells_in_x_dir];
		cellUpdateCounter[y] = new int[cells_in_x_dir];
		for (int x=0; x<cells_in_x_dir; x++)
		{
			gradientStrengths[y][x] = new float[gradientBinSize];
			cellUpdateCounter[y][x] = 0;

			for (int bin=0; bin<gradientBinSize; bin++)
				gradientStrengths[y][x][bin] = 0.0;
		}
	}

	// compute gradient strengths per cell
	int descriptorDataIdx = 0;
	int cellx = 0;
	int celly = 0;

	// 9 cells per block ...
	for (celly=0; celly<cells_in_y_dir; celly++)
	{
		for (cellx=0; cellx<cells_in_x_dir; cellx++)
		{
			for (int bin=0; bin<gradientBinSize; bin++)
			{
				float gradientStrength = vc_descriptors[ descriptorDataIdx ];
				descriptorDataIdx++;
				gradientStrengths[cellx][celly][bin] = gradientStrength;
			} // for (all bins)
		}
	}

	// draw cells
	for (celly=0; celly<cells_in_y_dir; celly++)
	{
		for (cellx=0; cellx<cells_in_x_dir; cellx++)
		{
			int drawX = cellx * cellSize;
			int drawY = celly * cellSize;
			int mx = drawX + cellSize/2;
			int my = drawY + cellSize/2;
			cv::rectangle(mt_ioimage, cv::Point(x_upperleft + (int)(drawX*zoomFac), y_upperleft + (int)(drawY*zoomFac)), cv::Point(x_upperleft + (int)((drawX+cellSize)*zoomFac), y_upperleft + (int)((drawY+cellSize)*zoomFac)), cv::Scalar(100,100,100), linewidth);

			// draw in each cell all gradient strengths
			for (int bin=0; bin<gradientBinSize; bin++)
			{
				float currentGradStrength = gradientStrengths[celly][cellx][bin];
				// no line to draw?
				if (currentGradStrength==0)
					continue;
				float currRad = bin * radRangeForOneBin + radRangeForOneBin/2;
				float dirVecX = cos( currRad );
				float dirVecY = sin( currRad );
				float maxVecLen = (float)(cellSize/2.f);

				// compute line coordinates
				float x1 = mx - dirVecX * currentGradStrength * maxVecLen * linescale;
				float y1 = my - dirVecY * currentGradStrength * maxVecLen * linescale;
				float x2 = mx + dirVecX * currentGradStrength * maxVecLen * linescale;
				float y2 = my + dirVecY * currentGradStrength * maxVecLen * linescale;

				// draw gradient visualization
				cv::line(mt_ioimage, cv::Point(x_upperleft + (int)(x1*zoomFac), y_upperleft + (int)(y1*zoomFac)), cv::Point(x_upperleft + (int)(x2*zoomFac), y_upperleft + (int)(y2*zoomFac)), cv::Scalar(0,255,0), linewidth);
			} // for (all bins)
		} // for (cellx)
	} // for (celly)

	// free memory allocated by helper data structures!
	for (int y=0; y<cells_in_y_dir; y++)
	{
		for (int x=0; x<cells_in_x_dir; x++)
		{
			delete[] gradientStrengths[y][x];
		}
		delete[] gradientStrengths[y];
		delete[] cellUpdateCounter[y];
	}
	delete[] gradientStrengths;
	delete[] cellUpdateCounter;

}

// adapted from http://www.juergenwiki.de/old_wiki/doku.php?id=public:hog_descriptor_computation_and_visualization
cv::Mat get3x3HOGVisualizationAtPoint(const cv::Mat& mt_src, const std::vector<float>& vc_descriptors, const std::vector<cv::Point>& vc_locations, const cv::Point& position, const int& patchsize, const int& numbins)
{
	assert(patchsize % 3 == 0);

	//patch size
	const int DIMX = patchsize, DIMY = patchsize;
	int cellSize        = patchsize / 3; // CHANGED BY MAT
	int gradientBinSize = numbins; // CHANGED BY MAT
	float radRangeForOneBin = (float)(CV_PI/(float)gradientBinSize); // dividing 180 into no of bins, how large (in rad) is one bin?

	//visualization scale, to see the lines better
	float linescale = 2.5;

	//zooming on output image
	float zoomFac = 2;
    cv::Mat visu;
    resize(mt_src, visu, cv::Size(cv::saturate_cast<int>(mt_src.cols*zoomFac), cv::saturate_cast<int>(mt_src.rows*zoomFac)));

    //patch position
    int x_upperleft = (position.x - cv::saturate_cast<int>(patchsize/2)) * zoomFac;
    int y_upperleft = (position.y - cv::saturate_cast<int>(patchsize/2)) * zoomFac;

	// prepare data structure: 9 orientation / gradient strenghts for each cell
	int cells_in_x_dir = DIMX / cellSize, cells_in_y_dir = DIMY / cellSize;
	float*** gradientStrengths = new float**[cells_in_y_dir];
	int** cellUpdateCounter   = new int*[cells_in_y_dir];
	for (int y=0; y<cells_in_y_dir; y++)
	{
		gradientStrengths[y] = new float*[cells_in_x_dir];
		cellUpdateCounter[y] = new int[cells_in_x_dir];
		for (int x=0; x<cells_in_x_dir; x++)
		{
			gradientStrengths[y][x] = new float[gradientBinSize];
			cellUpdateCounter[y][x] = 0;

			for (int bin=0; bin<gradientBinSize; bin++)
				gradientStrengths[y][x][bin] = 0.0;
		}
	}

	// compute gradient strengths per cell
	int descriptorDataIdx = 0;
	int cellx = 0;
	int celly = 0;

	// 9 cells per block ...
	for (celly=0; celly<cells_in_y_dir; celly++)
	{
		for (cellx=0; cellx<cells_in_x_dir; cellx++)
		{
			for (int bin=0; bin<gradientBinSize; bin++)
			{
				float gradientStrength = vc_descriptors[ descriptorDataIdx ];
				descriptorDataIdx++;
				gradientStrengths[cellx][celly][bin] += gradientStrength;
			} // for (all bins)
		}
	}

	// draw cells
	for (celly=0; celly<cells_in_y_dir; celly++)
	{
		for (cellx=0; cellx<cells_in_x_dir; cellx++)
		{
			int drawX = cellx * cellSize, drawY = celly * cellSize;
			int mx = drawX + cellSize/2, my = drawY + cellSize/2;
			cv::rectangle(visu, cv::Point(x_upperleft + (int)(drawX*zoomFac), y_upperleft + (int)(drawY*zoomFac)), cv::Point(x_upperleft + (int)((drawX+cellSize)*zoomFac), y_upperleft + (int)((drawY+cellSize)*zoomFac)), cv::Scalar(100,100,100), 1);

			// draw in each cell all gradient strengths
			for (int bin=0; bin<gradientBinSize; bin++)
			{
				float currentGradStrength = gradientStrengths[celly][cellx][bin];
				// no line to draw?
				if (currentGradStrength==0)
					continue;
				float currRad = bin * radRangeForOneBin + radRangeForOneBin/2;
				float dirVecX = cos( currRad ), dirVecY = sin( currRad );
				float maxVecLen = (float)(cellSize/2.f);

				// compute line coordinates
				float x1 = mx - dirVecX * currentGradStrength * maxVecLen * linescale;
				float y1 = my - dirVecY * currentGradStrength * maxVecLen * linescale;
				float x2 = mx + dirVecX * currentGradStrength * maxVecLen * linescale;
				float y2 = my + dirVecY * currentGradStrength * maxVecLen * linescale;

				// draw gradient visualization
				cv::line(visu, cv::Point(x_upperleft + (int)(x1*zoomFac), y_upperleft + (int)(y1*zoomFac)), cv::Point(x_upperleft + (int)(x2*zoomFac), y_upperleft + (int)(y2*zoomFac)), cv::Scalar(0,255,0), 1);
			} // for (all bins)
		} // for (cellx)
	} // for (celly)

	// free memory allocated by helper data structures!
	for (int y=0; y<cells_in_y_dir; y++)
	{
		for (int x=0; x<cells_in_x_dir; x++)
		{
			delete[] gradientStrengths[y][x];
		}
		delete[] gradientStrengths[y];
		delete[] cellUpdateCounter[y];
	}
	delete[] gradientStrengths;
	delete[] cellUpdateCounter;

	return visu;
}

// adapted from http://www.juergenwiki.de/old_wiki/doku.php?id=public:hog_descriptor_computation_and_visualization
cv::Mat get3x3HOGVisualization(const cv::Mat& mt_src, const std::vector<float>& descriptorValues, const int& numbins)
{
	assert(mt_src.cols == mt_src.rows);
	assert((mt_src.cols % 3 == 0)&&(mt_src.rows % 3 == 0));

	//size is window size
    const int DIMX = mt_src.size().width, DIMY = mt_src.size().height;
    int cellSize        = mt_src.rows/3; // CHANGED BY MAT
    int gradientBinSize = numbins; // CHANGED BY MAT
    float radRangeForOneBin = (float)(CV_PI/(float)gradientBinSize); // dividing 180 into no of bins, how large (in rad) is one bin?

    //visualization scale, to see the lines better
    float linescale = 2.5;

    //return matrix
    float zoomFac = 1;
    cv::Mat visu;
    resize(mt_src, visu, cv::Size( (int)(mt_src.cols*zoomFac), (int)(mt_src.rows*zoomFac) ) );

    // prepare data structure: 9 orientation / gradient strenghts for each cell
    int cells_in_x_dir = DIMX / cellSize, cells_in_y_dir = DIMY / cellSize;
    float*** gradientStrengths = new float**[cells_in_y_dir];
    int** cellUpdateCounter   = new int*[cells_in_y_dir];
    for (int y=0; y<cells_in_y_dir; y++)
    {
        gradientStrengths[y] = new float*[cells_in_x_dir];
        cellUpdateCounter[y] = new int[cells_in_x_dir];
        for (int x=0; x<cells_in_x_dir; x++)
        {
            gradientStrengths[y][x] = new float[gradientBinSize];
            cellUpdateCounter[y][x] = 0;

            for (int bin=0; bin<gradientBinSize; bin++)
                gradientStrengths[y][x][bin] = 0.0;
        }
    }

    // compute gradient strengths per cell
    int descriptorDataIdx = 0, cellx = 0, celly = 0;

    // 9 cells per block ...
    for (celly=0; celly<cells_in_y_dir; celly++)
    {
		for (cellx=0; cellx<cells_in_x_dir; cellx++)
		{
			for (int bin=0; bin<gradientBinSize; bin++)
			{
				float gradientStrength = descriptorValues[ descriptorDataIdx ];
				descriptorDataIdx++;
				gradientStrengths[cellx][celly][bin] += gradientStrength;
			} // for (all bins)
		}
    }

    // draw cells
    for (celly=0; celly<cells_in_y_dir; celly++)
    {
        for (cellx=0; cellx<cells_in_x_dir; cellx++)
        {
            int drawX = cellx * cellSize, drawY = celly * cellSize;
            int mx = drawX + cellSize/2, my = drawY + cellSize/2;
            cv::rectangle(visu, cv::Point((int)(drawX*zoomFac), (int)(drawY*zoomFac)), cv::Point((int)((drawX+cellSize)*zoomFac), (int)((drawY+cellSize)*zoomFac)), cv::Scalar(100,100,100), 1);

            // draw in each cell all gradient strengths
            for (int bin=0; bin<gradientBinSize; bin++)
            {
                float currentGradStrength = gradientStrengths[celly][cellx][bin];
                // no line to draw?
                if (currentGradStrength==0)
                    continue;
                float currRad = bin * radRangeForOneBin + radRangeForOneBin/2;
                float dirVecX = cos( currRad ), dirVecY = sin( currRad );
                float maxVecLen = (float)(cellSize/2.f);

                // compute line coordinates
                float x1 = mx - dirVecX * currentGradStrength * maxVecLen * linescale;
                float y1 = my - dirVecY * currentGradStrength * maxVecLen * linescale;
                float x2 = mx + dirVecX * currentGradStrength * maxVecLen * linescale;
                float y2 = my + dirVecY * currentGradStrength * maxVecLen * linescale;
                // draw gradient visualization
                cv::line(visu, cv::Point((int)(x1*zoomFac),(int)(y1*zoomFac)), cv::Point((int)(x2*zoomFac),(int)(y2*zoomFac)), cv::Scalar(0,255,0), 1);
            } // for (all bins)
        } // for (cellx)
    } // for (celly)

    // free memory allocated by helper data structures!
    for (int y=0; y<cells_in_y_dir; y++)
    {
        for (int x=0; x<cells_in_x_dir; x++)
        {
            delete[] gradientStrengths[y][x];
        }
        delete[] gradientStrengths[y];
        delete[] cellUpdateCounter[y];
    }
    delete[] gradientStrengths;
    delete[] cellUpdateCounter;

    return visu;
} // get_hogdescriptor_visu

hogdescriptor::hogdescriptor() {
	// TODO Auto-generated constructor stub

}

hogdescriptor::~hogdescriptor() {
	// TODO Auto-generated destructor stub
}

} /* namespace hog */
