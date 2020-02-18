/*
 * ROIAnnotation.cpp
 *
 *  Created on: 2017
 *      Author: aed
 */

#include "roi_annotation.h"

namespace an
{

/** Initialize*/

ROIAnnotation::ROIAnnotation(const int & p_frames)
{
    frames = p_frames;
}

ROIAnnotation::~ROIAnnotation()
{
	// TODO Auto-generated destructor stub
}

void ROIAnnotation::resizeData()
{
	vc_Object_labels.resize(frames);//Each frame may have multiple = (5) objects
	for (std::vector<std::vector<std::string> >::iterator vec_vec_str_it = vc_Object_labels.begin(); vec_vec_str_it != vc_Object_labels.end(); ++vec_vec_str_it)
	{
		vec_vec_str_it->resize(5);
	}

	vc_Objects.resize(frames);//A rectangle for each of the (5) objects
	for (std::vector<std::vector<cv::Rect> >::iterator vec_vec_rec_it = vc_Objects.begin(); vec_vec_rec_it != vc_Objects.end(); ++vec_vec_rec_it)
	{
		vec_vec_rec_it->resize(5);
	}

	vc_Attributes.resize(frames);//Each of the (5) objects may have (6) attributes = class labels
	for (std::vector<std::vector<std::vector<std::string> > >::iterator vec_vec_vec_str_it = vc_Attributes.begin(); vec_vec_vec_str_it != vc_Attributes.end(); ++vec_vec_vec_str_it)
	{
		vec_vec_vec_str_it->resize(5);
		for (std::vector<std::vector<std::string> >::iterator vec_vec_str_it = vec_vec_vec_str_it->begin(); vec_vec_str_it != vec_vec_vec_str_it->end(); ++vec_vec_str_it)
		{
			vec_vec_str_it->resize(6);
		}
	}
std::cout<<"Data resized.."<<std::endl;
}



int ROIAnnotation::loadAllROIsFromYml(const std::string& p_file_name)
{
	cv::FileStorage fsANN(p_file_name, cv::FileStorage::READ);
	if (!fsANN.isOpened())
	{
		std::cout<<"An annotation file with the name you specified does not exist."<<std::endl<<"Do you want to create a new one (y/n) ? "<<std::endl;
		char decision;
		std::cin >> decision;
		switch(decision)
		{
			case 'y':
				this->resizeData();
				this->saveAllROIsToYml(p_file_name);
				std::cout<<"New empty annotation file saved"<<std::endl;
				break;
			case 'n':
				return 1;
				break;
		}
	}
	else
	{
		int total_frames = (int)fsANN["total_frames"];
		if (total_frames != frames)
		{
			std::cout<<"The total amount of frames reported in the annotation file does not match with the total amount of frames in the video file.Segmentation faults may occur...Are you sure you want to continue (y/n)?"<<std::endl;
			char decision;
			std::cin >> decision;
			switch(decision)
			{
				case 'y':
					frames = total_frames;
					break;
				case 'n':
					return 1;
					break;
			}
		}

		//Read: std::vector<std::vector<std::string> > vc_Object_labels
		cv::FileNode node = fsANN["Object_labels"];
		if (node.type() != cv::FileNode::SEQ)
		{
			std::cerr << "Object_labels is not a sequence! FAIL" << std::endl;
			return 1;
		}
		if (!vc_Object_labels.empty())
		{
			vc_Object_labels.erase(vc_Object_labels.begin(), vc_Object_labels.end());
		}

		std::vector<std::string> vc_Object_labels_in_each_frame;
		// Go through the node
		for (cv::FileNodeIterator node_node_it = node.begin(); node_node_it != node.end(); ++node_node_it)
		{
			cv::FileNode node_node(*node_node_it);
			if (node_node.type() != cv::FileNode::SEQ)
			{
				std::cerr << "Object_labels in each frame is not a sequence! FAIL" << std::endl;
				return 1;
			}
			if (!vc_Object_labels_in_each_frame.empty())
			{
				vc_Object_labels_in_each_frame.erase(vc_Object_labels_in_each_frame.begin(), vc_Object_labels_in_each_frame.end());
			}
			for (cv::FileNodeIterator node_it = node_node.begin(); node_it != node_node.end(); ++node_it)
			{
				vc_Object_labels_in_each_frame.push_back(*node_it);
			}
			vc_Object_labels.push_back(vc_Object_labels_in_each_frame);
		}

		//Read std::vector<std::vector<cv::Rect> > vc_Objects
		node = fsANN["Objects"];
		if (node.type() != cv::FileNode::SEQ)
		{
			std::cerr << "Objects is not a sequence! FAIL" << std::endl;
			return 1;
		}
		if (!vc_Objects.empty())
		{
			vc_Objects.erase(vc_Objects.begin(), vc_Objects.end());
		}

		std::vector<cv::Rect> vc_Objects_in_each_frame;
		// Go through the node
		for (cv::FileNodeIterator node_node_it = node.begin(); node_node_it != node.end(); ++node_node_it)
		{
			cv::FileNode node_node(*node_node_it);
			if (node_node.type() != cv::FileNode::SEQ)
			{
				std::cerr << "Objects in each frame is not a sequence! FAIL" << std::endl;
				return 1;
			}
			if (!vc_Objects_in_each_frame.empty())
			{
				vc_Objects_in_each_frame.erase(vc_Objects_in_each_frame.begin(), vc_Objects_in_each_frame.end());
			}
			for (cv::FileNodeIterator node_it = node_node.begin(); node_it != node_node.end(); ++node_it)
			{
				vc_Objects_in_each_frame.push_back(cv::Rect((int)(*node_it)["x"],(int)(*node_it)["y"],(int)(*node_it)["w"],(int)(*node_it)["h"]));
			}
			vc_Objects.push_back(vc_Objects_in_each_frame);
		}

		//Read std::vector<std::vector<std::vector<std::string> > > vc_Attributes
		node = fsANN["Attributes"];
		if (node.type() != cv::FileNode::SEQ)
		{
			std::cerr << "Attributes is not a sequence! FAIL" << std::endl;
			return 1;
		}
		if (!vc_Attributes.empty())
		{
			vc_Attributes.erase(vc_Attributes.begin(), vc_Attributes.end());
		}

		std::vector<std::string> vc_Attributes_for_one_object;
		std::vector<std::vector<std::string> > vc_Attributes_in_each_frame;
		// Go through the node
		for (cv::FileNodeIterator node_node_it = node.begin(); node_node_it != node.end(); ++node_node_it)
		{
			cv::FileNode node_node(*node_node_it);
			if (node_node.type() != cv::FileNode::SEQ)
			{
				std::cerr << "Attributes in each frame is not a sequence! FAIL" << std::endl;
				return 1;
			}
			if (!vc_Attributes_in_each_frame.empty())
			{
				vc_Attributes_in_each_frame.erase(vc_Attributes_in_each_frame.begin(), vc_Attributes_in_each_frame.end());
			}
			for (cv::FileNodeIterator node_it = node_node.begin(); node_it != node_node.end(); ++node_it)
			{
				cv::FileNode node_node_node(*node_it);
				if (node_node.type() != cv::FileNode::SEQ)
				{
					std::cerr << "Attributes for one object is not a sequence! FAIL" << std::endl;
					return 1;
				}
				if (!vc_Attributes_for_one_object.empty())
				{
					vc_Attributes_for_one_object.erase(vc_Attributes_for_one_object.begin(), vc_Attributes_for_one_object.end());
				}
				for(cv::FileNodeIterator node_it_x = node_node_node.begin(); node_it_x != node_node_node.end(); ++node_it_x)
				{
					vc_Attributes_for_one_object.push_back(*node_it_x);
				}
				vc_Attributes_in_each_frame.push_back(vc_Attributes_for_one_object);
			}
			vc_Attributes.push_back(vc_Attributes_in_each_frame);
		}
		fsANN.release();
	}

	return 0;
}

int ROIAnnotation::saveAllROIsToYml(const std::string& p_file_name)
{
	//Open file_name
	cv::FileStorage fsANN(p_file_name, cv::FileStorage::WRITE);

	//Write total number of frames
	fsANN << "total_frames" << frames;

	//Write std::vector<std::vector<std::string> > vc_Object_labels
	fsANN << "Object_labels" << "[";
	for (std::vector<std::vector<std::string> >::iterator vec_vec_str_it = vc_Object_labels.begin(); vec_vec_str_it != vc_Object_labels.end(); ++vec_vec_str_it)
	{
		fsANN << "[:";
		for (std::vector<std::string>::iterator vec_str_it = vec_vec_str_it->begin(); vec_str_it != vec_vec_str_it->end(); ++vec_str_it)
		{
		fsANN << *vec_str_it;
		}
		fsANN << "]";
	}
	fsANN << "]";

	//Write std::vector<std::vector<cv::Rect> > vc_Objects
	fsANN << "Objects" << "[";
	for (std::vector<std::vector<cv::Rect> >::iterator vec_vec_rec_it = vc_Objects.begin(); vec_vec_rec_it != vc_Objects.end(); ++vec_vec_rec_it)
	{
		fsANN << "[:";
		for (std::vector<cv::Rect>::iterator vec_rec_it = vec_vec_rec_it->begin(); vec_rec_it != vec_vec_rec_it->end(); ++vec_rec_it)
		{
		fsANN << "{:" << "x" << vec_rec_it->x << "y" << vec_rec_it->y << "w" << vec_rec_it->width << "h" << vec_rec_it->height << "}";
		}
		fsANN << "]";
	}
	fsANN << "]";

	//Write std::vector<std::vector<std::vector<std::string> > > vc_Attributes
	fsANN << "Attributes" << "[";
	for (std::vector<std::vector<std::vector<std::string> > >::iterator vec_vec_vec_str_it = vc_Attributes.begin(); vec_vec_vec_str_it != vc_Attributes.end(); ++vec_vec_vec_str_it)
	{
		fsANN << "[:";
		for (std::vector<std::vector<std::string> >::iterator vec_vec_str_it = vec_vec_vec_str_it->begin(); vec_vec_str_it != vec_vec_vec_str_it->end(); ++vec_vec_str_it)
		{
			fsANN << "[:";
			for (std::vector<std::string>::iterator vec_str_it = vec_vec_str_it->begin(); vec_str_it != vec_vec_str_it->end(); ++vec_str_it)
			{
			fsANN << *vec_str_it;
			}
			fsANN << "]";
		}
		fsANN << "]";
	}
	fsANN << "]";

	fsANN.release();
	return 1;
}


int ROIAnnotation::saveOneROIToTxt(const std::string& p_file_name, const int& objectNumber)
{
	std::ofstream myfile;

			myfile.open ((char*)p_file_name.c_str(), std::ios::trunc); //Any current content is discarded
			if (!myfile.is_open())
			  {
				std::cerr<<"ERROR: Unable to open output file!"<<std::endl;
				exit(EXIT_FAILURE);
			  }
		//Write data header
			myfile << "% Object ID:" << objectNumber <<std::endl<< "% Title: ROI file for one object" <<std::endl<< "% Owner: aed" <<std::endl << "% Date and time saved: " << currentDateTime()  <<std::endl << "% Total_frames:" << frames <<std::endl<< "% Columns: FrameNo x y w h" <<std::endl<< "% x and y are for top left corner" <<std::endl<< "@@DATA@@" <<std::endl;

			//Write std::vector<std::vector<cv::Rect> > vc_Objects
			for(int currentFrame = 0 ; currentFrame < frames; currentFrame++)
			{
				myfile <<currentFrame<< "," << vc_Objects[currentFrame][objectNumber-1].x << "," << vc_Objects[currentFrame][objectNumber-1].y << "," << vc_Objects[currentFrame][objectNumber-1].width << "," << vc_Objects[currentFrame][objectNumber-1].height << std::endl;
			}
			myfile.close();
			return 1;
}


int ROIAnnotation::saveAllROIsToTxt(const std::string& p_file_name)
{
	std::ofstream myfile;

			myfile.open ((char*)p_file_name.c_str(), std::ios::trunc); //Any current content is discarded
			if (!myfile.is_open())
			  {
				std::cerr<<"ERROR: Unable to open output file!"<<std::endl;
				exit(EXIT_FAILURE);
			  }
		//Write data header
			myfile << "% Title: All ROIs for 5 objects" <<std::endl << "% Owner: aed" <<std::endl << "% Date and time saved: " << currentDateTime() <<std::endl << "% Total_frames:" << frames <<std::endl << "% Columns: FrameNo,x(obj1),y(obj1),w(obj1),h(obj1),x(obj2),y(obj2),w(obj2),h(obj2),x(obj3),y(obj3),w(obj3),h(obj3),x(obj4),y(obj4),w(obj4),h(obj4),x(obj5),y(obj5),w(obj5),h(obj5)" <<std::endl << "% x and y are for top left corner" <<std::endl << "% Convention:" <<std::endl << "% obj1: ?" <<std::endl << "% obj2: ?" <<std::endl << "% obj3: ?" <<std::endl << "% obj4: ?" <<std::endl << "% obj5: ?" <<std::endl << "@@DATA@@" <<std::endl;

			//Write std::vector<std::vector<cv::Rect> > vc_Objects
			for(int currentFrame = 0 ; currentFrame < frames; currentFrame++)
			{
				myfile << currentFrame;
				for(int objectNumber = 1 ; objectNumber <= 5 ; objectNumber++)
				{
					myfile << "," << vc_Objects[currentFrame][objectNumber-1].x << "," << vc_Objects[currentFrame][objectNumber-1].y << "," << vc_Objects[currentFrame][objectNumber-1].width << "," << vc_Objects[currentFrame][objectNumber-1].height;
				}
				myfile << std::endl;
			}
			myfile.close();
			return 1;
}


// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
// from  http://stackoverflow.com/questions/997946/c-get-current-time-and-date
const std::string ROIAnnotation::currentDateTime() {
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
    return buf;
}
}
