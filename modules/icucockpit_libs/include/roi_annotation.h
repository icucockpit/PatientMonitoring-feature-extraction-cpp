/*
 * Annotation.h
 *
 *  Created on: 2017
 *      Author: aed
 */

#ifndef ROIANNOTATION_H_
#define ROIANNOTATION_H_

#include <fstream>
#include <iostream>
#include <string>
#include <cmath> // for floor etc.

//OpenCV:
#include <opencv2/core/core.hpp>


namespace an
{


class ROIAnnotation {
public:
	/** PUBLIC data members*/
	std::vector<std::vector<std::string> > vc_Object_labels;//Each frame may have multiple objects
	std::vector<std::vector<cv::Rect> > vc_Objects;//A rectangle for each object
	std::vector<std::vector<std::vector<std::string> > > vc_Attributes;//Each object may have multiple attributes

	/** Initialize*/
	ROIAnnotation(const int & p_frames);
	virtual ~ROIAnnotation();

    /** Functions */
	void resizeData();

	int loadAllROIsFromYml(const std::string& p_file_name);
    int saveAllROIsToYml(const std::string& p_file_name);
    int saveOneROIToTxt(const std::string& p_file_name, const int& objectNumber);
    int saveAllROIsToTxt(const std::string& p_file_name);

protected:

    const std::string currentDateTime();

private:
    int frames;

};



}
#endif /* ROIANNOTATION_H_ */
