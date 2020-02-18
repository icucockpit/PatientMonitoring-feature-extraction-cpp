/*
 * xml_annotation.h
 *
 * Created on: 2017
 *     Author: aed
 */

#ifndef XML_ANNOTATION_H_
#define XML_ANNOTATION_H_

#include <string>
#include <vector>

//OpenCV:
//#include <opencv2/core/core.hpp>

//#include "boost/date_time/posix_time/posix_time_types.hpp" //no i/o just types
#include "boost/date_time/posix_time/posix_time.hpp" //include all types plus i/o

namespace anxml
{

struct xml_event {
	boost::posix_time::ptime start_time, end_time;
	boost::posix_time::time_duration duration;
	std::string event_label;
	bool has_start_time = false, has_duration = false, has_end_time = false, has_event_label = false;

};

bool loadXmlFile(const std::string& p_file_name, std::vector<xml_event>& vc_xml_events);
int findEvent(const boost::posix_time::ptime& time_point, const std::vector<xml_event>& vc_xml_events, std::vector<std::string>& vc_labels);
void printAllLabels(const std::vector<std::string>& vc_labels);

}
#endif /* XML_ANNOTATION_H_ */

