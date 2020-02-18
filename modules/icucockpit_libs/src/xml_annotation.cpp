/*
 * xml_annotation.cpp
 *
 * Created on: 2017
 *     Author: aed
 */

#include "xml_annotation.h"

#include <iostream>
#include <string>     // std::string, std::stod
#include <ctime>
#include <iomanip>
#include <sstream>

#include "libxml/parser.h" //needs libxml2-dev
#include "libxml/tree.h" //needs libxml2-dev

namespace anxml
{

bool loadXmlFile(const std::string& p_file_name, std::vector<xml_event>& vc_xml_events)
{
	vc_xml_events.clear();

	xmlDocPtr doc; /* Pointer to the parsed document */
    xmlNodePtr cur = NULL; /* Declare a node pointer */
    xmlNodePtr cur_1 = NULL; /* Declare a node pointer */
    xmlNodePtr cur_2 = NULL; /* Declare a node pointer */

    doc = xmlParseFile(p_file_name.c_str());
    /* check if parsing suceeded */
    if (doc == NULL)
    {
        std::cerr<< "Failed to parse " <<p_file_name<<std::endl;
        return false;
    }

    cur = xmlDocGetRootElement(doc);
	if (cur == NULL)
	{
		std::cerr<< "Empty xml file: " <<p_file_name<<std::endl;
		xmlFreeDoc(doc);
		return false;
	}

	//std::cout<< "Root node: " <<cur->name<<std::endl;

	if (xmlStrcmp(cur->name, (const xmlChar *) "Events"))
	{
		std::cerr<<"Document of the wrong type, root node != Events"<<std::endl;
		xmlFreeDoc(doc);
		return false;
	}

	xmlChar *key;
	cur = cur->xmlChildrenNode;
	while (cur != NULL) {
		if ((!xmlStrcmp(cur->name, (const xmlChar *)"Event")))
		{
			//std::cout<< "Children node: " <<cur->name<<std::endl;
			xml_event temp_event;
			cur_1 = cur->xmlChildrenNode;
			while (cur_1 != NULL)
			{
				if ((!xmlStrcmp(cur_1->name, (const xmlChar *)"Location")))
				{
					//std::cout<< "Children node 1: " <<cur_1->name<<std::endl;
					cur_2 = cur_1->xmlChildrenNode;
					while (cur_2 != NULL)
					{
						if ((!xmlStrcmp(cur_2->name, (const xmlChar *)"StartTime")))
						{
							//std::cout<< "Children node 2: " <<cur_2->name<<std::endl;
							key = xmlNodeListGetString(doc, cur_2->xmlChildrenNode, 1);
							if(key != NULL)
							{
								temp_event.start_time = boost::posix_time::time_from_string(reinterpret_cast<const char*>(key));
								temp_event.has_start_time = true;
								xmlFree(key);
							}
							else
							{
								std::cerr<<"ERROR: Xml element \"StartTime\" is empty."<<std::endl;
							}

						}
						if ((!xmlStrcmp(cur_2->name, (const xmlChar *)"Duration")))
						{
							key = xmlNodeListGetString(doc, cur_2->xmlChildrenNode, 1);
							if(key != NULL)
							{
								std::string duration = reinterpret_cast<const char*>(key);
								std::string::size_type sz;     // alias of size_t
								temp_event.duration = boost::posix_time::seconds(std::stod(duration, &sz));
								temp_event.has_duration = true;
								xmlFree(key);
							}
							else
							{
								std::cerr<<"ERROR: Xml element \"Duration\" is empty."<<std::endl;
							}
						}
						if ((!xmlStrcmp(cur_2->name, (const xmlChar *)"EventLabel")))
						{
							key = xmlNodeListGetString(doc, cur_2->xmlChildrenNode, 1);
							if(key != NULL)
							{
								temp_event.event_label = reinterpret_cast<const char*>(key);
								temp_event.has_event_label = true;
								xmlFree(key);
							}
							else
							{
								std::cerr<<"ERROR: Xml element \"EventLabel\" is empty."<<std::endl;
							}
						}
						cur_2 = cur_2->next;
					}
					if (temp_event.has_start_time && temp_event.has_duration && temp_event.has_event_label)
					{
						temp_event.end_time = boost::posix_time::ptime(temp_event.start_time + temp_event.duration);
						temp_event.has_end_time = true;
					}

				}
				cur_1 = cur_1->next;
			}
			//Add to vector
			if (temp_event.has_start_time && temp_event.has_duration && temp_event.has_event_label && temp_event.has_end_time)
			{
				vc_xml_events.push_back(temp_event);
			}
			else
			{
				std::cerr<<"ERROR: There is an incomplete event in the xml annotation file."<<std::endl;
				vc_xml_events.clear();
				xmlFreeDoc(doc);
				return false;
			}
		}
	cur = cur->next;
	}

	//cleanup
	xmlFreeDoc(doc);
    /*
     * Cleanup function for the XML library.
     */
    //xmlCleanupParser();

	std::cout<< "Annotation data loaded from: " <<p_file_name<<std::endl;
	return true;
}

int findEvent(const boost::posix_time::ptime& time_point, const std::vector<xml_event>& vc_xml_events, std::vector<std::string>& vc_labels)
{
	vc_labels.clear();
	int events_found = 0;
	//iterate through vector
	for (std::vector<xml_event>::const_iterator vc_iter = vc_xml_events.begin(); vc_iter != vc_xml_events.end(); ++vc_iter)
		{
			if((*vc_iter).has_start_time && (*vc_iter).has_end_time && (*vc_iter).has_event_label)
			{
				if (((*vc_iter).start_time <= time_point) && (time_point <= (*vc_iter).end_time))
				{
					vc_labels.push_back((*vc_iter).event_label);
					events_found++;
				}
			}
		}
	return events_found;
}

void printAllLabels(const std::vector<std::string>& vc_labels)
{
	for (std::vector<std::string>::const_iterator vc_iter = vc_labels.begin(); vc_iter != vc_labels.end(); ++vc_iter)
		std::cout <<"Label n: "<< (*vc_iter) << std::endl;
}

}
