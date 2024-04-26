#include "ConfigReader.h"

// C++ Header File(s)
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;
using namespace cppsecrets;

#define NULL_PTR 0

ConfigReader* ConfigReader::m_pInstance = NULL_PTR;

ConfigReader::ConfigReader()
{
   m_ConfigSettingMap.clear();
}

ConfigReader::~ConfigReader()
{
   m_ConfigSettingMap.clear();
}

ConfigReader* ConfigReader::getInstance()
{
   // No need to use double re-check lock mechanism here
   // because this getInstance() will call at the time of
   // initialization only and mostly, at the time of
   // initialization, there will be only one thread.

   if(NULL_PTR == m_pInstance)
   {
      m_pInstance = new ConfigReader;
   }
   return m_pInstance;
}

bool ConfigReader::getValue(std::string tag, unsigned int& value)
{
   map<string, string>::iterator it ;
   it = m_ConfigSettingMap.find(tag);
   if(it != m_ConfigSettingMap.end())
   {
      value = static_cast<unsigned int>(abs(atoi((it->second).c_str())));
      return true;
   }
   return false;
}

bool ConfigReader::getValue(std::string tag, int& value)
{
   map<string, string>::iterator it ;
   it = m_ConfigSettingMap.find(tag);
   if(it != m_ConfigSettingMap.end())
   {
      value = atoi((it->second).c_str());
      return true;
   }
   return false;
}

bool ConfigReader::getValue(std::string tag, double& value)
{
   map<string, string>::iterator it ;
   it = m_ConfigSettingMap.find(tag);
   if(it != m_ConfigSettingMap.end())
   {
      value = stod((it->second).c_str());
      return true;
   }
   return false;
}

bool ConfigReader::getValue(std::string tag, std::string& value)
{
   map<string, string>::iterator it ;
   it = m_ConfigSettingMap.find(tag);
   if(it != m_ConfigSettingMap.end())
   {
      value = it->second;
      return true;
   }
   return false;
}

bool ConfigReader::getValue(std::string tag, std::vector<double>& value)
{
   map<string, string>::iterator it;
   it = m_ConfigSettingMap.find(tag);
   if(it != m_ConfigSettingMap.end())
   {
      std::string line = it->second;
      // Find any commas and replace them with whitespace
      size_t comma = 0;
      while (true) {
        /* Locate the substring to replace. */
        comma = line.find(",", comma);
        if (comma == std::string::npos) break;

        /* Make the replacement. */
        line.replace(comma, 1, " ");

        /* Advance comma forward so the next iteration doesn't pick it up as well. */
        comma += 1;
      }

      double entry;
      std::stringstream lineValue(line);
      while(lineValue >> entry){
        value.push_back(entry);
      }
      return true;
   }
   return false;
}

bool ConfigReader::getValue(std::string tag, std::vector<std::string>& value)
{
   map<string, string>::iterator it;
   it = m_ConfigSettingMap.find(tag);
   if(it != m_ConfigSettingMap.end())
   {
      std::string line = it->second;
      // Find any commas and replace them with whitespace
      size_t comma = 0;
      while (true) {
        /* Locate the substring to replace. */
        comma = line.find(",", comma);
        if (comma == std::string::npos) break;

        /* Make the replacement. */
        line.replace(comma, 1, " ");

        /* Advance comma forward so the next iteration doesn't pick it up as well. */
        comma += 1;
      }

      std::string entry;
      std::stringstream lineValue(line);
      while(lineValue >> entry){
        value.push_back(entry);
      }
      return true;
   }
   return false;
}

bool ConfigReader::parseFile(string fileName)
{
   ifstream inputFile;
   inputFile.open(fileName.c_str());
   string delimeter = "=";
   int initPos = 0;

   if (inputFile.fail())
   {
      cout << "Unable to find defaultConfig file" << endl;
      return false;
   }

   string line;
   while(getline(inputFile, line))
   {
      // Remove comment Lines
      size_t found = line.find_first_of('#');
      string configData = line.substr(0, found);

      // Remove ^M from configData
      configData.erase(std::remove(configData.begin(), configData.end(), '\r'), configData.end());

      if(configData.empty())
         continue;

      const size_t length = configData.find(delimeter);

      string tag, value;

      if (length!=string::npos)
      {
         tag   = configData.substr(initPos, length);
         value = configData.substr(length+1);
      }

      // Trim white spaces
      tag   = reduce(tag);
      value = reduce(value);
      
      if(tag.empty() || value.empty())
         continue;

      // Check if any of the tags is repeated more than one times
      // it needs to pick the latest one instead of the old one.

      // Search, if the tag is already present or not
      // If it is already present, then delete an existing one

      std::map<std::string, std::string>::iterator itr = m_ConfigSettingMap.find(tag);
      if(itr != m_ConfigSettingMap.end())
      {
         m_ConfigSettingMap.erase(tag);
      }

      m_ConfigSettingMap.insert(std::pair<string, string>(tag, value));
   }
   return true;
}

std::string ConfigReader::trim(const std::string& str, const std::string& whitespace)
{
   size_t strBegin = str.find_first_not_of(whitespace);
   if (strBegin == std::string::npos)
      return "";

   size_t strEnd = str.find_last_not_of(whitespace);
   size_t strRange = strEnd - strBegin + 1;

   return str.substr(strBegin, strRange);
}

std::string ConfigReader::reduce(const std::string& str,
      const std::string& fill,
      const std::string& whitespace)
{
   // trim first
   string result = trim(str, whitespace);

   // replace sub ranges
   size_t beginSpace = result.find_first_of(whitespace);
   while (beginSpace != std::string::npos)
   {
      size_t endSpace = result.find_first_not_of(whitespace, beginSpace);
      size_t range = endSpace - beginSpace;

      result.replace(beginSpace, range, fill);

      size_t newStart = beginSpace + fill.length();
      beginSpace = result.find_first_of(whitespace, newStart);
   }

   return result;
}

void ConfigReader::dumpFileValues()
{
   map<string, string>::iterator it;
   for(it=m_ConfigSettingMap.begin(); it!=m_ConfigSettingMap.end(); ++it)
   {
      cout << it->first << " = " << it->second << endl;
   }
}