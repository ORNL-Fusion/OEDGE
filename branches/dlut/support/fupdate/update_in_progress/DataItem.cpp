#include <iostream>
#include <sstream>

using namespace std;

#include "Local_Types.h"
#include "DataItem.h"
#include "DataFile.h"


// C++ for DataItem related classes


// -------------------------------------------------
// DataItem 
// -------------------------------------------------



DataItem::DataItem()
{
  tag = 'NONE';
  tagtype = none;
  desc = '?';
  value_string='?'; 
  datatype = errd;
  isarray = 0;
  copystatus = 0;
};

DataItem::DataItem(dtypes dtype)
{
  tag = 'NONE';
  tagtype = none;
  desc = '?';
  value_string='?'; 
  datatype = dtype;
  isarray = 0;
  copystatus = 0;
};

DataItem::DataItem(string description)
{
  tag = 'NONE';
  tagtype = none;
  desc = description;
  value_string='?'; 
  datatype = errd;
  isarray = 0;
  copystatus = 0;

};


DataItem::DataItem(string description,dtypes dtype)
{
  tag = 'NONE';
  tagtype = none;
  desc = description;
  value_string='?'; 
  datatype = dtype;
  isarray = 0;
  copystatus = 0;
};


DataItem::DataItem(string description,string data,dtypes dtype)
{
  tag = 'NONE';
  tagtype = none;
  desc = description;
  value_string=data; 
  datatype = dtype;
  isarray = 0;
  copystatus = 0;
};


void DataItem::echodata()
{
  cout << tag << ":" << desc << endl;
  // cout << desc << endl;
}; 

int DataItem::actvalue()
{ 
  return (0);
}; 

int DataItem::compare(DataItem* item2)
{ 
  // return (desc.matches(item2->desc));
  return (desc==item2->desc);
}; 

int DataItem::istype(dtypes dtype)
{

  return ((datatype==dtype) ? TRUE : FALSE);

};

void DataItem::setarray(int num)
{
  isarray=num;
};

int DataItem::getarray()
{
  return isarray;
};

void DataItem::settype(dtypes dtype)
{
  datatype=dtype;
};

dtypes DataItem::gettype()
{
  return datatype;
};



DataItem* DataItem::NewDataItem(dtypes newditem)
{
  DataItem* tmpdata;
  string tmp_string1,tmp_string2;

  tmp_string1='?';
  tmp_string2='?';

  tmpdata=NewDataItem(newditem,tmp_string1,tmp_string2);
   
  return tmpdata;

};


DataItem* DataItem::NewDataItem(dtypes newditem, string d)
{
  DataItem* tmpdata;
  string tmp_string1;

  tmp_string1='?';

  tmpdata=NewDataItem(newditem,d,tmp_string1);
   
  return tmpdata;

};



DataItem* DataItem::NewDataItem(dtypes newditem, string d, string z)
{
  DataItem* tmpdata;

  string comp_string;

  comp_string = '?';

  switch(newditem)           
    {
              
    case errd:      
    case comn:
    case comd:   
      
      if (d==comp_string)     {  tmpdata = new DataItem(newditem)    ;}
      else if (z==comp_string){  tmpdata = new DataItem(d,newditem)  ;}
      else            {  tmpdata = new DataItem(d,z,newditem);} ;

      break;
      
    case strg:     
    case stra:
    case strd:
    case strp:

      if (d==comp_string)     {  tmpdata = new string_DataItem(newditem)    ;} 
      else if (z==comp_string){  tmpdata = new string_DataItem(d,newditem)  ;}
      else            {  tmpdata = new string_DataItem(d,z,newditem);};

      break;
                              
    case intd:       /* 1 integer */
			      
      if (d==comp_string)     {  tmpdata = new Int_DataItem(newditem)    ;}
      else if (z==comp_string){  tmpdata = new Int_DataItem(d,newditem)  ;}
      else            {  tmpdata = new Int_DataItem(d,z,newditem);};

      break;
      
    case reald:       /* 1 real*/
      
      if (d==comp_string)     {  tmpdata = new Real_DataItem(newditem)    ;} 
      else if (z==comp_string){  tmpdata = new Real_DataItem(d,newditem)  ;}
      else            {  tmpdata = new Real_DataItem(d,z,newditem);};

      break;
      
    case real2d:       /* 2 real */
      
      if (d==comp_string)     {  tmpdata = new Real2_DataItem(newditem)    ;}
      else if (z==comp_string){  tmpdata = new Real2_DataItem(d,newditem)  ;}
      else            {  tmpdata = new Real2_DataItem(d,z,newditem);};

      break;
      
    case real3d:       /* 3 real*/
      
      if (d==comp_string)     {  tmpdata = new Real3_DataItem(newditem)    ;}
      else if (z==comp_string){  tmpdata = new Real3_DataItem(d,newditem)  ;}
      else            {  tmpdata = new Real3_DataItem(d,z,newditem);};

      break;
      
    default:
      
      if (d==comp_string)     {   tmpdata = new DataItem(newditem)    ;}
      else if (z==comp_string){   tmpdata = new DataItem(d,newditem)  ;}
      else            {   tmpdata = new DataItem(d,z,newditem);};

      break; 
      
    };
  
  return tmpdata; 

};



DataItem* DataItem::NewDataItem(DataItem* original)
{
  DataItem* tmpdata; 
  dtypes tmpdatatype;

  tmpdata = DataItem::NewDataItem();

  *tmpdata=*original;

  return tmpdata;
 
};


  
DataItem* DataItem::analyse_line(string data)
{
        
  DataItem* tmpdata;
  string d,z;
  int num;
  dtypes ind; 

  // cout << "ANALYSE:" << data << endl;

  if (data.find("'",0)!=string::npos)

    /* Data line - all else are comments or data*/

    {

      num = data.find("'",1);

      d = data.substr(0,num+1);    
      z = data.substr(num+1,data.size());            

      ind = idtype(z);

      // cout << "Creating Input Line:" << num << ":" << ind << ":" << d << endl;

      tmpdata = DataItem::NewDataItem(ind,d,z);


    } 
  else if (data.find("$",0)!=string::npos)
    {
      /* Comment line*/
               
      // cout << "Creating comment:" << data.find("$",0) << ":" << data << endl;

      tmpdata = DataItem::NewDataItem(comd,data);

    }   
  else if ((data.find("c",0)!=string::npos)||(data.find("C",0)!=string::npos))
    {
      /* NIMBUS Comment line - invalid DIVIMP input */
               
      // cout << "Creating NIMBUS comment:" << data.find("c",0) << ":" << data << endl;

      tmpdata = DataItem::NewDataItem(comn,data);

    }   
  else 
    {
      /* Line containing just data OR an error condition */

      // cout << "Creating data line:" << data.find("$",0) << ":" << data << endl;

      tmpdata = DataItem::NewDataItem(errd,data);
      
    };

  // Check if the input line contains a tag - store the tag type and identifier if present 
  // Tag types are * for optional input which can appear in any order and + for standard input that must appear in the same order 
  // Comparisons will look for tags first if present and then compare the subsequent 20 characters of the description string when 
  // looking for matches

  tmpdata->extracttag();

  return tmpdata;  
          
};


void DataItem::extracttag()
{
  
  // Look at the description string and extract the tag if one is present

  // Tags are defined by the first character being a ' - following by a * for optional or a + for required - the next 3 characters are the tag

  // Default is no tag

  tagtype = none;
  // Use four characters so that this can never match a valid tag
  tag="NONE";

  if (desc.substr(0,1)=="'") 
    {
      if (desc.substr(1,1)=="*") 
	{
	  tagtype = opt;
	  tag = desc.substr(2,3);

	}
      else if (desc.substr(1,1)=="+") 
	{
	  tagtype=req;
	  tag = desc.substr(2,3);
      
	};

    };

};


tagtypes DataItem::gettagtype()
{

  return(tagtype);

};

string DataItem::gettag()
{

  return(tag);

};




dtypes DataItem::idtype(string data)
{

  dtypes retval; 
  int a,count;

  /* Determines the type of data being passed to it
     0 = error or blank string
     1 = character string
     2 = integer
     3 = 1 real
     4 = 2 reals
     5 = 3 reals  */


  if (data.find("'",0)!=string::npos)
    {
      /* Contains charater string */

      retval = strg;

    } 
  else if (data.find(".",0)!=string::npos)
    {
      /* Real numbers ... */

      count = 0;
      a=data.find(".",0);

      // Count number of decimal characters to determine quantity of real numbers on the line

      while (a!=string::npos) 
	{
	count = count + 1;
	a = a+1;
        a=data.find(".",a);
	}; 


      switch (count)
	{

	case 0:

           retval = intd;
           break;

	case 1:

           retval = reald;
           break;

        case 2:

           retval = real2d;
           break;

        case 3:

           retval = real3d;
           break;

        default:

           retval = intd;
           break;
	   
	};

    }
  else if (data.find_first_of("0123456789")!=string::npos) 
    {
      /* contains numeric values */
      retval = intd;
    }
  else
    {
      retval = errd;
    };           


  return (retval);
 

};




int DataItem::operator == (DataItem item2)
{

  // If the tags are the same then the lines are intended to represent the same
  // input data so test the tags first - if the tags are different then they are not the same line
  
  if (tagtype!=none&&item2.tagtype!=none)
    { 
      if (tag==item2.tag) 
	{
          return TRUE;
	}
      else
	return FALSE;
    }
  else if (desc.find("' '",0)!=string::npos&&item2.desc.find("' '",0)!=string::npos) 
    { 
      // If the description is a short empty string then compare the 
      // value_string for equality testing.
      // Strip any trailing spaces and compare the resulting substrings

      int len1 = value_string.length() -1 ;
      int len2 = item2.value_string.length() -1;

      while (len1>=0&&value_string[len1]==' ') len1--;
 
      while (len2>=0&&item2.value_string[len2]==' ') len2--;

      if (len1!=len2) return FALSE;        

      string desc1 = value_string.substr(0,len1+1);
      string desc2 = item2.value_string.substr(0,len2+1); 

      return ((desc1==desc2) ? TRUE : FALSE );

      /* return ((value_string.matches(item2.value_string)) ? TRUE : FALSE );*/
    }
  else
    {
      // Otherwise - compare the description 
      // Strip any trailing spaces as well as any tags at the beginning and compare the resulting descriptions

      int len1 = desc.length() -1 ;
      int len2 = item2.desc.length() -1;

      int start1 = 0;
      int start2 = 0;

      // Tags take up the first 5 characters plus a space typically - so the description will start in string position 6
      if (tagtype!=none) start1 = 6;
      if (item2.tagtype!=none) start2 = 6;

      while (len1>=0&&desc[len1]==' ') len1--;
 
      while (len2>=0&&item2.desc[len2]==' ') len2--;

      if (len1!=len2) return FALSE;        

      string desc1 = desc.substr(start1,len1-start1+1);
      string desc2 = item2.desc.substr(start2,len2-start2+1); 

      return ((desc1==desc2) ? TRUE : FALSE );
    };

}; 




// -------------------------------------------------
// string_DataItem -
// -------------------------------------------------



string_DataItem::string_DataItem() : DataItem()
{
};

string_DataItem::string_DataItem(dtypes dtype) : DataItem(dtype)
{
};

string_DataItem::string_DataItem(string description,dtypes dtype) : DataItem(description,dtype)
{
};

string_DataItem::string_DataItem(string description,string data, dtypes dtype) : DataItem(description,data,dtype)
{
};


void string_DataItem::echodata()
{
  cout << tag << ":" << desc << value_string << endl;

}; 



// -------------------------------------------------
// Int_DataItem - 
// -------------------------------------------------



Int_DataItem::Int_DataItem() : DataItem()
{
};

Int_DataItem::Int_DataItem(dtypes dtype) : DataItem(dtype)
{
};

Int_DataItem::Int_DataItem(string description,dtypes dtype) : DataItem(description,dtype)
{
};

Int_DataItem::Int_DataItem(string description,string data, dtypes dtype) : DataItem(description,data,dtype)
{
  value=atoi(data.c_str());
};


void Int_DataItem::echodata()
{
  cout <<  tag << ":" << desc << value_string << endl;
}; 

int Int_DataItem::actvalue()
{
   return(value);
};


// -------------------------------------------------
// Real_DataItem - 
// -------------------------------------------------



Real_DataItem::Real_DataItem() : DataItem()
{
};

Real_DataItem::Real_DataItem(dtypes dtype) : DataItem(dtype)
{
};

Real_DataItem::Real_DataItem(string description,dtypes dtype) : DataItem(description,dtype)
{
};

Real_DataItem::Real_DataItem(string description,string data, dtypes dtype) : DataItem(description,data,dtype)
{
  value=atof(data.c_str());
};


void Real_DataItem::echodata()
{
  cout <<  tag << ":" << desc << value_string << endl;
}; 



// -------------------------------------------------
// Real2_DataItem - 
// -------------------------------------------------



Real2_DataItem::Real2_DataItem() : DataItem()
{
};

Real2_DataItem::Real2_DataItem(dtypes dtype) : DataItem(dtype)
{
};

Real2_DataItem::Real2_DataItem(string description,dtypes dtype) : DataItem(description,dtype)
{
};

Real2_DataItem::Real2_DataItem(string description, string data,dtypes dtype) : DataItem(description,data,dtype)
{
  istringstream in_string(data);
  in_string >> value1 >> value2;

};


void Real2_DataItem::echodata()
{
  cout <<  tag << ":" << desc << value_string << endl;
}; 



// -------------------------------------------------
// Real3_DataItem - 
// -------------------------------------------------



Real3_DataItem::Real3_DataItem() : DataItem()
{
};

Real3_DataItem::Real3_DataItem(dtypes dtype) : DataItem(dtype)
{
};

Real3_DataItem::Real3_DataItem(string description,dtypes dtype) : DataItem(description,dtype)
{
};


Real3_DataItem::Real3_DataItem(string description,string data,dtypes dtype) : DataItem(description,data,dtype)
{
  istringstream in_string(data);
  in_string >> value1 >> value2 >> value3;
};


void Real3_DataItem::echodata()
{
  cout <<  tag << ":" << desc << value_string << endl;
}; 







