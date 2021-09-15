#include <iostream.h>
#include <strstream.h>
#include <fstream.h> 
#include <String.h>
#include <stdlib.h>
#include <math.h>


#include "Local_Types.h"
#include "DataItem.h"
#include "DataFile.h"


// C++ for DataItem related classes


// -------------------------------------------------
// DataItem 
// -------------------------------------------------



DataItem::DataItem()
{
  desc = '?';
  value_string='?'; 
  datatype = errd;
  isarray = 0;
};

DataItem::DataItem(dtypes dtype)
{
  desc = '?';
  value_string='?'; 
  datatype = dtype;
  isarray = 0;
};

DataItem::DataItem(String description)
{
  desc = description;
  value_string='?'; 
  datatype = errd;
  isarray = 0;

};


DataItem::DataItem(String description,dtypes dtype)
{
  desc = description;
  value_string='?'; 
  datatype = dtype;
  isarray = 0;
};


DataItem::DataItem(String description,String data,dtypes dtype)
{
  desc = description;
  value_string=data; 
  datatype = dtype;
  isarray = 0;
};


void DataItem::echodata()
{
  cout << desc << endl;
  // cout << desc << endl;
}; 

int DataItem::actvalue()
{ 
  return (0);
}; 

int DataItem::compare(DataItem* item2)
{ 
  return (desc.matches(item2->desc));
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


DataItem* DataItem::NewDataItem(dtypes newditem, String d='?', String z='?')
{
  DataItem* tmpdata;


  switch(newditem)           
    {
              
    case errd:      
    case comn:
    case comd:   
      
      if (d=='?')     {  tmpdata = new DataItem(newditem)    ;}
      else if (z=='?'){  tmpdata = new DataItem(d,newditem)  ;}
      else            {  tmpdata = new DataItem(d,z,newditem);} ;

      break;
      
    case strg:     
    case stra:
    case strd:
    case strp:

      if (d=='?')     {  tmpdata = new String_DataItem(newditem)    ;} 
      else if (z=='?'){  tmpdata = new String_DataItem(d,newditem)  ;}
      else            {  tmpdata = new String_DataItem(d,z,newditem);};

      break;
                              
    case intd:       /* 1 integer */
			      
      if (d=='?')     {  tmpdata = new Int_DataItem(newditem)    ;}
      else if (z=='?'){  tmpdata = new Int_DataItem(d,newditem)  ;}
      else            {  tmpdata = new Int_DataItem(d,z,newditem);};

      break;
      
    case reald:       /* 1 real*/
      
      if (d=='?')     {  tmpdata = new Real_DataItem(newditem)    ;} 
      else if (z=='?'){  tmpdata = new Real_DataItem(d,newditem)  ;}
      else            {  tmpdata = new Real_DataItem(d,z,newditem);};

      break;
      
    case real2d:       /* 2 real */
      
      if (d=='?')     {  tmpdata = new Real2_DataItem(newditem)    ;}
      else if (z=='?'){  tmpdata = new Real2_DataItem(d,newditem)  ;}
      else            {  tmpdata = new Real2_DataItem(d,z,newditem);};

      break;
      
    case real3d:       /* 3 real*/
      
      if (d=='?')     {  tmpdata = new Real3_DataItem(newditem)    ;}
      else if (z=='?'){  tmpdata = new Real3_DataItem(d,newditem)  ;}
      else            {  tmpdata = new Real3_DataItem(d,z,newditem);};

      break;
      
    default:
      
      if (d=='?')     {   tmpdata = new DataItem(newditem)    ;}
      else if (z=='?'){   tmpdata = new DataItem(d,newditem)  ;}
      else            {   tmpdata = new DataItem(d,z,newditem);};

      break; 
      
    };
  
  return tmpdata; 

};



DataItem* DataItem::NewDataItem(DataItem* original)
{
  DataItem* tmpdata; 
  dtypes tmpdatatype;

  tmpdatatype = original->datatype; 

  tmpdata = DataItem::NewDataItem(original->datatype,original->desc,original->value_string);

  return tmpdata;
 
};


  
DataItem* DataItem::analyse_line(String data)
{
        
  DataItem* tmpdata;
  String d,z;
  int num;
  dtypes ind; 

  if (data.contains("'",0)!=0)

    /* Data line - all else are comments or data*/

    {
      num = data.index("'",1);

      d = data.before(num+1);    
      z = data.after(num);             

      ind = idtype(z);

      tmpdata = DataItem::NewDataItem(ind,d,z);
    } 
  else if (data.contains("$",0)!=0)
    {
      /* Comment line*/
               
      tmpdata = DataItem::NewDataItem(comd,data);

    }   
  else if ((data.contains("c",0)!=0)||(data.contains("C",0)!=0))
    {
      /* NIMBUS Comment line - invalid DIVIMP input */
               
      tmpdata = DataItem::NewDataItem(comn,data);

    }   
  else 
    {
      /* Line containing just data OR an error condition */

      tmpdata = DataItem::NewDataItem(errd,data);
      
    };

  return tmpdata;  
          
};



dtypes DataItem::idtype(String data)
{

  dtypes retval; 

  /* Determines the type of data being passed to it
     0 = error or blank string
     1 = character string
     2 = integer
     3 = 1 real
     4 = 2 reals
     5 = 3 reals  */

  if (data.index("'",0)!=-1)
    {
      /* Contains charater string */

      retval = strg;

    } 
  else if (data.index(".",0)!=-1)
    {
      /* Real numbers ... */

      switch (data.freq("."))
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
  else if (data.contains(RXint)!=0) 
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

  if (desc.contains("' '",0)!=0&&item2.desc.contains("' '",0)!=0) 
    { 
      // If the description is a short empty string then compare the 
      // value_string for equality testing.
      // Strip any trailing spaces and compare the resulting substrings

      int len1 = value_string.length() -1 ;
      int len2 = item2.value_string.length() -1;

      while (len1>=0&&value_string[len1]==' ') len1--;
 
      while (len2>=0&&item2.value_string[len2]==' ') len2--;

      if (len1!=len2) return FALSE;        

      String desc1 = value_string.before(len1+1);
      String desc2 = item2.value_string.before(len2+1); 

      return ((desc1.matches(desc2)) ? TRUE : FALSE );

      /* return ((value_string.matches(item2.value_string)) ? TRUE : FALSE );*/
    }
  else
    {
      // Otherwise - compare the description 
      // Strip any trailing spaces and compare the resulting substrings

      int len1 = desc.length() -1 ;
      int len2 = item2.desc.length() -1;

      while (len1>=0&&desc[len1]==' ') len1--;
 
      while (len2>=0&&item2.desc[len2]==' ') len2--;

      if (len1!=len2) return FALSE;        

      String desc1 = desc.before(len1+1);
      String desc2 = item2.desc.before(len2+1); 

      return ((desc1.matches(desc2)) ? TRUE : FALSE );
    };

}; 




// -------------------------------------------------
// String_DataItem -
// -------------------------------------------------



String_DataItem::String_DataItem() : DataItem()
{
};

String_DataItem::String_DataItem(dtypes dtype) : DataItem(dtype)
{
};

String_DataItem::String_DataItem(String description,dtypes dtype) : DataItem(description,dtype)
{
};

String_DataItem::String_DataItem(String description,String data, dtypes dtype) : DataItem(description,data,dtype)
{
};


void String_DataItem::echodata()
{
  cout << desc << value_string << endl;

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

Int_DataItem::Int_DataItem(String description,dtypes dtype) : DataItem(description,dtype)
{
};

Int_DataItem::Int_DataItem(String description,String data, dtypes dtype) : DataItem(description,data,dtype)
{
  value=atoi((const char *)data);
};


void Int_DataItem::echodata()
{
  cout << desc << value_string << endl;
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

Real_DataItem::Real_DataItem(String description,dtypes dtype) : DataItem(description,dtype)
{
};

Real_DataItem::Real_DataItem(String description,String data, dtypes dtype) : DataItem(description,data,dtype)
{
  value=atoff((const char *)data);
};


void Real_DataItem::echodata()
{
  cout << desc << value_string << endl;
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

Real2_DataItem::Real2_DataItem(String description,dtypes dtype) : DataItem(description,dtype)
{
};

Real2_DataItem::Real2_DataItem(String description, String data,dtypes dtype) : DataItem(description,data,dtype)
{
  istrstream in_string(data,data.length());
  in_string >> value1 >> value2;

};


void Real2_DataItem::echodata()
{
  cout << desc << value_string << endl;
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

Real3_DataItem::Real3_DataItem(String description,dtypes dtype) : DataItem(description,dtype)
{
};


Real3_DataItem::Real3_DataItem(String description,String data,dtypes dtype) : DataItem(description,data,dtype)
{
  istrstream in_string(data,data.length());
  in_string >> value1 >> value2 >> value3;
};


void Real3_DataItem::echodata()
{
  cout << desc << value_string << endl;
}; 







