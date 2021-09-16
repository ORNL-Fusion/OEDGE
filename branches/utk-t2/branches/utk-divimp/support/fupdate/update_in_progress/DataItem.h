class DataItem
{

protected:
   
  string tag;
  tagtypes tagtype;
  string desc;
  string value_string;
  dtypes datatype; 
  int isarray;
  int copystatus;

public:

  DataItem();
  DataItem(dtypes);
  DataItem(string );
  DataItem(string,dtypes);
  DataItem(string,string,dtypes);

  virtual void echodata();
  virtual int actvalue(); 

  int compare(DataItem* );

  int istype(dtypes);

  void settype(dtypes);
  dtypes gettype();

  void setarray(int);  
  int  getarray();  

  string gettag();
  tagtypes gettagtype();
  void extracttag();


  //static DataItem* NewDataItem(dtypes,string d='?',string z='?');

  static DataItem* NewDataItem(dtypes);
  static DataItem* NewDataItem(dtypes,string);
  static DataItem* NewDataItem(dtypes,string,string);

  static DataItem* NewDataItem(DataItem *);

  static DataItem* analyse_line(string);

  static dtypes idtype(string);

  int operator == (DataItem);

};



class string_DataItem : public DataItem
{
public:

  string_DataItem();
  string_DataItem(dtypes);
  string_DataItem(string,dtypes);
  string_DataItem(string,string, dtypes);

  void echodata();

};


class Int_DataItem : public DataItem 
{

protected:

  int value;

public:

  Int_DataItem();
  Int_DataItem(dtypes);
  Int_DataItem(string,dtypes);
  Int_DataItem(string,string, dtypes);

  void echodata();
  int actvalue();

};
  

class Real_DataItem : public DataItem 
{

protected:

  float value;

public:


  Real_DataItem();
  Real_DataItem(dtypes);
  Real_DataItem(string,dtypes);
  Real_DataItem(string,string, dtypes);

  void echodata();
};


class Real2_DataItem : public DataItem 
{

protected:

  float value1,value2;

public:


  Real2_DataItem();
  Real2_DataItem(dtypes);
  Real2_DataItem(string,dtypes);
  Real2_DataItem(string, string,dtypes);

  void echodata();
 };  


class Real3_DataItem : public DataItem 
{

protected:

  float value1,value2,value3;

public:

  Real3_DataItem();
  Real3_DataItem(dtypes);
  Real3_DataItem(string,dtypes);
  Real3_DataItem(string, string,dtypes);

  void echodata();


 }; 









