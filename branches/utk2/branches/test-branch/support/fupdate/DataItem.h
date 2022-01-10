class DataItem
{

protected:
   
  String desc;
  String value_string;
  dtypes datatype; 
  int isarray;

public:

  DataItem();
  DataItem(dtypes);
  DataItem(String );
  DataItem(String,dtypes);
  DataItem(String,String,dtypes);

  virtual void echodata();
  virtual int actvalue(); 

  int compare(DataItem* );

  int istype(dtypes);

  void settype(dtypes);
  dtypes gettype();

  void setarray(int);  
  int  getarray();  

  static DataItem* NewDataItem(dtypes,String d='?',String z='?');

  // static DataItem* NewDataItem(dtypes,String);
  // static DataItem* NewDataItem(dtypes,String,String);

  static DataItem* NewDataItem(DataItem *);
  static DataItem* analyse_line(String);

  static dtypes idtype(String);

  int operator == (DataItem);

};



class String_DataItem : public DataItem
{
public:

  String_DataItem();
  String_DataItem(dtypes);
  String_DataItem(String,dtypes);
  String_DataItem(String,String, dtypes);

  void echodata();

};


class Int_DataItem : public DataItem 
{

protected:

  int value;

public:

  Int_DataItem();
  Int_DataItem(dtypes);
  Int_DataItem(String,dtypes);
  Int_DataItem(String,String, dtypes);

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
  Real_DataItem(String,dtypes);
  Real_DataItem(String,String, dtypes);

  void echodata();
};


class Real2_DataItem : public DataItem 
{

protected:

  float value1,value2;

public:


  Real2_DataItem();
  Real2_DataItem(dtypes);
  Real2_DataItem(String,dtypes);
  Real2_DataItem(String, String,dtypes);

  void echodata();
 };  


class Real3_DataItem : public DataItem 
{

protected:

  float value1,value2,value3;

public:

  Real3_DataItem();
  Real3_DataItem(dtypes);
  Real3_DataItem(String,dtypes);
  Real3_DataItem(String, String,dtypes);

  void echodata();


 }; 









