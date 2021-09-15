class LineItem
{
protected:

  string rawdata;
  DataItem* procdata;

public:

  LineItem(string);


  int operator == (LineItem);

};
