#include <iostream>
#include <sstream>

using namespace std;

#include "Local_Types.h"
#include "DataItem.h"
#include "DataFile.h"
#include "LineItem.h"

// -------------------------------------------------
// LineItem 
// -------------------------------------------------

LineItem::LineItem()
{

  rawdata = '?';
  procdata = NULL;

};

LineItem::LineItem(string line)
{

  rawdata = line;
  procdata = DataItem::analyse_line(rawdata[i]);  

};
