#include <iostream>
#include <fstream>

using namespace std;

#include "Local_Types.h"
#include "DataItem.h"
#include "DataFile.h"


// C++ code for class DataFile

DataFile::DataFile()
{ 
  num_lines=0; 
  filetype = 0;
  fileindex = 0;
  for (int i = 0; i<max_lines;i++)
  {
    procdata[i]=NULL;   
  }; 
  
};

DataFile::~DataFile()
{
  for (int i = 0; i<num_lines;i++)
  {
    if (procdata[i]!=NULL) delete procdata[i];   
  }; 
};



void DataFile::read(char* filename)
{
  ifstream infile;
  char* s;
  num_lines = 0;

  cout << "READ:" << filename << "\n"; 

  infile.open(filename);

  while (!infile.eof()&&!infile.fail())
    {
      infile.getline(buffer,max_length);

      if (num_lines<max_lines)  
	{       
	  rawdata[num_lines] = buffer; 
	  num_lines++;

	};

    };
  
  infile.close();

  // Set fileindex to the end of the file - that is where any new lines would be added
  fileindex = num_lines;

  cout << "FINISHED READ: RETURN CODE=" << infile.fail() << endl;  

};


void DataFile::print() 
{
  for (int i = 0; i<num_lines;i++)
    {
      procdata[i]->echodata(); 
    };
};


void DataFile::write(char* filename)
{
  ofstream outfile;
  char* s;

  outfile.open(filename);
  
  for (int i = 0; i <num_lines;i++)
    {   
      outfile << rawdata[i] << "\n";
    }; 

  outfile.close();
};

void DataFile::split_file()
{
  int isarray=0,startarray=0,numarray=0,cnt=0,arrcnt=1,i,j;
                   
  for (i=0;i<num_lines;i++)
    {
      // cout << "ANALYSING:" << i << ":" << rawdata[i] << endl;

      procdata[i] = DataItem::analyse_line(rawdata[i]);

      
      // Analyse the processed line to determine
      // if it is part of an Array DataItem 


      if ((procdata[i]->istype(stra)==TRUE)||(procdata[i]->istype(strg)==TRUE))
	{
	  // Deal with pre-existing arrays before starting this NEW one

          if (isarray==1) 
	    {
	      cnt = 0;      
	      procdata[startarray]->settype(strd);

	    }
	  else if (isarray==2)
	    {
	      if (cnt==(procdata[numarray]->actvalue()))
		{
		  for (j = startarray;j < i;j++)
		    procdata[j]->setarray(arrcnt);

		  arrcnt++;
	      
		  if (cnt==0) 
		    {
		      if ((procdata[startarray]->istype(strg)==TRUE))
			{
			  procdata[startarray]->settype(strp); 
			};
		    }
		  else
		    {		
		      procdata[startarray]->settype(stra);
		    };

		  isarray = 0;
		  cnt = 0;
		}
	      else
		{
		  /* This should be an error condition that is not encountered where
		     the number of data items read isn't equal to the number specified
		     on the integer input line - unless the number of items read is 0 for
		     a non-zero integer value */
		  
		  isarray = 0;
		  cnt = 0;
		  procdata[startarray]->settype(strd);
		};
	      
	      
	    };

	  isarray = 1;
	  startarray=i;
	}
      else if (isarray==1&&(procdata[i]->istype(intd)==TRUE))
	{
	  isarray=2; 
	  numarray=i;
	  cnt=0; 
	}
      else if (isarray==2&&(procdata[i]->istype(errd)==TRUE))
	{
	  cnt = cnt + 1 ;
	}
      else if (isarray==2&&(procdata[i]->istype(comd)==FALSE))
	{
	  if (cnt==(procdata[numarray]->actvalue()))
	    {
	      for (j = startarray;j < i;j++)
		procdata[j]->setarray(arrcnt);

              arrcnt++;
	      
              if (cnt==0) 
		{
		  if ((procdata[startarray]->istype(strg)==TRUE))
		    {
		      procdata[startarray]->settype(strp); 
		    };
		}
	      else
                {		
		  procdata[startarray]->settype(stra);
		};

              isarray = 0;
              cnt = 0;
	    }
	  else
	    {
              /* This should be an error condition that is not encountered where
                 the number of data items read isn't equal to the number specified
                 on the integer input line - unless the number of items read is 0 for
                 a non-zero integer value */

	      isarray = 0;
	      cnt = 0;
              procdata[startarray]->settype(strd);
	    };
	  
	}    
      else if (isarray==1&&(procdata[i]->istype(comd)==FALSE)) 
	{
	  isarray = 0;
          cnt = 0;      
          procdata[startarray]->settype(strd);
	};
             
      cout << "PROC:" << i << ":" << isarray << ":" << cnt << ":" << arrcnt << ":" procdata[i]->gettype() << ":";
      procdata[i]->echodata();
      


    }; 




  
};
  

int DataFile::search(DataItem test,int n1=0)
{
  // Check to see if the file contains a matching DataItem - description - not value.
  // Check the range of values listed.

  // If a match is found return it's index.

  for (int i=n1;i<num_lines;i++)
    {
      if(*procdata[i]==test) return i;
    };

  return -1;

};

void DataFile::add_item(int dest,int src, DataFile srcfile );
{

  if ( 

  rawdata[index]=raw;
  procdata[index]=DataItem::NewDataItem(proc);
  
};

void DataFile::add_item(DataFile datasource, int dataindex=-1);
{
  if (dataindex=-1) dataindex=fileindex;


	      rawdata[dataindex]  = datasource.rawdata[cntold]; 
	      
	      procdata[dataindex] = DataItem::NewDataItem(oldf.procdata[cntold]);


};

void DataFile::setfileindex(int newindex)
{
  fileindex = newindex;

};

int DataFile::getfileindex()
{
  getfileindex = fileindex;

};






DataFile* DataFile::update(DataFile& oldf)
{

  // Compares the current to the specified and generates an updated output
       
  DataFile* outp;

  int i; 
  int found_index;
  int arrnum,ind;
  int cntnew = 0;
  int cntold = 0;
  int cntout = 0; 
  int tmpcnt = 0;
  dtypes newdtype,olddtype;

  outp = new DataFile();           

  cout << "Starting Update ... " << num_lines << " " << oldf.num_lines <<  endl;


  // Scan OLD for tagged optional input not present in TEMPLATE - add to new at the beginning for easy review

  while (cntold < oldf.num_lines)
    {

      // Optional input found
      if (oldf.procdata[cntold]->tagtype==opt)
	{
	  // Optional input not found in TEMPLATE
	  if (search(*oldf.procdata[cntold])==-1)
	    {

	      // Add optional input to new file
	      outp->add_item(cntnew,cntold,oldf));


	    };


	};


    };


  // Scan through the template - copy all required lines - check to see if these lines exist in the OLD file and use their
  // data values instead. 

  while (cntnew < num_lines)
    {

      
      







    };



//      cout << "NEW:"<<cntnew << ":" << rawdata[cntnew] << endl;  
//      procdata[cntnew]->echodata();
//      cout << "OLD:"<<cntold << ":" << oldf.rawdata[cntold] << endl ;
//      oldf.procdata[cntold]->echodata();


//      if(cntout!=0) 
//        {

//          cout << "OUT:"<< cntout-1 << ":" << outp->rawdata[cntout-1] << endl ;
//          cout << "OUT:" << 0 << ":" << outp->rawdata[0] << endl;
//          outp->procdata[cntout-1]->echodata();
//        } ;

      // Check to see if the OLD file counter is pointing at
      // a comment - if it is - copy this to the output and update
      // the OLD and OUT counters.

      dtypes tmp1 = procdata[cntnew]->gettype();
      dtypes tmp2 = oldf.procdata[cntold]->gettype();



      


      // Search old for a match to the current line or array
      found_index=oldf.search(*procdata[cntnew]);

      // Line not found in OLD
      if (found_index==-1) 
	{
	  // Copy input from TEMPLATE TO NEW




	  cntnew = next_item(cntnew);
	}

      else
	{
	  // Match found in OLD
	  // Copy line from OLD to NEW




	  cntnew = next_item(cntnew);
	}




      
      // Comments in old

      if (oldf.procdata[cntold]->istype(comd)) 
	{
	  
	  // Search through NEW to see if comments will already be included. 
	  // If they are NEW comments in old copy comments over from new to output
	  
	  if (search(*oldf.procdata[cntold])==0) 
	    {
	      outp->rawdata[cntout]  = oldf.rawdata[cntold]; 
	      
	      // Create and assign output object
	      
	      outp->procdata[cntout] = DataItem::NewDataItem(oldf.procdata[cntold]);

	      cntout++;
	      
	    };
	     
	  cntold++; 
	  
	}

      // Comments in new 
    
      else if (procdata[cntnew]->istype(comd)) 

	{

	  // Copy comments from new directly into output

	  outp->rawdata[cntout]  = rawdata[cntnew]; 
	  
	  // Create and assign output object
	  
	  outp->procdata[cntout] = DataItem::NewDataItem(procdata[cntnew]);

	  cntout++;
          cntnew++;
             
	}

      // Reached NIMBUS code in the new input file - all else processed.
      // Discard OLD until nimbus reached then copy OLD NIMBUS to output.

      else if (procdata[cntnew]->istype(comn)) 

	{
//              cout << "CNT1:" << cntnew << " " << cntold << " " << cntout << endl;
//              cout << "CNT2:" << num_lines << " " << oldf.num_lines << " " << outp->num_lines << endl;

              // Discard up to NIMBUS input
	  
              

              while(oldf.procdata[cntold]->istype(comn)==FALSE&&
                    cntold<=oldf.num_lines)
		{	
		  cout << "N-DISCARDED: ";
		  oldf.procdata[cntold]->echodata();
		  cntold++;
		};
      
	      // If at end of old - assume it had no NIMBUS data - use new
	      
              if (cntold>oldf.num_lines) 
		{
		  
		  while (cntnew<=num_lines)
		    {
		      // Copy NIMBUS input comments from new directly into output
		      
		      outp->rawdata[cntout]  = rawdata[cntnew]; 
		      
		      // Create and assign output object
		  
		      outp->procdata[cntout] = DataItem::NewDataItem(procdata[cntnew]);
		  
		      cntout++;
		      cntnew++;
		    };
		}

	      else

		{

		  while (cntold<=oldf.num_lines)
		    {
		      // Copy NIMBUS input comments from old directly into output
		      
		      outp->rawdata[cntout]  = oldf.rawdata[cntold]; 
		      
		      // Create and assign output object
		  
		      outp->procdata[cntout] = DataItem::NewDataItem(oldf.procdata[cntold]);
		  
		      cntout++;
		      cntold++;
		    };
 
//		  cout << "Exiting ... " << endl;
//                  cout << "OUT:" << cntout << ":" << outp->rawdata[0] << endl;
		  cntnew = num_lines;
		};


	      
	}
      
    
      // Does new line equal old line - use old.

      else if (*procdata[cntnew]==*oldf.procdata[cntold]) 

	{

	  // If new line and old have the same description - use old
	  // Check to see if the matching line is a string data item - then 
	  // need to check for possible arrays ...
	  
	  // Check for string data values - these match.

	  if (procdata[cntnew]->istype(stra)||
	      procdata[cntnew]->istype(strp)||
	      procdata[cntnew]->istype(strd)) 
	    { 
	      
	      newdtype = procdata[cntnew]->gettype();
	      olddtype = oldf.procdata[cntold]->gettype();
		   
	      if (newdtype==stra&&(olddtype==stra||olddtype==strp)) 
		{
		  // Check to see if the integer descriptions are also the same

		  if (*procdata[cntnew+1]==*oldf.procdata[cntold+1])
		    {
		      // Same - use OLD

		      arrnum = oldf.procdata[cntold]->getarray();
		      ind = arrnum;
		  
		      while (ind==arrnum) 
			{
		      
			  outp->rawdata[cntout]  = oldf.rawdata[cntold]; 
			  outp->procdata[cntout] = DataItem::NewDataItem(oldf.procdata[cntold]);
			  cntout++;
			  cntold++;
			  ind = oldf.procdata[cntold]->getarray();
			};
		      
		      arrnum = procdata[cntnew]->getarray();
		      ind = arrnum;
		      
		      while (ind==arrnum) 
			{
			  cntnew++;
			  ind = procdata[cntnew]->getarray();
			};
		      
		    }
		  else
		    {
		      // Different - use NEW


		      arrnum = procdata[cntnew]->getarray();
		      ind = arrnum;
		  
		      while (ind==arrnum) 
			{
		      
			  outp->rawdata[cntout]  = rawdata[cntnew]; 
			  outp->procdata[cntout] = DataItem::NewDataItem(procdata[cntnew]);
			  cntout++;
			  cntnew++;
			  ind = procdata[cntnew]->getarray();
			};
		      
		      arrnum = oldf.procdata[cntold]->getarray();
		      ind = arrnum;
		      
		      while (ind==arrnum) 
			{
			  cntold++;
			  ind = oldf.procdata[cntold]->getarray();
			};

		    };

		}

	      else if (newdtype==strp&&olddtype==stra)
		{

		  // Use the OLD array entry - either stra or strp - pop off 
		  // the new file entries.
		  

		  arrnum = oldf.procdata[cntold]->getarray();
		  ind = arrnum;
		  
		  while (ind==arrnum) 
		    {
		      
		      outp->rawdata[cntout]  = oldf.rawdata[cntold]; 
		      outp->procdata[cntout] = DataItem::NewDataItem(oldf.procdata[cntold]);
		      cntout++;
		      cntold++;
		      ind = oldf.procdata[cntold]->getarray();
		    };
		  
		  arrnum = procdata[cntnew]->getarray();
		  ind = arrnum;
		  
		  while (ind==arrnum) 
		    {
		      cntnew++;
		      ind = procdata[cntnew]->getarray();
		    };
		}
	      else if ((newdtype==strp||newdtype==strd)&&(olddtype==strp||olddtype==strd))
		{
		       
		  // Copy rawdata from old directly into output
		  
		  outp->rawdata[cntout]  = oldf.rawdata[cntold]; 
		  
		  // Create and assign output object
		  
		  outp->procdata[cntout] = DataItem::NewDataItem(oldf.procdata[cntold]);
		  
		  cntout++;
		  cntnew++;
		  cntold++;
		}
	      else if (newdtype==strd&&olddtype==stra)
		{
		  
		  // Copy rawdata from new directly into output
		  
		  outp->rawdata[cntout]  = rawdata[cntnew]; 
		  
		  // Create and assign output object
		  
		  outp->procdata[cntout] = DataItem::NewDataItem(procdata[cntnew]);
		  
		  cntout++;
		  cntnew++;
		  
		  // Eliminate OLD array
		  
		  arrnum = oldf.procdata[cntold]->getarray();
		  ind = arrnum;
		  
		  while (ind==arrnum) 
		    {
		      cntold++;
		      ind = oldf.procdata[cntold]->getarray();
		    };
		}
	      else if (newdtype==stra&&olddtype==strd)
		{
		  // Use new array and discard old strd
		  
		  arrnum = procdata[cntnew]->getarray();
		  ind = arrnum;
		  
		  while (ind==arrnum) 
		    {
		      
		      outp->rawdata[cntout]  = rawdata[cntnew]; 
		      outp->procdata[cntout] = DataItem::NewDataItem(procdata[cntnew]);
		      cntout++;
		      cntnew++;
		      ind = procdata[cntnew]->getarray();
		    };
		  
		  
		  cntold++;
		  
		  
		};
	      
	    }
	 

	  else
	    {
	      
	      // Copy matching line from OLD

	      // IF it is part of an array - copy over all remaining elements
	      //                             of the array as well    


	      arrnum = oldf.procdata[cntold]->getarray();
	      ind = arrnum;

	      if (ind!=0) 
		{ 
		      
		  while (ind==arrnum) 
		    {
			  
		      outp->rawdata[cntout]  = oldf.rawdata[cntold]; 
		      outp->procdata[cntout] = DataItem::NewDataItem(oldf.procdata[cntold]);
		      cntout++;
		      cntold++;
		      ind = oldf.procdata[cntold]->getarray();
		    };
		}
	      else
		{
		  
		  // Copy rawdata from old directly into output
		  
		  outp->rawdata[cntout]  = oldf.rawdata[cntold]; 
		  
		  // Create and assign output object
		  
		  outp->procdata[cntout] = DataItem::NewDataItem(oldf.procdata[cntold]);
		  
		  cntout++;
		  cntold++;
		  
		};
	      
	      cntnew++;

	      
	    };
	  
	} 
      
      // Does not match - search for matching line - if one exists.

      else 
	
	{
	  // Search through Old file for a line matching the new one 
	  // This will also serve to skip past lines whose descriptions
	  // may have changed and replace them with updated versions. 
	  
	  
	  tmpcnt = oldf.search(*procdata[cntnew]);
	  
	  if (tmpcnt==0) 
	    {
	      
	      cout << "ADDED    : ";

	      if (procdata[cntnew]->istype(stra))
		{
		  
		  arrnum = procdata[cntnew]->getarray();
		  ind = arrnum;
		  
		  while (ind==arrnum) 
		    {
		      
		      outp->rawdata[cntout]  = rawdata[cntnew]; 
		      outp->procdata[cntout] = DataItem::NewDataItem(procdata[cntnew]);

                      outp->procdata[cntout]->echodata();

		      cntout++;
		      cntnew++;
		      ind = procdata[cntnew]->getarray();
		    };

		}

	      else
		{
		  outp->rawdata[cntout]  = rawdata[cntnew]; 
		  
		  // Create and assign output object
		  
		  outp->procdata[cntout] = DataItem::NewDataItem(procdata[cntnew]);
		  
                  outp->procdata[cntout]->echodata();

		  cntout++;
		  cntnew++;
		  
		};
	    }
	  
	  else
	    {


	      // If a match is found for optional input lines then the intervening old input 
	      // lines can not be discarded - they must be kept and checked.
	      // In fact optional input lines should never be discarded. 

	      if (tagtype==opt) 
		{





		}

	      else
		{
              cout << "DISCARDED: ";
	      
	      for (i=cntold;i<tmpcnt;i++)
		{
		  oldf.procdata[i]->echodata();
                  if (i!=(tmpcnt-1)) cout << "           ";
		};
	      
	      cntold = tmpcnt;
	      
	      // Check for string data values - these match.
	      
	      if (procdata[cntnew]->istype(stra)||
		  procdata[cntnew]->istype(strp)||
		  procdata[cntnew]->istype(strd)) 
		{ 
		  
		  newdtype = procdata[cntnew]->gettype();
		  olddtype = oldf.procdata[cntold]->gettype();
		  
		  if ((newdtype==stra&&(olddtype==stra||olddtype==strp))||
		      (newdtype==strp&&olddtype==stra))  
		    {
		      // Use the OLD array entry - either stra or strp - pop off 
		      // the new file entries.
		      
		      
		      arrnum = oldf.procdata[cntold]->getarray();
		      ind = arrnum;
		      
		      while (ind==arrnum) 
			{
			  
			  outp->rawdata[cntout]  = oldf.rawdata[cntold]; 
			  outp->procdata[cntout] = DataItem::NewDataItem(oldf.procdata[cntold]);
			  cntout++;
			  cntold++;
			  ind = oldf.procdata[cntold]->getarray();
			};
		      
		      arrnum = procdata[cntnew]->getarray();
		      ind = arrnum;
		      
		      while (ind==arrnum) 
			{
			  cntnew++;
			  ind = procdata[cntnew]->getarray();
			};
		    }
		  else if ((newdtype==strp||newdtype==strd)&&(olddtype==strp||olddtype==strd))
		    {
		      
		      // Copy rawdata from old directly into output
		      
		      outp->rawdata[cntout]  = oldf.rawdata[cntold]; 
		      
		      // Create and assign output object
		      
		      outp->procdata[cntout] = DataItem::NewDataItem(oldf.procdata[cntold]);
		      
		      cntout++;
		      cntnew++;
		      cntold++;
		    }
		  else if (newdtype==strd&&olddtype==stra)
		    {
		      
		      // Copy rawdata from new directly into output
		      
		      outp->rawdata[cntout]  = rawdata[cntnew]; 
		      
		      // Create and assign output object
		      
		      outp->procdata[cntout] = DataItem::NewDataItem(procdata[cntnew]);
		      
		      cntout++;
		      cntnew++;
		      
		      // Eliminate OLD array
			   
		      arrnum = oldf.procdata[cntold]->getarray();
		      ind = arrnum;
		      
		      while (ind==arrnum) 
			{
			  cntold++;
			  ind = oldf.procdata[cntold]->getarray();
			};
		    }
		  else if (newdtype==stra&&olddtype==strd)
		    {
		      // Use new array and discard old strd
		      
		      arrnum = procdata[cntnew]->getarray();
		      ind = arrnum;
		      
		      while (ind==arrnum) 
			{
			  
			  outp->rawdata[cntout]  = rawdata[cntnew]; 
			  outp->procdata[cntout] = DataItem::NewDataItem(procdata[cntnew]);
			  cntout++;
			  cntnew++;
			  ind = procdata[cntnew]->getarray();
			};
		      
		      
		      cntold++;
		      
		      
		    };
		  
		}
	      else
		{
		  
		  // IF it is part of an array - copy over all remaining elements
		  //    of the array as well       *ADD* 


		  arrnum = oldf.procdata[cntold]->getarray();
		  ind = arrnum;

		  if (ind!=0) 
		    { 
		       
		      while (ind==arrnum) 
			{
			  
			  outp->rawdata[cntout]  = oldf.rawdata[cntold]; 
			  outp->procdata[cntout] = DataItem::NewDataItem(oldf.procdata[cntold]);
			  cntout++;
			  cntold++;
			  ind = oldf.procdata[cntold]->getarray();
			};
		    }
		  else
		    {
		      
		      // Copy rawdata from old directly into output
		      
		      outp->rawdata[cntout]  = oldf.rawdata[cntold]; 
		       
		      // Create and assign output object
		       
		      outp->procdata[cntout] = DataItem::NewDataItem(oldf.procdata[cntold]);
		  
		      cntout++;
		      cntold++;
		      
		    };
		  
		  cntnew++;
		   
		};
	      
	    };
	};
      
    };
  
  outp->num_lines = cntout -1;

  return (outp);

};























































