

class DataFile
{


protected:

   static const int max_lines=5000;
   static const int max_length=1000;
   int num_lines;
   char buffer[max_length];   
   string rawdata[max_lines];
   DataItem *procdata[max_lines]; 
   int filetype;     
   int fileindex;   

public:

   // Constructors

   DataFile();

   // Destructors - if any

   ~DataFile();

   // Functions

   void read(char*);
   void print();
   void write(char*);
   void split_file();
   void add_item(DataItem);
   void add_item(DataFile,int);
   int  search(DataItem,int);

   void setfileindex(int);
   int getfileindex();

   DataFile* update(DataFile&);


};

