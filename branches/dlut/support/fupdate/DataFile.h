

class DataFile
{


protected:

   const int max_lines=1000;
   const int max_length=132;
   int num_lines;
   char buffer[max_length];   
   String rawdata[max_lines];
   DataItem *procdata[max_lines]; 
   int filetype;     


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
   int  search(DataItem,int);

   DataFile* update(DataFile&);


};

