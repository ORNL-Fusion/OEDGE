function [data] = mybinload(filename,Nt,dims,datatype,endian)

   % Open the file
   fid = fopen(filename,'r',endian);

   % Allocates memory for the data
   data = zeros([Nt,dims]);

   % Loop on time steps to read the data
   for it = 1:Nt
      % Dumps the opening character
      dum = fread(fid,1,'uint64');
      % Reads the data and reshapes it
      data(it,:)= reshape(fread(fid,prod(dims),datatype),[1,dims]);
      % Dumps the closing character
      dum = fread(fid,1,'uint64');
   end

   % Close the file
   fclose(fid);

end