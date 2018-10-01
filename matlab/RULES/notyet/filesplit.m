function filesplit(filename,max_size)
% Split files into smaller pieces e.g. for disk-transport.
  
  if nargin < 2
    max_size = 1457664;
  end
  fr = fopen(filename,'r');
  if fr ~= -1
    count = max_size
    i = 1;
    while count == max_size
      [t,count] = fread(fr,max_size,'uchar');
      disp(count)
      fw = fopen(sprintf('%s%d',filename,i),'w');
      fwrite(fw,t,'uchar');
      fclose(fw);
      i = i+1;
    end
  end
  fclose(fr);
