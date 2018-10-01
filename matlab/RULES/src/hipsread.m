function I = hipsread(name);
%HIPSREAD Read a HIPS image as a matrix or vector.
%
%       I = hipsread(name)
%         I - the resulting matrix
%         name - the name of the hips image file
%
%       WARNING: This function is build over a reverse engineering of the
%       header.  There is NO GUARANTEE that this will work for you image!.
%       However, it has worked for my color images, and you may use the code
%       as an example for writing your own hips-hacks.  HIPS Images are read
%       and converted to double if necessary.
%
%       October 12, 2000:  Header parser changed to read Per Larsens
%       Ear-images.  Jon Sporring.
%
%                                          Jon Sporring, March 4, 1998

