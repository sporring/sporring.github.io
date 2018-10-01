function I = dikuread(name);
%DIKUREAD Read a DIKU greylevel image as a matrix or vector.
%
%       I = dikuread(name)
%         I - the resulting matrix
%         name - the name of the diku image file
%
%       Scalar DIKU Images are read and converted to double if
%       necessary.  Colour images and higher dimensional images are
%       not supported yet.
%
%                                          Jon Sporring, March 4, 1998
