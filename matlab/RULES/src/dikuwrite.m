function dikuwrite(I,name,type);
%DIKUWRITE Write a matrix or vector as a DIKU image.
%
%       dikuwrite(I,name,type)
%         I - a matrix to be saved
%         name - the name of the diku image file
%         type - an optional image type of 
%                ('grey', 'gray', 'short', 'int', 'float', 'double')
%
%       DIKU images are written (and perhaps converted) with
%       dikuwrite.  The default type is 'double'.  Only scalar images
%       (matrices) are supported (not color).
%
%                                          Jon Sporring, March 4, 1998
