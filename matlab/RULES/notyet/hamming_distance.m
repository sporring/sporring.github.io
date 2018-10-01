function z = hamming_distance(x,y)
%HAMMING_DISTANCE The hamming distance between two binary matrices
%
%       z = hamming_distance(x,y)
%         z - is the vector of Hamming distances
%         x & y - are two vectors of binary numbers (a matrix)
%
%       Hamming_distance returns the hamming distance
%	between two binary matrices, where each entry X(i,:)
%       is interpreted as a binary string of +1,-1.
%
%                                          Jon Sporring, January 1, 1994

z = sum((x~=y)')';
