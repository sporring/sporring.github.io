function J = stretch(I,l,u);
%STRETCH Do a linear re-mapping of the values in a matrix
%       
%       J = stretch(I,l,u)
%         J - the extended matrix
%         I - the original matrix
%         l - the lower boundary of the new value domain
%         u - the upper boundary of the new value domain
%
%       Stretch does a linear remapping of the values in I according to
%       J(x,y) = l + (I(x,y)-min(I))*(u-l)/(max(I)-min(I)).
%       Behaviour is unspecified if max(I)=min(I).
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, January 13, 1997

minI = min(min(I));

J = l+(I-minI)*(u-l)/(max(max(I))-minI);
