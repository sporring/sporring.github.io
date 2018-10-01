function m = spatialcmoment(I,n,m)
% SPATIALCMOMENT  Calculate spatial central moment of an image.
%

rn = (ones(size(I,2),1)*([1:size(I,1)]-spatialmoment(I,1,0)).^n)'; 
cm = ones(size(I,1),1)*([1:size(I,2)]-spatialmoment(I,0,1)).^m; 

m = sum(sum(I.*rn.*cm))/sum(I(:));
