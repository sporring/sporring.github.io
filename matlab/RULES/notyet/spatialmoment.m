function m = spatialmoment(I,n,m)
% SPATIALMOMENT  Calculate spatial moment of an image.
%

rn = (ones(size(I,2),1)*[1:size(I,1)].^n)'; 
cm = ones(size(I,1),1)*[1:size(I,2)].^m; 

m = sum(sum(I.*rn.*cm))/sum(I(:));
