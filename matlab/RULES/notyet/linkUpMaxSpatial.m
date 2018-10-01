function linkMatrix = linkUpMaxSpatial(top, bottom);

% linkMatrix = linkUpMaxSpatial(top, bottom);

% top and down are integer "images" with the same dimensions;
% Pixels with same integer value belong to the same region.
% The function linkUpMaxSpatial makes a many-to-one linking 
% between  top regions and bottom region.
% A bottom region links to the top region with which it has 
% maximum overlap.


maxtop = max(max(top));
maxbottom = max(max(bottom));

overlap = sparse(maxtop, maxbottom);

% index = row# + (col#-1) * rowLength 
indexInOverlap = top + (bottom-1) * maxtop;
% reshape to vector
indexInOverlap = reshape(indexInOverlap, 1, prod(size(indexInOverlap)));

for i=indexInOverlap
  overlap(i) = overlap(i) + 1;
end;


% Select the maximum in each column
% that is for each bottom segment find the top segment with the greatest overlap
[Y I] = max(overlap);

linkMatrix = sparse(I, 1:maxbottom, ones(1,maxbottom), maxtop, maxbottom, maxbottom);



