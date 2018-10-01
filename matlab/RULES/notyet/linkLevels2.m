function [map, newlinkCell] =  linkLevels2(linkCell, top, bottom);

% function [map, linkCell] =  linkLevels2(linkCell, top, bottom);
%
%				Ole Fogh Olsen, March 1998

if(top == bottom)
  nrRegions = linkCell{top,bottom};
  map = sparse(1:nrRegions,1:nrRegions,ones(1,nrRegions));
  newlinkCell = linkCell;
  return;
end

% invers map are the transposed of the map
if(top > bottom)
  row = top;
  col = bottom;
else
  row = bottom;
  col = top;
end;

% Note row > bottom always

% Calculate map if not in store
% Speed up later calls 
% by keeping intermediate calculated maps
if( isempty(linkCell{row,col}) )    
  cur_row = row-1;
  while(isempty(linkCell{cur_row,col}))
    cur_row = cur_row -1;
  end
  for i=cur_row+1:row,
    linkCell{i,col} =  linkCell{i,i-1} * linkCell{i-1,col};
  end
end

if(top < bottom)
  map = transpose(linkCell{row,col});
else
  map = linkCell{row,col};
end

newlinkCell = linkCell;




