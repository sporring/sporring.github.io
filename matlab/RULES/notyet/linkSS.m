function linkCell = linkSS(SS, linkfunctionname);

% function linkCell = linkSS(SS, linkfunctionname);
% A linking structure for the scale space SS is created
% the function linkfunctionname is used to link neighbouring levels
% Only linking between neighbouring levels is made
%
%			Ole Fogh Olsen, March 1998


[m n l] = size(SS);
linkCell = cell(l, l);

for i=1 : l-1,
  linkCell{i+1,i} = feval(linkfunctionname, SS(:,:,i+1), SS(:,:,i));
  %
  [m n] = size(linkCell{i+1,i});
  linkCell{i,i} = [n];
end

[m n] = size(linkCell{l,l-1});
linkCell{l,l} = [m];

return;


