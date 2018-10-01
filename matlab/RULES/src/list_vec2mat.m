function m = list_vec2mat(v);
%LIST_VEC2MAT Transform a vector prepended with the dimensions into a matrix
%       
%       m = list_vec2mat(v)
%         m - the matrix
%         v - the vector representation prepended with the dimensions of m
%
%       list_vec2mat uses reshape to transform a vector prepended with the
%       dimensions into a matrix, as produced by list_mat2vec.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, May 23, 1997

[size(v),v(1:2)]

if size(v,2) < 3
  m = [];
else
  m = reshape(v(3:size(v,2)),v(1:2));
end
disp(sprintf('list_vec2mat: %d, %d -> %d, %d',size(v),size(m)));
