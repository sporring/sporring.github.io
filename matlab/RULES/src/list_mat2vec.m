function v = list_mat2vec(m);
%LIST_MAT2VEC Transform an N by M matric into a 2+N*M vector
%       
%       v = list_mat2vec(m)
%         v - the vector representation
%         m - the original matrix
%
%       list_mat2vec uses reshape to transform a matrix into a vector
%       and prepends the original dimensions for later use by list_vec2mat.
%       An empty matrix returns an empty vector.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, May 23, 1997

if isempty(m)
  v = [];
else
  v = [size(m),reshape(m,1,prod(size(m)))];
end
disp(sprintf('list_mat2vec: %d, %d -> %d, %d',size(m),size(v)));
