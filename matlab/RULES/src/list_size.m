function s = list_size(L);
%LIST_SIZE Returns the number of elements in the list
%       
%       s = list_size(L)
%         s - the number of elements in L.
%         L - the list.
%
%       list_size counts the number of elements in the list L.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, February 12, 1997

s = 0;
n = list_first(L);
while n ~= 0
  s = s+1;
  n = list_next(L,n);
end

