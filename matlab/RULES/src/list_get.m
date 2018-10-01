function X = list_get(L,n);
%LIST_GET Get the value of list-element n
%       
%       X = list_get(L,n)
%         X - the list element is a list of n dimensional vectors (an array).
%         L - the list.
%         n - a pointer to a list element.
%
%       list_get returns the list-element indicated by the list
%       pointer n.  The list is a 'contourc' type list.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, February 12, 1997

if list_valid_pointer(L,n)
  X = L(:,n+1:n+L(size(L,1),n));
else
  error(sprintf('The pointer ''%d'' is not a valid pointer',n));
end
