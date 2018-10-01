function L = list_delete(L,n);
%LIST_DELETE Removes element pointed to by n
%       
%       L = list_delete(L,n)
%         L - the list.
%         n - a pointer to a list element.
%
%       list_delete removes the element pointed to by n.  If n is
%       before the first element or after the last, the first or last
%       is removed respectively.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, February 12, 1997

if n < 1
  n = list_first(L);
else
  if n > size(L,2)
    n = list_last(L);
  end
end

if list_valid_pointer(L,n)
  L = [L(:,1:n-1), L(:,n+1+L(size(L,1),n):size(L,2))];
else
  error(sprintf('The pointer ''%d'' is not a valid pointer',n));
end
