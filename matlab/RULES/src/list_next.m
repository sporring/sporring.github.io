function m = list_next(L,n);
%LIST_NEXT Returns a pointer to the next element in L
%       
%       m = list_next(L,n)
%         m - a pointer to the next list element.
%         L - the list.
%         n - a pointer to the current list element.
%
%       list_next finds the pointer to the next element in the
%       list.  If n is the last element, 0 is returned.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, February 12, 1997

m = n+1+L(size(L,1),n);
if m > size(L,2)
  m = 0;
end
