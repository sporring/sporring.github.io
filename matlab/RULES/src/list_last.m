function n = list_last(L);
%LIST_LAST Returns a pointer to the first element in L
%       
%       n = list_last(L)
%         n - a pointer to the last list element.
%         L - the list.
%
%       list_last returns a pointer to the last list element in the
%       list.  If the list is empty, 0 is returned.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, February 12, 1997

n = list_first(L);
if n > 0
  while n ~= 0
    previous = n;
    n = list_next(L,n);
  end
  n = previous;
end

