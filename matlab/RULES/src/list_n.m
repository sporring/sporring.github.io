function n_L = list_n(L,n);
%LIST_N Get the pointer to the n'th elment of the list
%
%       n_L = list_n(L,n)
%         n_L - a pointer to a list element
%         L - the list.
%         n - an integer in the range 1..list_size(L).
%
%       list_n returns a pointer to the the n'th element of the list.
%       If n < 1, then the pointer to the first element is returned,
%       if n > list_size(L), then the pointer to the last element is
%       returned. 
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, April 2, 1997

n_L = list_first(L);
if n < 1
  disp('Warning (list_n): Element number less than 1');
  dbstack
else
  while (n_L > 0) & (n > 1)
    n_L = list_next(L, n_L);
    n = n - 1;
  end
  if n_L == 0
    disp('Warning (list_n): Element number too large');
    dbstack
    n_L = list_last(L);
  end
end
