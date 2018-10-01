function X = list_getn(L,n);
%LIST_GETN Get the n'th elment of the list
%       
%       X = list_getn(L,n)
%         X - the list element is a list of n dimensional vectors (an array).
%         L - the list.
%         n - an integer in the range 1..list_size(L).
%
%       list_getn returns the n'th element of the list.  The list is a
%       'contourc' type list.  If n < 1, then the first element is
%       return, if n > list_size(L), then the last is returned.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, April 2, 1997

X = list_get(L, list_n(L,n));
