function n = list_first(L);
%LIST_FIRST Returns a pointer to the first element in L
%       
%       n = list_first(L)
%         n - a pointer to the first list element.
%         L - the list.
%
%       list_first returns a pointer to the first list element in the
%       list.  If the list is empty, 0 is returned.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, February 12, 1997

n = (size(L,2) > 0);
