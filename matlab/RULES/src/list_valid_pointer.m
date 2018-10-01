function b = list_valid_pointer(L,n);
%LIST_VALID_POINTER Checks the validity of the pointer n in L
%       
%       b = list_valid_pointer(L,n)
%         b - 0 or 1 depending on the validity of n in L
%         L - the list.
%         n - a possibly pointer to a list element.
%
%       list_valid_pointer returns non-zero if the pointer is not a
%       valid list-element pointer
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, February 12, 1997

% A much more correct scheme could be implemented!.
b = (n <= size(L,2));
