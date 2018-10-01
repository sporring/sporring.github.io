function a = list_get_attributes(L,n);
%LIST_GET_ATTRIBUTES Get the attribute values of list-element n
%       
%       a = list_get_attributes(L,n)
%         a - the attributes of the list element.
%         L - the list.
%         n - a pointer to a list element.
%
%       list_get_attributes returns the element attributes indicated
%       by the list  pointer n.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, February 12, 1997

if list_valid_pointer(L,n)
  a = L(1:size(L,1)-1,n);
else
  error(sprintf('The pointer ''%d'' is not a valid pointer',n));
end
