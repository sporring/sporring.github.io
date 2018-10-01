function L = list_set_attributes(L,n,a);
%LIST_SET_ATTRIBUTES Set the attribute values of list-element n
%       
%       list_set_attributes(L,n,a)
%         L - the list.
%         n - a pointer to a list element.
%         a - the new attributes of the list element.
%
%       list_set_attributes replaces the current attribute values with
%       a for element n.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, February 12, 1997

if list_valid_pointer(L,n)
  L(1:size(L,1)-1,n) = a;
else
  error(sprintf('The pointer ''%d'' is not a valid pointer',n));
end
