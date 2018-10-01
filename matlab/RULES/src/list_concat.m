function newlist = list_concat(oldlist)
%LIST_CONCAT Concatenates all element of a list into one.
%
%       newlist = list_concat(oldlist)
%         newlist - the matrix of concatted elements
%         oldlist - a list of contourc type
%
%       list_concat returns a matrix of concatted elements from oldlist.
%       this can e.g. be used together with list_forall for easy evaluation
%       of operations etc.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, March 1, 1997

%       Note, a possible speed improvement: preallocation of newlist.

newlist = [];
n = list_first(oldlist);

while n > 0
  newlist = [newlist, list_get(oldlist,n)];
  n = list_next(oldlist,n);
end
