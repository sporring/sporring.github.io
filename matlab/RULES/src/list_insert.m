function L = list_insert(L,n,e,a);
%LIST_INSERT Inserts an element into list at position n.
%       
%       L = list_insert(L,n,e,a)
%         L - the list (after and before).
%         n - a pointer to a list element.
%         e - an element to be inserted.
%         a - the size(e,1)-1 attributes.
%
%       list_insert inserts the element e at position pointed to by n.
%       If n is before the first element or after the last, the
%       element is inserted first or last respectively.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, February 12, 1997

% Comment:  The particular list form chosen for contourc structures
% have a problem with empty elements.  I.e. if I have a list with row
% size n, one cannot make an empty element with equal number of rows.
% (Note, v.5 has the possibility with (full(sparse(n,0))).)  There are
% two concerns to address.  First, an already full list of abitrary
% numbers of rows should be able to hold an empty element, and
% secondly an empty list or one with only empty elements should allow
% for a 'late' definition of the number of rows.  The first is easily
% implemented, while the last is done by defering the decision to the
% first non-empty element in the list, and post-correcting the row
% size of the list.

% Reset point if pointing out of list
if n < 1
  n = 1;
else
  if n > size(L,2)
    n = size(L,2)+1;
  end
end

% Is this a list of only empty elements, and e is the first non-empty
if (size(e,1) > size(L,1)) & ~isempty(L)
  empty = sum(L(size(L,1),:));
  if empty > 0
    error('Trying to insert an element with more rows than the list.')
  else
    L = zeros(size(e,1),size(L,2));
  end
end

attributes = zeros(max([size(e,1),size(L,1)])-1,1);
if ~isempty(a)
  if size(a,1) == 1
    a = a';
  end
  cut = min(size(a,1),size(attributes,1));
  attributes(1:cut,1) = a(1:cut);
end
if isempty(e)
  L = [L(:,1:n-1) [attributes;size(e,2)] L(:,n:size(L,2))];
else
  L = [L(:,1:n-1) [attributes;size(e,2)] e L(:,n:size(L,2))];
end
