function new = remove_space(old)
%REMOVE_SPACE removes all extra spaces in a string
%
%       new = remove_space(old)
%         new - The new string with no double or more spaces
%         old - The old string
%
%       For some output operations it is usefull to have extra spaces
%       removed.  A single space in the beginning or end is not removed.

new = old;
j = 0;
space = 0;
for i = 1:size(old,2)
  if old(i) == ' '
    space = space + 1;
  else
    space = 0;
  end
  
  if space < 2
    j=j+1;
    new(j) = old(i);
  end
end

new = new(1:j);
