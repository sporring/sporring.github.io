function R = list_forall(L,string,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)
%LIST_FORALL executes a specified command for all elements of the list.
%
%       R = list_forall(L,string,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)
%         R - a list of return values.
%         L - the list whose elements are to be processed
%         string - the execution string
%         p1..p10 - the 10 optional parameters
%
%       This function evaluates the string by `eval' on each element
%       on the list.  The elements are called 'm', their attributes 'a',
%       the list index is 'i', the parameters are called 'p1'...'p10',
%       and the return value is assumed to be 'r'.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, April 7, 1997

R = [];

current = list_first(L);
i = 1;
while(current > 0)
  a = list_get_attributes(L,current);
  m = list_get(L,current);
  eval(string);
  if exist('r')
    R = list_insert(R,Inf,r,[]);
  end
  current = list_next(L,current);
  i = i+1;
end
