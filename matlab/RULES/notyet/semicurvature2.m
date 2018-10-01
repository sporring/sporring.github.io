function [cs,ss] = semicurvature2(C)
%SEMICURVATURE2 Calculates the semi-curvature function in the Frenet frame.
%
%       [cs,ss] = semicurvature2(C)
%         cs - the resulting list of curvature functions
%         ss - the resulting list of step functions
%         C - a list of lists of points on the isophote as given by contourc
%
%       Semicurvature2 calculates the angular difference between two
%       consecutive displacement-vectors divided by the step-length.
%       This can be seen to be an approximation of the curvature
%       function in the Frenet frame where the first frame vector is the
%       tangent and the second the normal in a right handed coordinate
%       system.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, February 5, 1997

cs = [];
ss = [];
a = list_first(C);
while a > 0
  % CONTOUR USES X,Y INSTEAD OF ROW,COL
  V = list_get(C,a);

  % CLOSE GAPPING LOOPS
  if max(abs(V(:,1)-V(:,size(V,2)))) ~= 0
    V = [V, V(:,1)];
  end

  [c,dV,s] = semicurvature(V);

  cs = list_insert(cs,0,c,[]);
  cs = list_insert(cs,0,V(:,size(V,2))'-V(:,size(V,2)-1)',[]);
  cs = list_insert(cs,0,V(:,1)',[]);
  ss = list_insert(ss,0,s,[]);

  a=list_next(C,a);
end
