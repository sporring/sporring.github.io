function [c,dV,s] = semicurvature(V)
%SEMICURVATURE Calculates the semi-curvature function in the Frenet frame.
%
%       [c,dV,s] = semicurvature(V)
%         c - the resulting curvature function
%         dV - the resulting tangent function
%         s - the resulting tangent length function
%         V - a list of points on the isophote as given by isophote
%
%       Semicurvature calculates the angular difference between two
%       consecutive displacement-vectors divided by the step-length.
%       This can be seen to be an approximation of the curvature
%       function in the Frenet frame where the first frame vector is the
%       tangent and the second the normal in a right handed coordinate
%       system.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, February 5, 1997

c = zeros(1,size(V,2)-1);
dV = [V(:,1)-V(:,size(V,2)-1),V(:,2:size(V,2))-V(:,1:size(V,2)-1)];
s = sqrt(sum(dV.^2));

c = acos(dot(dV(:,1:size(dV,2)-1),dV(:,2:size(dV,2)))./(s(1:size(dV,2)-1).*s(2:size(dV,2))))./s(2:size(dV,2));

for i=1:size(V,2)-1
  c(i) = sign(det([dV(:,i) dV(:,i+1)]))*c(i);
end
dV = dV(2:size(dV,2));
s = s(2:size(s,2));
