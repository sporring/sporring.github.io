function c = semicurvature3(V)
%SEMICURVATURE3 Calculates the semi-curvature function in the Frenet frame.
%
%       c = semicurvature3(V)
%         c - the resulting curvature function
%         V - a list of points on the isophote as given by isophote
%
%       Semicurvature3 calculates the angular difference between two
%       consecutive displacement-vectors divided by the step-length.
%       This can be seen to be an approximation of the curvature
%       function in the Frenet frame where the first frame vector is the
%       tangent and the second the normal in a right handed coordinate
%       system.
%
%       Copyright: IBM Almaden Reserach Center & Jon Sporring, February 5, 1997

c = zeros(1,size(V,2)-1);
V = [V(:,size(V,2)-3:size(V,2)-2),V,V(:,2)];
px = zeros(size(c,2),3);
py = zeros(size(c,2),3);

for i=1:size(V,2)-4
  px(i,:) = polyfit([-2:2],V(1,i:i+4),2);
  py(i,:) = polyfit([-2:2],V(2,i:i+4),2);
end
c = (2*(py(:,2).*px(:,1)-px(:,2).*py(:,1))./(px(:,2).^2+py(:,2).^2).^1.5)';
