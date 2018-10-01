function C = isocurvature(Lx, Ly, Lxx, Lxy, Lyy);
%ISOCURVATURE Calculate the isophote curvature from derivative images
%       
%       C = isocurvature(Lx, Ly, Lxx, Lxy, Lyy)
%         C - the curvature image
%         Lx, Ly - the first order derivatives in x and y directions
%         Lxx, Lxy, Lyy - the second derivatives
%
%       Isocurvature calculates the isophote curvature using
%       (2*Lx.*Ly.*Lxy-Lyy.*(Lx.^2)-Lxx.*(Ly.^2))./((Lx.^2+Ly.^2).^1.5).
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, January 13, 1997

C = (Lx.*(2*Ly.*Lxy-Lyy.*Lx)-Lxx.*(Ly.^2))./((Lx.^2+Ly.^2).^1.5);
