function zc = crossings(A,B)
%CROSSINGS finds the points where both A and B are zero
%
%       zc = crossings(A,B);
%         zc - a list of fractional indices where A==B==0
%         A,B - 2D real images.
%
%       crossing gives the result in the x-y coordinate system as dictated
%       by contourc.
%
%       Function has been greatly simplified by using zerocrossing.  
%       It now uses the x-y coordinate system as all other image functions.  July 4, 2014.
%
%       July 2014: Improved detection accuracy and speed by use of
%       interpolation.
%
%       Copyright: Jon Sporring, January 1, 1996

C = list_toCell(contourc(A,[0,0]));
zc = [];
for i=1:length(C)
    c = C{i};
    if(sum((c(:,1)-c(:,end)).^2)<1)
        % we add end point if this is a closed contour
        c(:,end+1)=c(:,1);
    end
    val = interp2(B, c(1,:), c(2,:));
    ind = zerocrossing(val);
    zc = [zc,[interp1(c(1,:),ind);interp1(c(2,:),ind)]];
end
