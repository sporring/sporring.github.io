function x = chebyshevspace(a,b,n)
%CHEBYSHEVSPACE Calculates the Chebyshevs sampling points.
%       
%       x = chebyshevspace(a,b,n)
%         x - A vector of Chebyshev sampling points
%         a - The left boundary of the interval
%         b - The right boundary of the interval
%         n - The total number of sampling points
%
%       This function calculates the Chebyshev sampling points known to be
%       near optimal when fitting a polynomial to a random curve.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, January 15, 1997


x = (a+b - (a-b)*cos((2*[1:n]-1)*pi/(2*n)))/2;

