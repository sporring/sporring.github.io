function J = ind_aos(I,iterations,stepsize,lambda,sigma);
%IND_AOS Isotropic Nonlinear Diffusion with Additive Operator Splitting (AOS)
%
%       J = ind_aos(I,iterations,stepsize,lambda,sigma);
%         J - the resulting double matrix
%         I - the double matrix original matrix
%         iterations - an optional iteration number
%             (0 <= iterations, default 1)
%         stepsize - an optional step size
%             (0 < stepsize, default 2)
%         lambda - the optional edge enhancement parameter
%             (0 < lambda, default 7)
%         sigma - the optional regularization parameter
%             (0 < sigma, default 1)
%
%       Copyright: Joachim Weickert & Jon Sporring, June 15, 1998
