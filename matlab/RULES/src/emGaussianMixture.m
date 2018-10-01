function [mu,C,Posteriori,Likelihood,Prior,l] = em(I,mu,C,Prior);
% EM - An expectation maximization algorithm for a Gaussian Mixture Model
%
% [mu,C,Posteriori,Likelihood,Prior,l] = em(I,mu,C,Prior)
%    I           - is an NxD data set, with N data points of dimension D
%    mu          - is an KxD set of mean vectors
%    C           - is a KxD^2 set of linearized covariance matrices
%    Prior       - is a K vector of prior probabilities
%    Posteriori  - is a NxK vector of posteriori probabilities
%    Likelihood  - is a NxK vector of likelihood probabilities
%    l           - is an estimated codelength
%
% The algorithm is an extension of vector quantification, which softly
% fitting each data point I:NxD into a mixture of Gaussians consisting of K
% Gaussians represented by a mean mu(i) and a covariance matrix
% reshape(C(i,:),[D,D]).
%
%                                Jon Sporring, DIKU, 17/11-2003.
  
  [Posteriori,Likelihood] = expectation(I,mu,C,Prior);
  [mu,C,Prior] = maximazation(I,mu,C,Posteriori);
  % l = nml(I,mu,C,Posteriori,Likelihood,Prior);
  l = mdl(I,mu,C,Posteriori,Likelihood,Prior);
  
function [Posteriori,Likelihood,Prior] = expectation(I,mu,C,Prior)
% Expectation-step

  Likelihood = zeros(size(I,1),size(mu,1));
  for i = 1:size(mu,1)
    % First we calculate the deviation from the mean in class i
    DI = I-repmat(mu(i,:),[size(I,1),1]);
    % Then we fetch the i'th covariance matrix
    c = reshape(C(i,:),[size(I,2),size(I,2)]);
    
    % To calculate the likelihood, we want to do the following:
    %   Likelihood = zeros(size(Likelihood,1),1);
    %   for j = 1:size(DI,1)
    %     Likelihood(j) = Likelihood(j)+DI(j,:)*inv(c)*DI(j,:)';
    %   end
    %   Likelihood = eps+exp(-Likelihood/2);
    % which in Matlab is fastest as:
    Likelihood(:,i) = eps+exp(-sum(DI.*(inv(c)*DI')',2)/2)/(2*pi)^(size(c,1)/2)/det(c)^(1/2);
    % Likelihood(:,i) = Likelihood(:,i)/sum(Likelihood(:,i));
    
    % The unnormalized posteriori is found by Bayes formula:
    Posteriori(:,i) = Prior(i)*Likelihood(:,i);
  end
  % The posteriori should be normalized wrt. the class
  Posteriori = Posteriori./repmat(sum(Posteriori,2),[1,size(Posteriori,2)]);

function [mu,C,Prior] = maximazation(I,mu,C,Posteriori)
% Maximization-step

  for i = 1:size(Posteriori,2)
    % New mean values are estimated as the marginal of the joint distribution
    % assuming uniform probability of position:
    mu(i,:) = sum(I.*repmat(Posteriori(:,i),[1,size(I,2)]))/sum(Posteriori(:,i));

    % New covariance values are also estimated as the marginal of the joint
    % distribution assuming uniform probability of position:
     DI = I-repmat(mu(i,:),[size(I,1),1]);
     c = zeros(size(I,2),size(I,2));
     for s = 1:size(I,2)
       for t = s:size(I,2)
         c(s,t) = sum(Posteriori(:,i).*DI(:,s).*DI(:,t))/sum(Posteriori(:,i)); 
         if t > s
           c(t,s) = c(s,t);
         end
       end
     end
     % In case the covariance matrix becomes very small, then the Gaussian model
     % becomes singular.  Hence, we regularized by adding a small
     % diagonal...
     c = 0.1*eye(size(c))+c;
     C(i,:) = c(:)';
  end
  
  % Finally, the prior is re-estimated as the marginal distribution under
  % uniform position:
  Prior = sum(Posteriori,1)/size(I,1);

function l = nml(I,mu,C,Posteriori,Likelihood,Prior)
% Normalized Maximum Likelihood
  
  l = sum(-log(sum(repmat(Prior,[size(Likelihood,1),1]).*Likelihood,2)))/size(Likelihood,1)/size(I,2);

function l = mdl(I,mu,C,Posteriori,Likelihood,Prior)
% Minimum Description Length (Likelihood x Prior code)

  [val,ind] = max(Posteriori,[],2);
  i = sub2ind(size(Likelihood),1:size(Likelihood,1),ind');
  l = (sum(-log(Likelihood(i).*Prior(ind))))/size(Likelihood,1)/size(I,2);
