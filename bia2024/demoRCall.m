% This script gives an example on how RCall can be used to work with R and the spatstat package 
% from matlab. The following has been tested on an m1 mac. You may need to customize for your setup
% 
% Jon Sporring
% Department of Computer Science
% University of Copenhagen
% August 11, 2024
%
% 1. Install R, which to my experience works best directly from  https://cran.r-project.org/
% then start R and install some packes:
%   install.packages("spatstat")
%   install.packages("lazyeval")
%   install.packages("GET")
% In R you can get access to lots of manuals, e.g.,
%   library("spatstat")
%   ?ppp
% and there are also more general descriptions
%   library("GET")
%   vignette("pointpatterns")
%
% 2. Install RCall https://www.mathworks.com/matlabcentral/fileexchange/104945-rcall
% 
% 3. I'm a fan of export_fig, so the last calls prints the figure to pdf with that function.
% delet that part, if you don't need it.
%
% RCall is a minimalistic package which works by 
% a) storing values in a mat file with Rpush possibly replacing earlier values,
% b) sequentially records R-commands in a global environment file with Rrun, and
% c) when Rpull is called, then R is started, the values are read, the commands are executed, 
%    and the command list is cleared (but not the values)
  
% In the following performs repeated experiements of drawing random point patterns and plotting
% Ripley's K-function, L-function, or the normalized L-function.
M = 10; % Number of repeated trials
N = 100; % number of spatial sample points
W = 500; % width of window
Ktype = 3; % K-type function id. 1->Kest(r), 2->Lest(r), 3->Lest(r)-r
q = linspace(0,W,100); % sample points of K-like functions
K = 999; % number of simulations for envelope test

% setup default plot values
set(groot,'defaultLineLineWidth',2) 
set(0,'DefaultaxesLineWidth', 1.5) 
set(0,'DefaultaxesFontSize', 18) 
set(0,'DefaultaxesFontWeight', 'bold') 
set(0,'DefaultTextFontSize', 18) 
set(0,'DefaultaxesFontName', 'Times new Roman') 
set(0,'DefaultlegendFontName', 'Times new Roman')

% setup Rcall
libpath = '/opt/homebrew/lib/R/library';
path =  '/opt/homebrew/bin/R';
Rclear
Rinit([],path,libpath)
Rpush('M',M,'N',N,'W',W,'q',q,'K',K)

% repeated experiments and plotting
wX = [0,W];
wY = [0,W];
Rpush('wX',wX,'wY',wY)
figure(1)
clf
hold on;
for i = 1:10
    x = wX(1) + (wX(2)-wX(1)) * rand(N,1);
    y = wY(1) + (wY(2)-wY(1)) * rand(N,1);
    Rpush('x',x,'y',y)
    Rrun('library("spatstat")')
    Rrun('p <- ppp(x,y,window=owin(c(wX),c(wY)))')
    if Ktype == 1
        Rrun('k <- as.function(Kest(p))')
    else
        Rrun('k <- as.function(Lest(p))')
    end
    Rrun('res <- k(q)')
    res = Rpull('res');
    if Ktype == 3
        res = res - q';
    end
    plot(q,res,'k-')
end
idx = find(~isnan(res),1,'last');
q = q(1:idx);
switch Ktype
    case 1, fq = pi*q.^2;
    case 2, fq = q;
    otherwise, fq = zeros(size(q));
end
plot(q,fq,'k-','LineWidth',3)
hold off
axis tight
xlabel('r')
switch Ktype
    case 1
        ylabel('Kest(r)')
    case 2
        ylabel('Lest(r)')
    otherwise
        ylabel('Lest(r) - r')
end
set(gcf,"color",'w')
export_fig("Ktype.pdf")

% envelope test for Poisson process
Rrun('library("spatstat")')
Rrun('library("GET")')
Rrun('p <- ppp(x,y,window=owin(c(wX),c(wY)))')
Rrun('env <- envelope(p, nsim=K, savefuns = TRUE, fun = Lest, simulate = expression(runifpoint(ex = p)))')
Rrun('test <- global_envelope_test(env)')
Rrun('res <- list(attr(test,"p"),test$lo,test$hi,test$r,test$obs,test$central)')
res = Rpull('res');
p = res.sub1;
lo = res.sub2;
hi = res.sub3;
r = res.sub4;
obs = res.sub5;
central = res.sub6;
figure(2)
fill([r; flipud(r)], [lo; flipud(hi)], 0.75*ones(1,3), EdgeColor='none');
hold on; plot(r,[lo,hi], 'k--', 'LineWidth', 2); hold off;
hold on; plot(r,[central,obs], 'k-', 'LineWidth', 2); hold off
axis tight
title(sprintf('p-value = %f',p))
xlabel('r')
if Ktype == 1
    ylabel('Kest(r)')
else
    ylabel('Lest(r)')
end
set(gcf,"color",'w')
export_fig("envelope.pdf")
