function b = nlfilter3(varargin)
%NLFILTER3 General sliding-neighborhood operations in 3 dimensions.
%   B = NLFILTER(A,[M N P],FUN) applies the function FUN to each
%   M-by-N-by-O sliding block of the grayscale image A. FUN is a function
%   that accepts an M-by-N-by-P matrix as input and returns a scalar:
%
%       C = FUN(X)
%
%   FUN must be a FUNCTION_HANDLE.
%
%   C is the output value for the center pixel in the M-by-N-by-P block X.
%   NLFILTER calls FUN for each pixel in A. NLFILTER zero pads the
%   M-by-N-by-P block at the edges, if necessary.
%
%   Jon Sporring, 2015/11/15, University of Copenhagen

[a, nhood, fun, params, padval] = parse_inputs(varargin{:});

% Expand A
[ma,na,oa] = size(a);
aa = repmat(padval, size(a)+nhood-1);
aa(floor((nhood(1)-1)/2)+(1:ma),floor((nhood(2)-1)/2)+(1:na),floor((nhood(3)-1)/2)+(1:oa)) = a;

% Find out what output type to make.
rows = 0:(nhood(1)-1);
cols = 0:(nhood(2)-1);
depths = 0:(nhood(3)-1);
b = zeros(size(a));

% create a waitbar if we are able
if images.internal.isFigureAvailable()
    wait_bar = waitbar(0,'Applying neighborhood operation...');
else
    wait_bar = [];
end

% Apply fun to each neighborhood of a
for i=1:ma
    for j=1:na
        for k=1:oa
            x = aa(i+rows,j+cols,k+depths);
            b(i,j,k) = feval(fun,x,params{:});
            
            % udpate waitbar
            if ~isempty(wait_bar)
                waitbar(i/ma,wait_bar);
            end
        end
    end
end

close(wait_bar);



%%%
%%% Function parse_inputs
%%%
function [a, nhood, fun, params, padval] = parse_inputs(varargin)

switch nargin
    case {0,1,2}
        error(message('images:nlfilter3:tooFewInputs'))
    otherwise
        % NLFILTER(A, [M N O], 'fun')
        a = varargin{1};
        nhood = varargin{2};
        fun = varargin{3};
        params = cell(0,0);
        padval = 0;
end
