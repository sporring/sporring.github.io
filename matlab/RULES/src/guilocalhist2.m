function guilocalhist2(action)
  
  global I;
  global minI;
  global maxI;
  global Ifig;
  global filename;
  global cursor;
  global sigmacircle;
  global alphacircle;
  global H;
  global Hfig;
  
  if nargin < 1
    clear all
    close all
    
    defaultSigma = 1;
    defaultAlpha = 5;
    defaultBeta = 3;
    defaultLevels = 255;
    
    backgroundcolor = [0.7 0.7 0.7];
    smallheight = 14;
    largeheight = 25;
    smallwidth = 152;
    largewidth = 186;
    halfdeltawidth = round((largewidth-smallwidth)/2);
    separation = 10;
    margin = 20;
    
    h0 = figure('Color',[0.8 0.8 0.8], ...
      'Colormap',gray(256), ...
      'Units','points', ...
      'Position',[878 174 largewidth+2*margin 5*largeheight+10*smallheight+25*separation+2*margin], ...
      'Tag','Fig1', ...
      'ToolBar','none');
    dx = margin;
    h1 = uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',backgroundcolor, ...
      'FontSize',12, ...
      'FontWeight','demi', ...
      'HorizontalAlignment','left', ...
      'ListboxTop',0, ...
      'Position',[dx margin largewidth smallheight], ...
      'String','Status: Ready!', ...
      'Style','text', ...
      'Tag','StaticText6');
    dx = dx+smallheight+3*separation;
    h1 = uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',backgroundcolor, ...
      'Callback','guilocalhist2 quit', ...
      'FontSize',12, ...
      'FontWeight','demi', ...
      'ListboxTop',0, ...
      'Position',[margin+halfdeltawidth dx smallwidth largeheight], ...
      'String','Quit', ...
      'Tag','Pushbutton3');
    dx = dx+largeheight+2*separation;
    h1 = uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',backgroundcolor, ...
      'Callback','guilocalhist2 start', ...
      'FontSize',12, ...
      'FontWeight','demi', ...
      'ListboxTop',0, ...
      'Position',[margin+halfdeltawidth dx smallwidth largeheight], ...
      'String','Evaluate', ...
      'Tag','Pushbutton1');
    dx = dx+largeheight+2*separation;
    h1 = uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',backgroundcolor, ...
      'Callback','guilocalhist2 simulate', ...
      'FontSize',12, ...
      'FontWeight','demi', ...
      'ListboxTop',0, ...
      'Position',[margin+halfdeltawidth dx smallwidth largeheight], ...
      'String','Simulate', ...
      'Tag','Pushbutton4');
    dx = dx+largeheight+2*separation;
    h1 = uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',backgroundcolor, ...
      'Callback','guilocalhist2 keep', ...
      'FontSize',12, ...
      'FontWeight','demi', ...
      'ListboxTop',0, ...
      'Position',[margin+halfdeltawidth dx smallwidth largeheight], ...
      'String','Keep Figure', ...
      'Style','checkbox', ...
      'Tag','Checkbox1', ...
      'Value',1);
    dx = dx+largeheight+2*separation;
    h1 = uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',backgroundcolor, ...
      'Callback','guilocalhist2 levels', ...
      'ListboxTop',0, ...
      'Max',256, ...
      'Min',1, ...
      'Position',[margin+halfdeltawidth dx smallwidth smallheight], ...
      'SliderStep',[0.0039 0.1], ...
      'Style','slider', ...
      'Tag','Slider4', ...
      'Value', 256);
    dx = dx+smallheight+separation;
    h1 = uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',backgroundcolor, ...
      'FontSize',12, ...
      'FontWeight','demi', ...
      'HorizontalAlignment','left', ...
      'ListboxTop',0, ...
      'Position',[margin+halfdeltawidth dx smallwidth smallheight], ...
      'String','Levels: 256', ...
      'Style','text', ...
      'Tag','StaticText7');
    dx = dx+smallheight+2*separation;
    h1 = uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',backgroundcolor, ...
      'Callback','guilocalhist2 beta', ...
      'ListboxTop',0, ...
      'Position',[margin+halfdeltawidth dx smallwidth smallheight], ...
      'SliderStep',[0.01 0.1], ...
      'Style','slider', ...
      'Tag','Slider3');
    dx = dx+smallheight+separation;
    h1 = uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',backgroundcolor, ...
      'FontSize',12, ...
      'FontWeight','demi', ...
      'HorizontalAlignment','left', ...
      'ListboxTop',0, ...
      'Position',[margin+halfdeltawidth dx smallwidth smallheight], ...
      'String','Beta: 0.0', ...
      'Style','text', ...
      'Tag','StaticText3');
    dx = dx+smallheight+2*separation;
    h1 = uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',backgroundcolor, ...
      'Callback','guilocalhist2 alpha', ...
      'ListboxTop',0, ...
      'Position',[margin+halfdeltawidth dx smallwidth smallheight], ...
      'SliderStep',[0.01 0.1], ...
      'Style','slider', ...
      'Tag','Slider2');
    dx = dx+smallheight+separation;
    h1 = uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',backgroundcolor, ...
      'FontSize',12, ...
      'FontWeight','demi', ...
      'HorizontalAlignment','left', ...
      'ListboxTop',0, ...
      'Position',[margin+halfdeltawidth dx smallwidth smallheight], ...
      'String','Alpha: 0.0', ...
      'Style','text', ...
      'Tag','StaticText2');
    dx = dx+smallheight+2*separation;
    h1 = uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',backgroundcolor, ...
      'Callback','guilocalhist2 sigma', ...
      'ListboxTop',0, ...
      'Position',[margin+halfdeltawidth dx smallwidth smallheight], ...
      'SliderStep',[0.01 0.1], ...
      'Style','slider', ...
      'Tag','Slider1');
    dx = dx+smallheight+separation;
    h1 = uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',backgroundcolor, ...
      'FontSize',12, ...
      'FontWeight','demi', ...
      'HorizontalAlignment','left', ...
      'ListboxTop',0, ...
      'Position',[margin+halfdeltawidth dx smallwidth smallheight], ...
      'String','Sigma: 0.0', ...
      'Style','text', ...
      'Tag','StaticText1');
    dx = dx+smallheight+2*separation;
    h1 = uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',backgroundcolor, ...
      'FontSize',12, ...
      'FontWeight','demi', ...
      'HorizontalAlignment','left', ...
      'ListboxTop',0, ...
      'Position',[margin+halfdeltawidth dx smallwidth smallheight], ...
      'String','File:', ...
      'Style','text', ...
      'Tag','StaticText5');
    dx = dx+smallheight+separation;
    h1 = uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',backgroundcolor, ...
      'Callback','guilocalhist2 read', ...
      'FontSize',12, ...
      'FontWeight','demi', ...
      'ListboxTop',0, ...
      'Position',[margin+halfdeltawidth dx smallwidth largeheight], ...
      'String','Load New Image', ...
      'Tag','Pushbutton2');
    dx = dx+largeheight+3*separation;
    h1 = uicontrol('Parent',h0, ...
      'Units','points', ...
      'BackgroundColor',backgroundcolor, ...
      'FontSize',18, ...
      'FontWeight','bold', ...
      'ListboxTop',0, ...
      'Position',[margin dx largewidth largeheight], ...
      'String','Local Histograms', ...
      'Style','text', ...
      'Tag','StaticText4');
    % if nargout > 0, fig = h0; end
  else
    switch(action)
      case 'read'
        set(findobj('Tag','StaticText6'),'String','Status: Ready!');
        delete(cursor);
        cursor = [];
        delete(sigmacircle);
        sigmacircle = [];
        delete(alphacircle);
        alphacircle = [];
        H = [];
        Hfig = [];
        
        [fname,pname] = uigetfile('*');
        if fname ~= 0
          filename = fname;
          set(findobj('Tag','StaticText5'),'String',['File: ',filename]);
          if exist([pname,fname],'file') > 0
            try,
              I = double(imread([pname,filename]));
            catch,
              set(findobj('Tag','StaticText6'),'String','Status: Could not Read image!'),
            end
            if ndims(I) ~= 2,
              set(findobj('Tag','StaticText6'),'String','Status: Not a Gray image!');
              I = [];
            end
          else
            set(findobj('Tag','StaticText6'),'String','Status: File doesn''t exist!');
            I = [];
          end
        end
        if ~isempty(I)
          if isempty(Ifig) | ~ishandle(Ifig)
            Ifig = figure;
          else
            figure(Ifig);
          end
          Ih = findobj('Tag','Image1');
          if ~isempty(Ih)
            delete(Ih);
          end
          Ih = image(I);
          axis image;
          colormap(gray(256));
          % set(findobj('Tag','Slider1'),'Max',max(size(I)));
          % set(findobj('Tag','Slider2'),'Max',max(size(I)));
          % set(findobj('Tag','Slider3'),'Max',max(size(I)));
          set(Ih,'Tag','Image1');
          set(Ih,'ButtonDownFcn','guilocalhist2 image');
          [cursor,sigmacircle,alphacircle] = updatecircle('r','g','b');
        end
      case 'sigma'
        set(findobj('Tag','StaticText6'),'String','Status: Ready!');
        H = [];
        sigma = T(get(gcbo,'Value'));
        h = findobj('Tag','StaticText1');
        set(h,'String',sprintf('Sigma: %.2f',sigma));
        [cursor,sigmacircle,alphacircle] = updatecircle('r','g','b');
      case 'alpha'
        set(findobj('Tag','StaticText6'),'String','Status: Ready!');
        H = [];
        alpha = T(get(gcbo,'Value'));
        h = findobj('Tag','StaticText2');
        set(h,'String',sprintf('Alpha: %.2f',alpha));
        [cursor,sigmacircle,alphacircle] = updatecircle('r','g','b');
      case 'beta'
        set(findobj('Tag','StaticText6'),'String','Status: Ready!');
        H = [];
        beta = T(get(gcbo,'Value'));
        h = findobj('Tag','StaticText3');
        set(h,'String',sprintf('Beta: %.2f',beta));
      case 'levels'
        set(findobj('Tag','StaticText6'),'String','Status: Ready!');
        H = [];
        levels = round(get(gcbo,'Value'));
        h = findobj('Tag','StaticText7');
        set(h,'String',sprintf('Levels: %d',levels));
      case 'keep'
        % Nothing to do.  Value will be read when needed.
      case 'simulate'
        set(findobj('Tag','StaticText6'),'String','Status: Ready!');
        if isempty(I) | all(I==0)
          set(findobj('Tag','StaticText6'),'String','Status: Load image!');
        else
          set(findobj('Tag','StaticText6'),'String','Status: Working...');
          drawnow;
          sigma = T(get(findobj('Tag','Slider1'),'Value'));
          alpha = T(get(findobj('Tag','Slider2'),'Value'));
          beta = T(get(findobj('Tag','Slider3'),'Value'));
          levels = round(get(findobj('Tag','Slider4'),'Value'));
          H = localhist2(I*(levels-1)/255,sigma,beta,alpha,levels);
          H = H + 0.0001*rand(size(H));
          C = cumsum(H,3);
          C = C./repmat(C(:,:,end),[1,1,size(H,3)]);
          y = rand(size(H,1),size(H,2));
          J = zeros(size(y));
          for si = 1:size(J,1)
            for sj = 1:size(J,2)
              J(si,sj) = interpn(squeeze(C(si,sj,:)),1:size(H,3),y(si,sj));
            end
          end
          figure
          imagesc(J); colormap(gray); axis image; colorbar;
          title(sprintf('%s: Sigma=%.2f, Alpha=%.2f, Beta=%.2f, Levels=%d',filename,sigma,alpha,beta,levels));
          set(findobj('Tag','StaticText6'),'String','Status: Ready!');
        end
      case 'start'
        set(findobj('Tag','StaticText6'),'String','Status: Ready!');
        if isempty(I) | all(I==0)
          set(findobj('Tag','StaticText6'),'String','Status: Load image!');
        else
          set(findobj('Tag','StaticText6'),'String','Status: Working...');
          drawnow;
          sigma = T(get(findobj('Tag','Slider1'),'Value'));
          alpha = T(get(findobj('Tag','Slider2'),'Value'));
          beta = T(get(findobj('Tag','Slider3'),'Value'));
          levels = round(get(findobj('Tag','Slider4'),'Value'));
          H = localhist2(I*(levels-1)/255,sigma,beta,alpha,levels);
          set(findobj('Tag','StaticText6'),'String','Status: Ready!');
          [cursor,sigmacircle,alphacircle] = updatecircle('r','g','b');
        end
      case 'quit'
        set(findobj('Tag','StaticText6'),'String','Status: Closing Figures.');
        drawnow;
        close all;
      case 'image'
        set(findobj('Tag','StaticText6'),'String','Status: Ready!');
        [cursor,sigmacircle,alphacircle] = updatecircle('r','g','b');
      otherwise
        warning('Unknown Action.  Doing nothing.');
    end
  end
end

function x = T(x)
  
  global I;
  
  if isempty(I)
    x = 0;
  elseif x == 1
    x = inf;
  elseif x > 0
    x = 0.1*exp((log(max(size(I)))-log(0.1))*x);
  end
end

function [cursor,sigmacircle,alphacircle] = updatecircle(xcolor,sigmacolor,alphacolor)
  
  global I;
  global filename;
  global cursor;
  global sigmacircle;
  global alphacircle;
  global H;
  global Hfig;
  
  delete(cursor);
  cursor = [];
  delete(sigmacircle);
  sigmacircle = [];
  delete(alphacircle);
  alphacircle = [];
  
  h = findobj('Tag','Image1');
  ha = get(h,'Parent');
  x = get(ha,'CurrentPoint');
  if ~isempty(x)
    x = round(x(1,1:2));
    sigma = T(get(findobj('Tag','Slider1'),'Value'));
    alpha = T(get(findobj('Tag','Slider2'),'Value'));
    
    [cursor,sigmacircle,alphacircle] = setcircle(ha,x,'r',sigma,'b',alpha,'g');
    
    if ~isempty(H)
      if any(x(1,1:2) < 1) | any(x(1,1:2) > size(I))
        set(findobj('Tag','StaticText6'),'String','Status: Choose a point.');
      else
        beta = T(get(findobj('Tag','Slider3'),'Value'));
        levels = round(get(findobj('Tag','Slider4'),'Value'));
        keep = round(get(findobj('Tag','Checkbox1'),'Value'));
        if keep | isempty(Hfig) | ~ishandle(Hfig)
          Hfig = figure;
        else
          figure(Hfig);
        end
        minI = min(min(I));
        maxI = max(max(I));
        plot(linspace(minI,maxI,size(H,3)),squeeze(H(x(2),x(1),:)));
        set(gca,'xlim',[minI,maxI]);
        title(sprintf('%s: Sigma=%.2f, Alpha=%.2f, Beta=%.2f, Levels=%d, X=(%d,%d)',filename,sigma,alpha,beta,levels,x(1),x(2)));
        xlabel('Intensity');
        drawnow;
      end
    end
  end
end

function [cursor,sigmacircle,alphacircle] = setcircle(handle,x,xcolor,sigma,sigmacolor,alpha,alphacolor)
  
  set(handle,'NextPlot','add');
  axes(handle);
  tmp = [2:4];
  cursor = plot(tmp+x(1),ones(size(tmp))*x(2),xcolor, ...
    x(1)-tmp,ones(size(tmp))*x(2),xcolor, ...
    ones(size(tmp))*x(1),tmp+x(2),xcolor, ...
    ones(size(tmp))*x(1),x(2)-tmp,xcolor);
  theta = linspace(0,2*pi,60);
  sigmacircle = plot(x(1)+2*sigma*cos(theta),x(2)+2*sigma*sin(theta),sigmacolor);
  alphacircle = plot(x(1)+2*alpha*cos(theta),x(2)+2*alpha*sin(theta),alphacolor);
end

function H = localhist2(I,sigma,beta,alpha,N)
  
  % Checking parameters and setting default;
  if nargin < 1
    error('Image I must be supplied! Usage: H = localhist2(I,sigma,beta,alpha,N)');
  end
  if nargin < 2
    sigma = 1;
  end
  if nargin < 3
    beta = 1;
  end
  if nargin < 4
    alpha = 1;
  end
  if nargin < 5
    if isreal(I)
      N = ceil(max(max(I)));
    else
      warning(['You will get faster computation if you specify N, when Image is ' ...
        'given as its Fourier Transform.']);
      N = ceil(max(max(real(ifft2(I)))));
    end
  end
  
  % Smoothing image
  if sigma > 0
    if isreal(I)
      L = real(ifft2(scale2(fft2(I),sigma,0,0)));
    else
      L = real(ifft2(scale2(I,sigma,0,0)));
    end
  else
    if isreal(I)
      L = I;
    else
      L = real(ifft2(I));
    end
  end
  
  % Calculating soft isophote and sum under window
  H = zeros([size(I),N]);
  for i = 1:N
    % disp(i);
    % v = [1:N]-i;
    % MAP = exp(-v.^2/(eps+2*beta^2));
    % figure(1); plot(MAP); pause
    % J = ind2gray(L+1,[MAP' MAP' MAP']);
    % figure(1); plot(J(1,:)); pause
    J = exp(-(L-i).^2/(eps+2*beta^2));
    
    if alpha > 0
      H(:,:,i) = real(ifft2(scale2(fft2(J),alpha,0,0)));
      % figure(1); plot(H(1,:,i+1)); pause
    else
      H(:,:,i) = J;
    end
  end
  
end

function Is = scale2(I,s,dr,dc);
  %SCALE2  Gaussian Scale-Space using the Fourier Domain
  %
  %       FL = scale(I,s,dr,dc)
  %         I - the Fourier transform of a matrix or an image
  %         s - the standard deviation of the Gaussian (must be larger than 0)
  %             (default 1)
  %         dr - the derivative order in the direction of the rows
  %              (default 0)
  %         dc - the derivative order in the direction of the columns
  %              (default 0)
  %
  %       This is an implementation of the Gaussian scale space on matrices.
  %       The convolution is implemented in the Fourier domain and for that
  %       reason the number of rows and columns of the matrix must be powers
  %       of 2.
  %       Fractional valued dr and dc are possible, but be warned the result
  %       will probably be complex.
  %       The complexity of this algorithm is O(n) where n is the total number
  %       of elements of I.
  %
  %       To calculate an image of scale (variance) 2^2 use,
  %         L = real(ifft2(scale(fft2(I),2,0,0)));
  %       To derive an image once in the direction of rows at scale 1^2 do,
  %         L = real(ifft2(scale(fft2(I),1,1,0)));
  %
  %       Copyright: Jon Sporring, January 1, 1996
  
  if nargin < 4
    dc = 0;
  end
  if nargin < 3
    dr = 0;
  end
  if nargin < 2
    s = 1;
  end
  if nargin < 1
    error('No input!  Use: Is = scale(FI,s,dr,dc);');
  end
  
  if any([s,dr,dc] < 0)
    error('s, dr, and dc must be larger than zero');
  else
    if s == 0
      if (dr == 0) & (dc == 0)
        Is = I;
      else
        Is = zeros(size(I));
      end
    else
      rows = size(I,1);
      cols = size(I,2);
      
      G = zeros(rows,cols);
      DG = ones(rows,cols);
      if s == Inf
        G(1,1) = 1;
        Is = I.*G;
      else
        % Calculate the Fourier transform of a gaussian fct. and the
        % Differentiation matrix.
        if (rows > 1) & (cols > 1)
          if any(rem(size(I),2) ~= 0)
            error('Each image sides must be divisible by 2');
          else
            % 2 dimensional image
            x = [0:rows/2]'/(rows-1);
            y = [0:cols/2]/(cols-1);
            
            G(1:rows/2+1, 1:cols/2+1) = exp(-(repmat(((s*2*pi)*x).^2/2,[1,cols/2+1])+repmat(((s*2*pi)*y).^2/2,[rows/2+1,1])));
            G(rows/2+1:rows, 1:cols/2+1) = flipud(G(2:rows/2+1, 1:cols/2+1));
            G(:, cols/2+1:cols) = fliplr(G(:, 2:cols/2+1));
            
            if any([dr,dc] > 0)
              DG(1:rows/2+1, 1:cols/2+1) = repmat(x.^dr,[1,cols/2+1]).*repmat(y.^dc,[rows/2+1,1])*(sqrt(-1)*2*pi)^(dr+dc);
              DG(rows/2+1:rows, 1:cols/2+1) = (-1)^dr*flipud(DG(2:rows/2+1, 1:cols/2+1));
              DG(:,cols/2+1:cols) = (-1)^dc*fliplr(DG(:, 2:cols/2+1));
            end
          end
        else
          if rem(length(I),2) ~= 0
            error('Number of vector elements must be divisible by 2');
          else
            % 1 dimensional image
            [val,ind] = max([rows,cols]);
            x = [0:val/2]'/(val-1);
            G(1:val/2+1) = exp(-([0:val/2]'/(val-1)).^2*(s*2*pi)^2/2);
            G(val/2+1:val) = flipdim(G(2:val/2+1),ind);
            
            DG = (j*2*pi*x').^dr;
            DG(val/2+1:val) = flipdim(DG(2:val/2+1),ind);
          end
        end
        
        Is = I.*G.*DG;
      end
    end
  end
end