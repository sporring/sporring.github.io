function dexpr = dexpr(E,I,sigma)
% DEXPR Parses and calculates differential expressions of images 
%       using gaussian derivatives.
%
%       Idexpr = dexpr(E,I,s)
%         E - differential expression in L and spatial parameters x and y. (string)
%         I - input image. (double matrix)
%         s - the scale (standard deviation) of the gaussian derivatives. (double)
%
%       E must comply to the EBNF below.       
%
%       Examples:
%         Idexpr = dexpr('L',I,s); 
%            Image I at scale s.
%         
%         Idexpr = dexpr('(Lx^2 + Ly^2)^(1/2)',I,s); 
%             Gradient magnitude image of I at scale s.
%  
%       EBNF for E:
%         E  = T {["+"|"-" E}
%         T  = [PE|F] {["*"|"/"|"^"] T}
%         PE = "(" E ")"
%         F  = "REAL" | D
%         D  = "L" {["x"|"y"]}
%
%       Copyright: Martin Lillholm, March 13, 2001. 
%       Actual calculations of image derivative from Jon Sporrings scale.

  global tokens;
  global iTokens;
  global token;

  tokens  = lex(E);
  iTokens = 1;
  nextToken;
  dexpr   = readE(double(I),sigma);
  if ~in(token.type, ['EOF'])
    msg = [token.lexeme ' cannot succeed ' tokens(iTokens-2).lexeme '.'];
    error(msg);
  end
  
function E=readE(I,sigma)
  global token;

  E = readT(I,sigma);

  if in(token.type, ['add'; 'sub'])
    tmpToken = token;
    nextToken;
    tmp = readE(I,sigma);
    switch tmpToken.type
     case 'add'
      E = E + tmp;
     case 'sub'
      E = E - tmp;
    end
  end

function T = readT(I,sigma)
  global token;

  if in(token.type, ['lpa'])
    nextToken;
    T = readPE(I,sigma);
  else
    T = readF(I,sigma);
  end

  if in(token.type, ['mul'; 'div';'pow'])
    tmpToken = token;
    nextToken;
    tmp = readT(I,sigma);
    switch tmpToken.type
     case 'mul'
      T = T .* tmp;
     case 'div'
      T = T ./ tmp;
     case 'pow'
      T = T .^ tmp;
    end
  end
  
function PE=readPE(I,sigma)
  global token;

  PE = readE(I,sigma);
  if in(token.type, ['rpa'])
    nextToken;
  else
    error('Missing right parentheses.');
  end
  
function F=readF(I,sigma)
  global token;

  if in(token.type, ['num'])
    F = token.value;
    nextToken;
  else
    F = readD(I,sigma);
  end

function D=readD(I,sigma)
  global token;
  
  if in(token.type, ['ims'])
    D  = I;
    xs = 0;
    ys = 0;
    nextToken;
  else
    error('A factor should either be a real or a partial derivative of L.');
  end
  if in(token.type, ['var'])
    tmp = token.lexeme;
    xs = sum(upper(tmp)=='X');
    ys = sum(upper(tmp)=='Y');
    nextToken;
  end
  D  = real(ifft2(scale(fft2(I),sigma,ys,xs)));
  
function nextToken
  global tokens;
  global iTokens;
  global token;
  
  token = tokens(iTokens);
  if ~in(token.type, ['EOF'])
    iTokens = iTokens + 1;
  end
  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
% Lex'er and help functions

function tokens=lex(inStream)
  global dataIn;
  global iData;
  global ch;
  
  iData  = 1;
  dataIn = [inStream '¤'];
  
  [tokenName domain tokenType compound] = textread('tokens.def', '%s%s%s%d');

  for i=1:length(tokenType)
    tokenDefs(i) = struct('name', tokenName(i), 'domain', domain(i), 'type', tokenType(i), ...
			  'compound', compound(i));
  end
  
  numTokens = 0;
  ch = nextChar;
  while ch ~= '¤';
    tokenType = 0;
    for i=1:length(tokenDefs)
      if in(ch, tokenDefs(i).domain)
	tokenType = i;
      end
    end
    if tokenType==0
      msg = ['In input expression ' dataIn(1:end-1)  ...
             ' character ' ch ' does not belong to any known token type.'];
      error(msg);
    end
    numTokens = numTokens+1;
    tokens(numTokens) = readToken(tokenDefs(tokenType));
  end
  numTokens = numTokens + 1;
  tokens(numTokens) = struct('type', 'EOF', 'lexeme', '', 'value',0.0);

function token=readToken(tokenType)
  global ch;

  lexeme(1) = ch;
  i  = 2;
  ch = nextChar;
  
  if tokenType.compound
    while in(ch, tokenType.domain)
      lexeme(i) = ch;
      i = i + 1;
      ch = nextChar;
    end
  end
  
  token = struct('type', tokenType.name, 'lexeme', lexeme, 'value',0.0);
  if strcmp(tokenType.type, 'number') 
    token.value = str2double(lexeme);
  end
  
function nc = nextChar
  global iData;
  global dataIn;

  whitespace = ' \n';
  nc    = dataIn(iData);
  if nc ~= '¤' 
    iData = iData + 1;

    if in(nc, whitespace) 
      nc = nextChar;
    end
  end
  
function in=in(c, list)

  in = 0;
  if length(c) == 1
    if isempty(find(list==c))
      in = 0;
    else
      in = 1;
    end
  else
    for i=1:size(list,1)
      if strcmp(c, list(i,:))
	in = 1;
      end
    end
  end
  





