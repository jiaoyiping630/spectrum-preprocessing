function [m,ssq,p,q,w,t,u,b] = pls(x,y,lv,out)
%PLS Partial least squares regression via NIPALS algorithm
%  Inputs are the scaled predictor block (x), scaled
%  predicted block (y), the number of latent variables 
%  to be calculated (maxlv), and an optional variable (out) to
%  suppress intermediate output [out=0 suppresses output].
%  Outputs are the the matrix of regression vectors or
%  matrices (b), the fraction of variance used in the x and
%  y matrices (ssq), x loadings (p), y loadings (q), x
%  weights (w), x scores (t), y scores (u), and inner
%  relation coefficients (bin),
%
%  The regression matrices are ordered in b such that each
%  ny (number of y variables) rows correspond to the regression
%  matrix for that particular number of latent variables.
%
%I/O: [b,ssq,p,q,w,t,u,bin] = pls(x,y,maxlv,out);
%
%See also: CROSSVAL, PCR, PLSNIPAL, PLSDEMO, SIMPLS

%Copyright Eigenvector Research, Inc. 1991-98
%Modified NBG 4/96,9/96,12/97
%BMW Checked on MATLAB 5
%Modified BMW 12/98  added rank check

if nargin < 4
  out   = 1;
end
[mx,nx] = size(x);
[my,ny] = size(y);
if nx < lv
  error('No. of LVs must be <= no. of x-block variables')
end 
p       = zeros(nx,lv);
q       = zeros(ny,lv);
w       = zeros(nx,lv);
t       = zeros(mx,lv);
u       = zeros(my,lv);
b       = zeros(1,lv);
ssq     = zeros(lv,2);
ssqx    = sum(sum(x.^2)');
ssqy    = sum(sum(y.^2)');
olv = lv;
rankx = rank(x);
if rankx < olv
  lv = rankx;
  if out == 1
    disp('  ')
    sss = sprintf('Rank of X is %g, which is less than lv of %g',lv,olv);
    disp(sss);
    sss = sprintf('Calculating %g LVs only',lv);
    disp(sss);
  end
end
for i = 1:lv
  [pp,qq,ww,tt,uu] = plsnipal(x,y);
  b(1,i)   = uu'*tt/(tt'*tt);
  x        = x - tt*pp';
  y        = y - b(1,i)*tt*qq';
  ssq(i,1) = (sum(sum(x.^2)'))*100/ssqx;
  ssq(i,2) = (sum(sum(y.^2)'))*100/ssqy;
  t(:,i)   = tt(:,1);
  u(:,i)   = uu(:,1);
  p(:,i)   = pp(:,1);
  w(:,i)   = ww(:,1);
  q(:,i)   = qq(:,1);
end
if olv > lv
  ssq(lv+1:end,2) = ssq(lv,2);
end
ssqdif   = zeros(lv,2);
ssqdif(1,1) = 100 - ssq(1,1);
ssqdif(1,2) = 100 - ssq(1,2);
for i = 2:olv
  for j = 1:2
    ssqdif(i,j) = -ssq(i,j) + ssq(i-1,j);
  end
end
ssq = [(1:olv)' ssqdif(:,1) cumsum(ssqdif(:,1)) ssqdif(:,2) ...
  cumsum(ssqdif(:,2))];
if out ~= 0
  disp('  ')
  disp('       Percent Variance Captured by PLS Model   ')
  disp('  ')
  disp('           -----X-Block-----    -----Y-Block-----')
  disp('   LV #    This LV    Total     This LV    Total ')
  disp('   ----    -------   -------    -------   -------')
  format = '   %3.0f     %6.2f    %6.2f     %6.2f    %6.2f';
  for i = 1:olv
    tab = sprintf(format,ssq(i,:)); disp(tab)
  end
  disp('  ')
end
m = zeros(olv*ny,nx);
m(1:lv*ny,:) = conpred(b,w,p,q,lv);
if ny > 1
  for i = 2:olv
    j   = (i-1)*ny+1;
    i0  = j-ny;
    m(j:i*ny,:) = m(j:i*ny,:) + m(i0:(i-1)*ny,:);
  end
else
  m     = cumsum(m,1);
end

function m = conpred(b,w,p,q,lv)
%CONPRED Converts PLS models to regression vectors
%  The inputs are the inner-relation coefficients (b),
%  the x-block weights (w), the x-block loadings (p), the
%  y-block loadings (q) and the number of latent variables
%  to consider (lv). The output is a matrix (m) of the 
%  contribution of each latent variable to the final 
%  regression vector.
%
%  CONPRED works with either single or multiple variable y-block 
%  PLS models. If there is only 1 y-block variable each row of 
%  the output matrix corresponds to the contribution from each lv to
%  the y-block prediction.  If there are N y-block variables
%  each block of N rows corresponds to the contribution from
%  each lv to the prediction. See CONPRED1 for obtaining
%  final models.
%
%  The I/O format: m = conpred(b,w,p,q,lv);
%
%  See also: CONPRED1, PLS, PLSPRED

%  Copyright Eigenvector Research 1993-98
%  Modified BMW 5/94

[mq,nq] = size(q);
[mw,nw] = size(w);
if nw ~= lv
  if lv > nw
    s = sprintf('Original model has a maximum of %g LVs',nw);
    disp('  '), disp(s)
	s = sprintf('Calculating vectors for %g LVs only',nw);
	disp(s), disp('  ')
	lv = nw;
  else
    w = w(:,1:lv);
	q = q(:,1:lv);
	p = p(:,1:lv);
	b = b(:,1:lv);
  end
end
m = zeros(mq*lv,mw);
if mq == 1
  m = (w*inv(p'*w)*diag(b))';
else
  mp = (w*inv(p'*w)*diag(b))';
  for i = 1:lv
    mpp = mp(i,:);
    m((i-1)*mq+1:i*mq,:) = diag(q(:,i))*mpp(ones(mq,1),:);
  end
end