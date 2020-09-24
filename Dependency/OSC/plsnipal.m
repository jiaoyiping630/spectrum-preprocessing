function [p,q,w,t,u] = plsnipal(x,y)
%PLSNIPAL NIPALS algorithm for PLS
%  This program does the nipals algorithm for PLS
%  It is generally run as a subprogram of PLS, since it
%  calculates only one latent variable.
%
%I/O: [p,q,w,t,u] = plsnipal(x,y);
%
%See also: MODLGUI, PLS, SIMPLS

%Copyright Eigenvector Research, Inc. 1991-98
%Modified BMW 11/93

[my,ny] = size(y);
if ny > 1
  ssy = sum(y.^2);
  [ymax,yi] = max(ssy);
  u = y(:,yi);
else
  u = y(:,1);
end
conv = 1;
told = x(:,1);
count = 1.0;
%  Specify the conversion tolerance
while conv > 1e-10
  count = count + 1;
  w = (u'*x)';
  w = (w'/norm(w'))';
  t = x*w;
  if ny == 1
    q = 1;
    break
  end
  q = (t'*y)';
  q = (q'/norm(q'))';
  u = y*q;
  conv = norm(told - t);
  told = t;
  if count >= 50.0
    disp('Algorithm Failed to Converge after 50 Iterations')
    break;
  end
end
p = (t'*x/(t'*t))';
p_norm=norm(p);
t = t*p_norm;
w = w*p_norm;
p = p/p_norm;