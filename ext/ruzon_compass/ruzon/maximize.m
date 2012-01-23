function B = maximize(A, c, n)

% function B = maximize(A, c, n)
% maximize returns a matrix where each entry is the maximum of itself and
% all of its neighbors.  c can be 4, 8 (default) to specify connectivity.
% n is a recursive parameter; if n is 2, for instance, maximization will
% take place over a 5x5 neighborhood rather than a 3x3 one.
     
  if (nargin < 2) | (c ~= 4)
    c = 8;
  end
  
  if (nargin < 3)
    n = 1;
  end
  
  
  [h,w] = size(A);

  % Do the 4-connected part first
  B = max(A, [A(2:end,:); zeros(1,w)-Inf]);
  B = max(B, [zeros(1,w)-Inf; A(1:end-1,:)]);
  B = max(B, [A(:,2:end) zeros(h,1)-Inf]);
  B = max(B, [zeros(h,1)-Inf A(:,1:end-1)]);
  
  if (c == 8)
    B = max(B, [A(2:end,2:end) zeros(h-1,1)-Inf; zeros(1,w)-Inf]);
    B = max(B, [zeros(1,w)-Inf; zeros(h-1,1)-Inf A(1:end-1,1:end-1)]);
    B = max(B, [zeros(h-1,1)-Inf A(2:end,1:end-1); zeros(1,w)-Inf]);
    B = max(B, [zeros(1,w)-Inf; A(1:end-1,2:end) zeros(h-1,1)-Inf]);
  end
  
  if (n > 1)
    B = maximize(B, c, n - 1);
  end
  
