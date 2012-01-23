function Y = DrawLine(X,r,c,theta,length)

% Y = DrawLine(X,r,c,theta,length)
%
% This function takes a binary matrix X and adds to it a line of pixels
% numbered one, starting from the (r,c) entry and moving in direction theta
% for a distance length.  The distance between horizontally or vertically
% adjacent pixels is 1.

% First determine the quadrant we are in
quadrant = floor(theta/90) + 1;
if (quadrant < 1 | quadrant > 4)
     error('Theta must be in the range [0, 360)');
end

% Create matrix of angles (Z)
switch (quadrant)
     case 1,
       Xr = (-0.5:length+0.5);
       Yr = (length+0.5:-1:-0.5)';
     case 2,
       Xr = (-length-0.5:0.5) ;
       Yr = (length+0.5:-1:-0.5)';
     case 3,
       Xr = (-length-0.5:0.5);
       Yr = (0.5:-1:-length-0.5)';
     case 4,
       Xr = (-0.5:length+0.5);
       Yr = (0.5:-1:-length-0.5)';
end

Z = atan2(repmat(Yr,1,length+2),repmat(Xr,length+2,1));
% Fix inconsistencies at theta = 180;
if (quadrant == 2)
     Z(end,:) = Z(end,:) + 2 * pi;
end
if (quadrant == 3)
     Z(1,:) = Z(1,:) - 2 * pi;
end

% Create a mask of entries within "length" of the corner 
V = repmat(0:length,length+1,1);
M = rot90(V .* V + flipud(V' .* V') <= length * length, quadrant-1);

% Convert theta to work with the atan2 function
theta2 = theta * pi / 180;
if (theta2 >= pi)
     theta2 = theta2 - 2 * pi;
end

% W holds the line; it must be masked to the proper length
switch (quadrant)
     case 1,
       W = Z(1:end-1,1:end-1) > theta2 & Z(2:end,2:end) < theta2;
     case 2,
       W = Z(2:end,1:end-1) > theta2 & Z(1:end-1,2:end) < theta2;
     case 3,
       W = Z(1:end-1,1:end-1) < theta2 & Z(2:end,2:end) > theta2;
     case 4,
       W = Z(2:end,1:end-1) < theta2 & Z(1:end-1,2:end) > theta2;
end

W = W & M;

% Now W has to be OR'ed into X to produce Y; of course, W may not lie inside X
switch (quadrant)
     case 1,
       rows = min(length+1,r);
       cols = min(length+1,size(X,2)-c+1);
       X(r-rows+1:r,c:c+cols-1) = X(r-rows+1:r,c:c+cols-1) | ...
                                  W(end-rows+1:end,1:cols);
     case 2,
       rows = min(length+1,r);
       cols = min(length+1,c);
       X(r-rows+1:r,c-cols+1:c) = X(r-rows+1:r,c-cols+1:c) | ...
                                  W(end-rows+1:end,end-cols+1:end);
     case 3,
       rows = min(length+1,size(X,1)-r+1);
       cols = min(length+1,c);
       X(r:r+rows-1,c-cols+1:c) = X(r:r+rows-1,c-cols+1:c) | ...
                                  W(1:rows,end-cols+1:end);
     case 4,
       rows = min(length+1,size(X,1)-r+1);
       cols = min(length+1,size(X,2)-c+1);
       X(r:r+rows-1,c:c+cols-1) = X(r:r+rows-1,c:c+cols-1) | ...
                                  W(1:rows,1:cols);
end

Y = X;



