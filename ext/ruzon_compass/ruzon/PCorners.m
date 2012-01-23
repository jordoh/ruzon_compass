function [A, B] = PCorners(img, S, O, radius, Zlow, Zhigh)

% This function finds relative extrema in (x,y,theta) and does non-maximal
% suppression on compass operator strength for edges.  A holds the corners
% drawn in by the first and B the lines drawn by the second.  S and O are 
% assumed to be lists of matrices representing output at one scale for a set 
% of evenly spaced orientations between 0 and 180 (excluding 0).  S is
% assumed to have the 'plot' flag on; that is, each element of S should
% be M x N x 2*L, except the last one.
% 
% Here is a much-needed example:
% img = imread('img.gif');
% region = img(358:390,153:173,:);
% [S,O] = compassmex(img,2,1,[357 152 390 173],15:15:180,'plot');
% [X,Y] = PCorners(region,S,O,6);
%
% Note the offset between region and the rectangle passed into
% compassmex.  It is there to make the corner lines match up on the image.

L = length(S);
Y = 180 / L;
zfactor = 4;                     % Amount to zoom images by
img2 = imresize(img,zfactor);
C = zeros((size(S{1},1) - 1) .*zfactor, (size(S{1},2) - 1) .* zfactor);
B = zeros(size(C));
Sb = 0.31;              % Corner strength threshold (b = y-intercept)
Sm = 0.03;             % Corner strength increment (m = slope)
Pb = 1.97;              % Degree of match threshold
Pm = 0.00;             % Degree of match increment
Tb = 0.4;              % Edge strength threshold
Tm = 0.02;             % Edge strength increment
checklen = radius;

if (nargin < 5)
  Zlow = 0.2;
end
if (nargin < 6)
  Zhigh = 0.5;
end

if (size(C,1) ~= size(img2,1)) | (size(C,2) ~= size(img2,2))
  error('img2 must be one pixel smaller than S in each direction');
end

% Create Strength Array (Str) over all alpha's and threshold it.
% Also find the maximum over each 3x3 neighborhood at each theta
m = 0;
for i=1:L
  Str{i} = max(S{i},[],3);
  m = m + sum(sum((Str{i} > Sb + Sm * i)));
  Smax{i} = S{i};
  for j=1:size(Smax{i},3)
    Smax{i}(:,:,j) = maximize(Smax{i}(:,:,j));
  end
end

% Create Data Structure for Good Corners
E = zeros(m,8) - 1;
e = 0;

% Find strong edge points
K = NMS(Str{L},O{L},Zhigh,Zlow);

% Find initial corner candidates
for i=2:L-2

  [I,J] = find(Str{i} > Sb + Sm * i);
  for j=1:length(I)
    Strength = Str{i}(I(j),J(j));
    o = 1;
    Ori = O{i}(I(j),J(j),o) / Y + 1;  % Get wedge #
    while (Ori >= 0)
      if (Ori ~= floor(Ori)) % Orientation falls between wedges
	OriPlus = ceil(mod(Ori,2*L));
	OriMinus = mod(floor(Ori)-1,2*L)+1;
	Ori = 0;
      else
	OriPlus = mod(Ori,2*L)+1;
	OriMinus = mod(Ori-2,2*L)+1;
      end
      
      % Check if the point is a maximum over a 3x3x3 neighborhood, 
      % including adjacent theta values
      if (~Ori | Strength >= Smax{i}(I(j),J(j),Ori)) & ...
	(Strength >= Smax{i}(I(j),J(j),OriPlus)) & ...
	(Strength >= Smax{i}(I(j),J(j),OriMinus)) % & ...
	%(~Ori | Strength >= Smax{i-1}(I(j),J(j),Ori)) & ...
	%(Strength >= Smax{i-1}(I(j),J(j),OriPlus)) & ...
	%(Strength >= Smax{i-1}(I(j),J(j),OriMinus)) & ...
	%(~Ori | Strength >= Smax{i+1}(I(j),J(j),Ori)) & ...
	%(Strength >= Smax{i+1}(I(j),J(j),OriPlus)) & ...
	%(Strength >= Smax{i+1}(I(j),J(j),OriMinus))

	% Position (r,c), Right & Left Orientations, Strength
	E2 = [I(j) J(j) O{i}(I(j),J(j),o) ...
		  rem(O{i}(I(j),J(j),o) + Y * i, 360)];
	
	% Compute the dot product of the two sides of the corner (DR &
        % DL) with the edge vector most of the way down the sides.
	XR = round(E2(1) - checklen * sin(E2(3)*pi/180));
	YR = round(E2(2) + checklen * cos(E2(3)*pi/180));
	if (XR >= 1) & (XR <= size(Str{L},1)) & ...
		  (YR >= 1) & (YR <= size(Str{L},2))
	  DR = min(abs(E2(3)-O{L}(XR,YR)), ...
		   min(360 - abs(E2(3)-O{L}(XR,YR)), ...
		       min(abs(mod(E2(3)+180,360)-O{L}(XR,YR)), ...
			   360 - abs(mod(E2(3)+180,360)-O{L}(XR,YR)))));
	end
	XL = round(E2(1) - checklen * sin(E2(4)*pi/180));
	YL = round(E2(2) + checklen * cos(E2(4)*pi/180));
	if (XL >= 1) & (XL <= size(Str{L},1)) & ...
		  (YL >= 1) & (YL <= size(Str{L},2))
	  DL = min(abs(E2(4)-O{L}(XL,YL)), ...
		   min(360 - abs(E2(4)-O{L}(XL,YL)), ...
		       min(abs(mod(E2(4)+180,360)-O{L}(XL,YL)), ...
			   360 - abs(mod(E2(4)+180,360)-O{L}(XL,YL)))));
	end
	
	% Compute the degree of match of the corner and edges (P)
	if (XR >= 1) & (XR < size(Str{L},1)) & ...
		  (XL >= 1) & (XL < size(Str{L},1)) & ...
		  (YR >= 1) & (YR < size(Str{L},2)) & ...
		  (YL >= 1) & (YL < size(Str{L},2))
		P = cos(DL*pi/180) + cos(DR*pi/180);
		P2 = Str{L}(XL,YL) * cos(DL*pi/180) * ...
		     Str{L}(XR,YR) * cos(DR*pi/180);
		%[E2 Strength P]
		
	  if (E2(1) >= 1) & (E2(1) < size(Str{L},1)) & ...
	    (E2(2) >= 1) & (E2(2) < size(Str{L},2)) & (P > Pb + Pm * i)
	    if (E2(3) < E2(4))
	      Mid = (E2(3) + E2(4)) / 2;
	    elseif (E2(3) >= 360 - E2(4))
	      Mid = (E2(3) + E2(4) - 360) / 2;
	    else
	      Mid = (E2(3) + E2(4) + 360) / 2;
	    end
	    MidPerp = rem(Mid + 90, 360);
	    OuterLine = zeros(size(S{L},1),size(S{1},2));
	    InnerLine = OuterLine;
	    Len = radius * (sec(Y*i*(pi/180)/2) - tan(Y*i*(pi/180)/2));
	    OuterLine = DrawLine(OuterLine,E2(1),E2(2),Mid, ...
				 min(ceil(1.2*Len),ceil(Len+1)));
	    InnerLine = DrawLine(InnerLine,E2(1),E2(2),Mid, ...
				 max(ceil(0.8*Len),ceil(Len-1)));
	    CheckingArea = OuterLine & ~InnerLine;
	    AngleDiff = min(abs(O{L}(:,:,1)-MidPerp), ...
		   min(360 - abs(O{L}(:,:,1)-MidPerp), ...
		       min(abs(mod(O{L}(:,:,1)+180,360)-MidPerp), ...
			   360 - abs(mod(O{L}(:,:,1)+180,360)-MidPerp))));
	    if (all(all(Strength > Str{L} .* cos(AngleDiff*pi/180) .* ...
			CheckingArea)))
	    %if (max(max(Str{L} .* LC)) < Sb+Sm*L)
	      e = e + 1;
	      E(e,:) = [E2 Strength P Str{L}(XL,YL) Str{L}(XR,YR)];
	      %[E(e,:) Str{L}(XL,YL) cos(DL*pi/180) Str{L}(XR,YR) ...
	      %cos(DR*pi/180)]
	    end
	  end
	end
      end
      
      o = o + 1;
      if (o <= size(O{i},3))
	Ori = O{i}(I(j),J(j),o) / Y;
      else
	Ori = -1;
      end
    end  
  end
end

e
E(1:e,:);
% Determine which corners are close to which other corners
F = zeros(e);
AngleThresh = 10;
AngleSumThresh = 40;
for i=1:e
  for j=i+1:e
    if (norm(E(i,1:2) - E(j,1:2)) <= 3 * radius / 4)
      CWDiff = min(abs(E(i,3) - E(j,3)), 360 - abs(E(i,3) - E(j,3)));
      CCWDiff = min(abs(E(i,4) - E(j,4)), 360 - abs(E(i,4) - E(j,4)));
      if (CWDiff <= AngleThresh) | (CCWDiff <= AngleThresh) | ...
	    (CWDiff + CCWDiff <= AngleSumThresh)
	F(i,j) = 1;
      end
    end
  end
end
F = F + F' + eye(e);

% Compute transitive closure
for k=1:e
  Fnew = zeros(e);
  for i=1:e
    for j=1:e
      Fnew(i,j) = F(i,j) | (F(i,k) & F(k,j));
    end
  end
  F = Fnew;
end
    
% Partition the corner candidates
G = zeros(0,8);  % List of "true" corners
H = (1:e)';  % Index matrix
while (~isempty(F))
  marked = zeros(size(F,1),1);
  marked(1) = 1;
  for i=2:size(F,1)
    if (F(1,:) == F(i,:))
      marked(i) = 1;
    end
  end
  R = H(find(marked));
  E2 = E(R,:)
  [J, I] = max(E(R,5) + E(R,6)/2 + (E(R,7)+E(R,8))/2);
  G = [G; E(R(I),:)];
  F = F(find(~marked),:);
  H = H(find(~marked),:);
end

%G = E(1:e,:)   % To see all the corners
for j=1:size(G,1)
  [(G(j,1)-1)*zfactor (G(j,2)-1)*zfactor G(j,3:6) (G(j,7)+G(j,8))/2]
  if (G(j,1) > 0) & (G(j,1) <= size(img,1)) & (G(j,2) > 0) & ...
	(G(j,2) <= size(img,2))
    C = DrawLine(C,(G(j,1)-1)*zfactor+1,(G(j,2)-1)*zfactor+1,G(j,3), ...
		 radius*zfactor);
    C = DrawLine(C,(G(j,1)-1)*zfactor+1,(G(j,2)-1)*zfactor+1,G(j,4), ...
		 radius * zfactor);
  end
end


% Draw edges
[I,J] = find(K);
for j=1:length(I)
  if (K(I(j)+1,J(j)))
    B = DrawLine(B,(I(j)-1)*zfactor+1,(J(j)-1)*zfactor+1, 270, zfactor);
  end
  if (K(I(j),J(j)+1))
    B = DrawLine(B,(I(j)-1)*zfactor+1,(J(j)-1)*zfactor+1, 0, zfactor);
  end
  if (K(I(j)+1,J(j)+1))
    B = DrawLine(B,(I(j)-1)*zfactor+1,(J(j)-1)*zfactor+1, 315, ...
		 round(zfactor * sqrt(2)));
  end
  if (K(I(j)-1,J(j)+1))
    B = DrawLine(B,(I(j)-1)*zfactor+1,(J(j)-1)*zfactor+1, 45, ...
		 round(zfactor * sqrt(2)));
  end
end

figure
subplot(121);
A = img2 .* (C == 0) + C .* 255;
gimage(A);
title({['Partial Matching'];['R=' num2str(radius) ' S=' num2str(Sb) '+' ...
       num2str(Sm) '*i P=' num2str(Pb) '+' num2str(Pm) '*i M=S+P/2+E/2 CL=' ...
       num2str(checklen)]});
subplot(122);
B = img2 .* (B == 0) + B .* 255;
gimage(B);
title(['Edges extracted Zlow=' num2str(Zlow)]);

