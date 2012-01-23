function N = NMS(S, O, high, low)

% N = NMS(S, O, high, low)
%
% This function performs Non-Maximal Suppression on Compass Operator output
% (Strength (S) and Orientation (O)), along with hysteresis thresholding.
% The NMS implementation is similar to the Vista implementation.  I don't
% know how good that one is.  O can have many sheets, but only the
% topmost one (O(:,:,1)) is used.  
%
% Unspecified values for 'high' and/or 'low' default to zero.  If 'low' is 
% greater than 'high', a warning will be printed, and the result will be
% empty.  If either argument is a vector, N will be a cell array of binary 
% images.
%
% The output is a binary edge map as a logical array.
%
% Mark Ruzon, 1999-2004
  
if nargin < 3
  high = 0.0;
end

if nargin < 4
  low = 0.0;
end

M = zeros(size(S));
Z = O(:,:,1) >= 0;  % Remove -1's
O = Z .* rem(O(:,:,1) + 90, 180);  % Must rotate orientations 90 degrees

% B and D are opposite 4-neighbors of each pixel.  A is the 8-neighbor on
% one side of B or the other (depending on the orientation).  C has a similar
% relationship to D.

T = (O(2:end-1,2:end-1) >= 45) & (O(2:end-1,2:end-1) < 135);
Z = Z(2:end-1,2:end-1);
B = (T & Z) .* S(1:end-2,2:end-1) + (~T & Z) .* S(2:end-1,3:end);
D = (T & Z) .* S(3:end,2:end-1) + (~T & Z) .* S(2:end-1,1:end-2);

T1 = T & Z & (O(2:end-1,2:end-1) < 90);
T2 = T & Z & ~T1;
T3 = ~T & Z & (O(2:end-1,2:end-1) < 45);
T4 = ~T & Z & ~T3;
A = T1 .* S(1:end-2,3:end) + T2 .* S(1:end-2,1:end-2) + ...
    T3 .* S(1:end-2,3:end) + T4 .* S(3:end,3:end);
C = T1 .* S(3:end,1:end-2) + T2 .* S(3:end,3:end) + ...
    T3 .* S(3:end,1:end-2) + T4 .* S(1:end-2,1:end-2);

% U determines how much to weight the pair of pixels on each side
% Dbad keeps infinities from showing up as a result of tangents
Dbad = (O(2:end-1,2:end-1) == 90) | (O(2:end-1,2:end-1) == 0);
Ogood = Dbad .* -1 + ~Dbad .* O(2:end-1,2:end-1);
Utan = tan(Ogood*pi/180);
Ucot = 1 ./ Utan;
Utan(Dbad) = 0;
Ucot(Dbad) = 0;

U = Z .* abs(T .* Ucot + ~T .* Utan);

G1 = U .* A + (1 - U) .* B;
G2 = U .* C + (1 - U) .* D;

M(2:end-1,2:end-1) = (S(2:end-1,2:end-1) > G1) & (S(2:end-1,2:end-1) >= G2) ...
                     & Z;

% Now perform hysteresis thresholding
M2 = M .* S;
N = cell(length(high),length(low));

for k = 1 : length(low)
	N2 = bwlabel(M2 > low(k), 8); % forms connected components
	
	% remove edges of length 1 by counting connected 8-neighbors
	N3 = N2 > 0;
	N4 = zeros(size(N3));
	N4(1:end-1,1:end-1) = N4(1:end-1,1:end-1) + N3(2:end,2:end);
	N4(1:end-1,:) = N4(1:end-1,:) + N3(2:end,:);
	N4(1:end-1,2:end) = N4(1:end-1,2:end) + N3(2:end,1:end-1);
	N4(:,2:end) = N4(:,2:end) + N3(:,1:end-1);
	N4(2:end,2:end) = N4(2:end,2:end) + N3(1:end-1,1:end-1);
	N4(2:end,:) = N4(2:end,:) + N3(1:end-1,:);
	N4(2:end,1:end-1) = N4(2:end,1:end-1) + N3(1:end-1,2:end);
	N4(:,1:end-1) = N4(:,1:end-1) + N3(:,2:end);
	
	N2 = N2 .* (N4 > 0);
	
	for j = 1 : length(high)
		if (low(k) > high(j))
			warning(['High threshold ' num2str(high(j)) ' is lower than low threshold ' num2str(low(k))]);
			N{j,k} = [];
			continue;
		end

		% Form a list of only those components with a high value
		U = unique(N2 .* (M2 > high(j)));
		if U(1) == 0  % permissible because unique sorts the elements
			U = U(2:end);
		end
		N{j,k} = ismember(N2,U) > 0; % map that list back onto the image
	end
end

if prod(size(N)) == 1
	N = N{1,1}; % remove from cell structure
end

    



