% [S, O, A, U] = compassmex(I, R, spc, Is, Angles, wedges)
%
% compassmex: the Generalized Compass Operator for images
%
% This MEX-file computes various information related to edges and corners in 
% an image.
%
% Output Arguments:
% 
%   REQUIRED:
%     S -- the Strength of the Compass Operator, lying in the interval [0,1].
%          Its dimensions are H x W, where H is the height and W the width of
%          the subimage chosen.
%
%   OPTIONAL:
%     O -- the Orientation of the edge or corner, in degrees.  For edge 
%          detection, this number lies in [0,180) and is measured with respect
%          to the positive X-axis.  For corners, the interval is [0,360) and
%          refers to the orientation of the "right" side of the corner.  E.g.,
%          a 90 degree corner at orientation 45 subtends an angle whose rays
%          point in the 45 and 135 degree orientations.  Its dimensions are
%          H x W x D, where D is the maximum number of responses (currently 3).
%          at a single point in the image.  If a point (r,c) registers C < D
%          responses, then O(r,c,C+1) == -1.
%     A -- the Abnormality of the Compass Operator at that point.  This is 
%          the minimum value over all orientations.  It is H x W.
%     U -- the Uncertainty of the Compass Operator.  This is the number of
%          degrees over which, at each orientation, the compass operator gives
%          similar responses.  It is also H x W x D.
%
% Input Arguments:
% 
%   REQUIRED:
%     I -- an image of dimensions M x N x 3.  If I is of class uint8, the
%          values are treated as RGB values.  If I is of class double, the
%          values are assumed to be in the CIE-Lab color space.
%     R -- the standard deviation of the Gaussian used to weight pixels in
%          the Compass Operator (the operator's radius is 3*R).  If R is a 
%          vector, then the output arguments S, O, A, and U will be a row 
%          array of matrices, one for each value.
%   
%   OPTIONAL:
%     spc -- the spacing at each scale.  This is a vector of the same length
%            as R stating how many pixels apart applications of the Compass
%            Operator are to be placed.  If spc is a scalar, it applies to 
%            all scales.  Default: 1.
%     Is -- subimage to apply the Compass Operator to.  If I is a scalar, the
%           entire image is used.  If I is of the form [r c], then the Compass
%           Operator is applied only at that point.  If I is of the form
%           [TR LC BR RC] (top row, left column, etc.), then only that subimage
%           will be used.  NOTE: if the subimage is too close to the image 
%           border, so that the maximum radius would cause the circle to fall
%           outside the image, the subimage is trimmed accordingly. 
%           Default: 1 (the entire image).
%     Angles -- the angles subtended by the Generalized Compass Operator.  If
%               Angles is a vector, the Compass Operator will detect corners
%               and edges at different angles, and the output arguments will
%               have multiple rows (if R and Angles are both vectors, the 
%               output will be arrays of matrices).  Default: 180 (edges only).
%     nwedges -- the number of wedges (orientations) in one-quarter of the
%                circle.  nwedges defines the set of possible orientations and
%                angles, so the values in the Angles vector must match.
%                Default: 6 (90/6 = 15 degree increments)
%     maxclusters -- the maximum number of clusters to be created by the
%                clustering algorithm.  Default: 10.

