function gimage(img)

% function gimage(img)
%
% Takes a greyscale or color image and displays it

if (size(img,3) == 1)
  image(cat(3,repmat(img,[1 1 3])));
elseif (size(img,3) == 3)
  image(img);
else
  error('gimage display error');
end

axis image;
