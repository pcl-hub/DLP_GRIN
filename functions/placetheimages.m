function fullimage = placetheimages(image,qty,spacings,offcenter,fullsize,pixelsize)

% spacings is a 1x2 variable, unit: mm
% where the 1st and the 2nd elements are the spacing between...
% print items in Dimension 1 (y) and Dimension 2 (x), respectively 
%
% offcenter is also a 1x2 varible, unit: mm

if isempty(fullsize); fullsize = [1080 1920]; end
if isempty(offcenter); offcenter = [0 0]; end % unit: mm
if isempty(pixelsize); pixelsize = 64.8; end % unit: um

fullimage = zeros(fullsize);
imagesize = size(image);

c = 1000/pixelsize; % coeffienct for mm to pixel conversion

pts = pointcoordinates(qty);
spacings = ones(qty,1)*spacings;
pts = ceil(pts.*spacings*c);

pixel1 = ceil(fullsize/2-imagesize/2+offcenter*c);

for i = 1:qty
    yrange = pixel1(1)+pts(i,1)+(1:imagesize(1));
    xrange = pixel1(2)+pts(i,2)+(1:imagesize(2));
    fullimage(yrange,xrange) = image;
end



function pts = pointcoordinates(qty)

switch qty
    case 1
        pts = [0 0];
    case 2
        pts = [-1 1; 1 -1]*sqrt(2)/4;
    case 3
        pts = [ sqrt(3)/4 -1/2;...
                sqrt(3)/4 1/2;...
               -sqrt(3)/4 0 ];
    case 4
        pts = [-1 1; 1 1; 1 -1; -1 -1]*1/2;
    case 5
        pts = [1 0; sin(0.1*pi) cos(0.1*pi);...
                    sin(0.1*pi) -cos(0.1*pi);...
                    -cos(0.2*pi) sin(0.2*pi);...
                    -cos(0.2*pi) -sin(0.2*pi)];
    case 6
        pts = [sqrt(3) -1; sqrt(3) 1; 0 2;...
               -sqrt(3) 1; -sqrt(3) -1; 0 -2]*1/2;
    case 7
        pts = [sqrt(3) -1; sqrt(3) 1; 0 2;...
               -sqrt(3) 1; -sqrt(3) -1; 0 -2; 0 0]*1/2;
    case 8
        pts = [1 -1; 1 0; 1 1; 0 -1; 0 1; -1 -1; -1 0; -1 1];
    case 9
        pts = [1 -1; 1 0; 1 1; 0 -1; 0 0; 0 1; -1 -1; -1 0; -1 1];
end
