function I = imreadgray(file) 
% IMREADGRAY  Reads an image as gray-scale
%   I=IMREADGRAY(FILE) reads the image FILE and converts the result to
%   a gray scale image (in double format normalized with range in
%   [0,1]).

I=double(imread(file))/256 ;

if(size(I,3) > 1)
  I = rgb2gray( I ) ;
end
