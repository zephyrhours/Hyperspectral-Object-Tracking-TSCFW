function hypershow(im,bands)
if (nargin<2)
imshow(hyper2im(im),'border','tight','initialmagnification','fit');
else
imshow(hyper2im(im,bands),'border','tight','initialmagnification','fit');
end
% imshow(hyper2im(im));