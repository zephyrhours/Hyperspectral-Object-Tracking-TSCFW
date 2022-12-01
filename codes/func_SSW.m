function result= func_SSW(patch,sigma_d,sigma_e)
% Spatail Speetral for Hyperspectral Target Tracking
% Author: Zephyr Hou
% Time: 2021-11-25
%
%% Main Function
sz =size(patch) ;
% patch = imresize(patch, [sz(1)/cell_size, sz(2)/cell_size]);
% sz =size(patch);

[rs, cs] = ndgrid((1:sz(1)) - floor(sz(1)/2), (1:sz(2)) - floor(sz(2)/2));
w1 =  sqrt((rs.^2 + cs.^2))/sigma_d;       % spatial weighted

sz=sz(1:2);
cen_pix = patch(round(sz(1)/2),round(sz(2)/2),:);              % bands x 1
w2 = sqrt(sum((patch-cen_pix).^2,3))/sigma_e;   % spectral weighted

% w1=(w1-min(w1(:)))/(max(w1(:))-min(w1(:)));
% w2=(w2-min(w2(:)))/(max(w2(:))-min(w2(:)));

result = exp(-0.5*(w1+w2));

end

