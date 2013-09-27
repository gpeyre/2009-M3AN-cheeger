% test for cheeger optimization by projecting on TV ball

path(path, 'images/');
path(path, 'toolbox/');

if not(exist('name'))
    name = 'L';
    name = 'ellipse-fat';
    name = 'ellipse-thin';
    name = 'square-tube-2';
    name = 'square-tube-3';
    name = 'square-tube-1';
    name = 'polygon-8';
    name = 'polygon-10';
    name = 'polygon-12';
    name = 'pacman';
    name = 'square-hole';
    name = 'chicken';
    name = 'bird';
    name = 'camel';
    name = 'device';
    name = 'apple';
    name = 'square';
end
clear options;

% load the image
n = 100;

options.radius = .9;
mask = rescale( load_image(name,n, options) );
if mask(1)==1
    mask = 1-mask;
end
mask = double(mask>0.5);


options.tv_norm = 'l2';
options.tv_norm = 'l1';
options.tv_norm = 'linf';

weight_l2 = 'ramp';
weight_l2 = 'none';
weight_tv = 'none';
w2 = ones(n); w1 = w2;
[Y,X] = meshgrid(1:n,1:n);
if strcmp(weight_l2, 'ramp')
	w2 = rescale(X,1,5).^2;
end
options.weight_l2 = w2;
options.weight_tv = w1;

% make a total variation of 1 for the reference
tv = compute_total_variation(mask, options);

% reduce TV norm by projection
options.bound = 'per';
options.niter = 40000;
options.niter = 40000*2;
options.niter = 10000;
options.mask = mask;

tv_divide = [10];
m = length(tv_divide);

clf;
imageplot(-mask, 'Original', 1,m+1,1);
for i=1:m
    tau = tv/tv_divide(i);
    [M1,err_tv,err_l2] = perform_tv_projection(mask,tau,options);
    tv1 = compute_total_variation(M1, options);
    err = (tv1-tau)/tv1;
    disp( ['Final TV error: ' num2str( err ) '.'] );
    options.x = M1;
    imageplot(-M1, ['TV/' num2str(tv_divide(i)) ', (err=' num2str(err*100, 2) '%)'], 1,m+1,i+1);
end

rep = 'results/';
if not(exist(rep))
    mkdir(rep);
end
saveas(gcf, [rep name '-cheeger.png'], 'png');
repimg = 'results/imgs/';
if not(exist(repimg))
    mkdir(repimg);
end
warning off;
imwrite(rescale(-mask), [repimg name '-original.png'], 'png');
imwrite(rescale(-M1), [repimg name '-cheeger.png'], 'png');
warning on;