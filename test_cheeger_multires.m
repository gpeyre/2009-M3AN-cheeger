% test for cheeger optimization by projecting on TV ball
% Use a multi-grid like algorithm.

path(path, 'images/');
path(path, 'toolbox/');

% clear options;
if not(exist('name'))
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
    name = 'bird';
    name = 'camel';
    name = 'device';
    name = 'L';
    name = 'square';
    name = 'constant';
    name = 'apple';
    name = 'chicken';
    name = 'square';
end

options.bound = 'per';  % boundary conditions for gradient
options.order = 2;      % order for gradient

% load the image
n = 400;

options.radius = .9;
options.sigma = 0;
switch name
    case { 'chicken' 'camel' 'bird' 'apple'}
        options.sigma = 8;
end
mask = rescale( load_image(name,n,options) );
mask = double(mask>0.5);
if mask(1)==1
    mask = 1-mask;
end

% some weight
weight_l2 = 'ramp';
weight_l2 = 'periodic';
weight_l2 = 'constant';
weight_tv = 'periodic';
weight_tv = 'bump';
weight_tv = 'constant';
options.freq = 5;
options.center = [.9 -.66];
options.bump_size = [.1,.1];
w1 = rescale(load_image(weight_tv, n, options),1,10);
w2 = rescale(load_image(weight_l2, n, options),1,12);
options.weight_l2 = w2;
options.weight_tv = w1;

if not(isfield(options, 'tv_norm'))
options.tv_norm = 'linf';
options.tv_norm = 'l1';
options.tv_norm = 'l2';
end

% make a total variation of 1 for the reference
tv = compute_total_variation(mask, options);

% reduce TV norm by projection
options.niter = 40000;
options.niter = 8000;
options.mask = mask;


tv_divide = 10;
m = length(tv_divide);

nresol = 5;
resol = linspace(0.4, 1, nresol);
t = linspace(1,0,nresol).^3;
niter = round( rescale(t,200,7000) ); % 7000

M1 = mask;
for i=1:length(resol)
    % switch of resolution
    nn = round(resol(i)*[n n]); % new dimension
    M1 = perform_image_resize(M1, nn);
    mask1 = double( perform_image_resize(mask, nn )>0.5 );
    M1 = M1.*mask1;
    if i==1
        M1 = mask1;
    end
    % resize also the weights
    options.weight_l2 = perform_image_resize(w2,nn);
    options.weight_tv = perform_image_resize(w1,nn);

    tau = nn(1)/n* tv/tv_divide;
    options.mask = mask1;
    options.niter = niter(i);
    options.x = M1;
    [M1,err_tv,err_l2] = perform_tv_projection(mask1,tau,options);
    tv1 = compute_total_variation(M1, options);
    err = (tv1-tau)/tv1;
    disp( ['TV Error, step ' num2str(i) '/' num2str(length(resol)) ': ' num2str( err*100 ) '%.'] );
end

MM = rescale(M1);
MM = sin(MM*pi/2).^2;
c = perform_contour_extraction( rescale(MM), .5)*(size(MM,1)-1)+1;

rep = 'results/';
if not(exist(rep))
    mkdir(rep);
end
repimg = 'results/imgs/';
if not(exist(repimg))
    mkdir(repimg);
end
nameadd = [name '-' options.tv_norm];
if not(strcmp(weight_l2, 'constant')) || not(strcmp(weight_tv, 'constant'))
    nameadd = [name '-' weight_l2 '-' weight_tv ];
end 

clf;
if strcmp(weight_l2,'constant')==0
    imageplot(-rescale(w2,1,2).*mask, 'L2 weight', 1,2,1);
elseif strcmp(weight_tv,'constant')==0
    imageplot(-rescale(w1,1,2).*mask, 'TV weight', 1,2,1);
else
    imageplot(-mask, 'Original', 1,2,1);
end
subplot(1,2,2);
hold on;
imagesc(-MM); axis image; axis off;
h = plot(c(2,:),c(1,:));
set(h, 'LineWidth', 2);
hold off;
axis ij;
% imageplot(M1, ['TV/' num2str(tv_divide(i)) ', (err=' num2str(err*100, 2) '%)'], 1,m+1,i+1);
title(['Cheeger, norm ' options.tv_norm]);
saveas(gcf, [rep nameadd '-cheeger.png'], 'png');

% plot cheeger along
clf; hold on;
imagesc(-MM); axis image; axis off;
h = plot(c(2,:),c(1,:));
set(h, 'LineWidth', 2);
hold off; axis ij;
saveas(gcf, [repimg nameadd '-cheeger-curve.png'], 'png');


warning off;
imwrite(rescale(-mask), [repimg name '-original.png'], 'png');
imwrite(rescale(-MM), [repimg nameadd '-cheeger.png'], 'png');
if strcmp(weight_l2,'constant')==0
    imwrite(rescale(-rescale(w2,1,2).*mask), [repimg nameadd '-l2.png'], 'png');
end
if strcmp(weight_tv,'constant')==0
    imwrite(rescale(-rescale(w1,1,2).*mask), [repimg nameadd '-tv.png'], 'png');
end
warning on;