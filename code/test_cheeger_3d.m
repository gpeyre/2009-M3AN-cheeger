% cheeger in 3D

path(path, 'images/');
path(path, 'toolbox/');

name = 'twospheres';
name = 'spheres-tube';
name = 'cubes-tube';
name = 'cone';
name = 'multi-cones';
name = 'cube';
n = 60;


options.bound = 'per';  % boundary conditions for gradient
options.order = 2;      % order for gradient

options.tv_norm = 'l1';
options.tv_norm = 'linf';
options.tv_norm = 'l2';

[mask,M0] = load_3d_shape(name, n, options);

clf;
plot_isosurface(M0);

% make a total variation of 1 for the reference
tv = compute_total_variation(mask);

% reduce TV norm by projection
% clear options;
options.bound = 'per';
options.niter = 1000;
options.mask = mask;

tv_divide = [4];
m = length(tv_divide);

options.x = [];
for i=1:m
    tau = tv/tv_divide(i);
    [M1,err_tv,err_l2] = perform_tv_projection(mask,tau,options);
    tv1 = compute_total_variation(M1, options);
    err = (tv1-tau)/tv1;
    disp( ['Final TV error: ' num2str( err ) '.'] );
    options.x = M1;
    clf;
    rho = mean(M1(:));
    rho = max(M1(:))/2;
%    imageplot(M1, ['TV/' num2str(tv_divide(i)) ', (err=' num2str(err*100, 2) '%)'], 1,m+1,i+1);
end

rep = 'results/cheeger-3d/';
if not(exist(rep))
    mkdir(rep);
end

% you should modify eta for scaling
eta = n/5; 
eta = 0;
ax = [1+eta n-eta 1+eta n-eta 1+eta n-eta];
rho = .1;
clf;
subplot(1,2,1);
plot_isosurface(M0,.5); axis(ax);
title('Shape');
subplot(1,2,2);
plot_isosurface(M1, rho); axis(ax);
title('Cheeger');
saveas(gcf, [rep name '-' options.tv_norm '-cheeger-3d.png'], 'png');

clf;
plot_isosurface(M0, .5); axis(ax);
saveas(gcf, [rep name '-shape.png'], 'png');
clf;
plot_isosurface(M1, rho); axis(ax);
saveas(gcf, [rep name '-' options.tv_norm '-cheeger.png'], 'png');