% test for error decay of cheeger selection

name = 'square';

n = 100;
options.radius = .9;

options.bound = 'per';  % boundary conditions for gradient
options.order = 2;      % order for gradient

rep = 'results/error/';
if not(exist(rep))
    mkdir(rep);
end

% exact curve
[c0,C0] = compute_exact_cheeger(name,n,options);

% load the image
mask = rescale( load_image(name,n, options) );
if mask(1)==1
    mask = 1-mask;
end
mask = double(mask>0.5);

options.tv_norm = 'l2';
weight_tv = 'none';

% make a total variation of 1 for the reference
tv = compute_total_variation(mask, options);

% reduce TV norm by projection
options.mask = mask;

niter = 50000*8;
ncomp = 60; % computing hausdorff distance
tv_divide = linspace(1,1,ncomp)*100;
m = length(tv_divide);


C0 = C0 / compute_total_variation(C0, options) * tv/max(tv_divide);

options.niter = round(niter/ncomp);

tv_error = [];
err_inf = [];
err_l2 = [];

% initial guess
options.x = [];
options.x = mask;

clf;
imageplot(-mask, 'Original', 1,m+1,1);
for i=1:m
    tau = tv/tv_divide(i);
    [M1,err_tv] = perform_tv_projection(mask,tau,options);
    tv1 = compute_total_variation(M1, options);
    
    % extract cheeger
    MM = rescale(M1);
    MM = sin(MM*pi/2).^2;
    c = perform_contour_extraction( rescale(MM), .5)*(size(MM,1)-1)+1;
    c(1,:) = rescale(c(1,:)); c(2,:) = rescale(c(2,:));
    % compute errors
    tv_error(end+1) = (tv1-tau)/tv1;
    % [err_l2(end+1),err_inf(end+1)] = compute_hausdorff_distance(c,c0,-1);
    err_l2(end+1) = norm( C0-M1, 'fro' )/n;
    err_inf(end+1) = mmax( abs(C0-M1) );
    
    % display the two curves
    lw = 2;
    clf;
    hold on;
    h = plot(c(1,:), c(2,:), 'k');
    set(h, 'LineWidth', 2);
    h = plot(c0(1,:), c0(2,:), 'k--');
    set(h, 'LineWidth', 2);
    hold off;
   
    clf;
    imageplot(M1);
    drawnow;
    
    options.x = M1;
end

% explore the level sets
tlist = linspace(mmax(M1)*.02,mmax(M1),60);
errc = [];
for i=1:length(tlist)
    progressbar(i,length(tlist));
%    c = perform_contour_extraction( M1, tlist(i))*(size(M1,1)-1)+1;
    errc(end+1) = norm( double(M1>tlist(i)) - rescale(C0), 'fro' );  % compute_hausdorff_distance(c,c0,-1);   
end
[tmp,i] = min(errc);
c = perform_contour_extraction( M1, tlist(i)); % *(size(M1,1)-1)+1;

% display the two curves
lw = 2;
clf;
hold on;
h = plot(c(1,:), c(2,:), 'k');
set(h, 'LineWidth', lw);
h = plot(c0(1,:), c0(2,:), 'k--');
set(h, 'LineWidth', lw);
hold off;
saveas(gcf, [rep name '-cheeger-progression.eps'], 'epsc');


clf;
h = plot(log10(err_l2), 'k');
set(h, 'LineWidth', lw);
axis tight;
saveas(gcf, [rep name '-error-decay.eps'], 'epsc');


[c0,C0] = compute_exact_cheeger(name,256,options);
warning off;
imwrite(1-C0*.6, [rep name '-cheeger-exact.png'], 'png');
warning on;