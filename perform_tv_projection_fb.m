function [M1,err,tv] = perform_tv_projection_fb(M, tv0, options)

options.null = 0;

niter = getoptions(options, 'niter', 1000);
gradop = getoptions(options, 'gradop', @grad);
divop = getoptions(options, 'divop', @div);
mu = getoptions(options, 'mu', 0.25);
disp = getoptions(options, 'display', 1);

n = size(M,1);
u = zeros(n,n,2);
err = [];
tv = [];
for i=1:niter
    progressbar(i,niter);
    % descent step
    M1 =  M - feval( divop, u, options);
    u = u - mu*feval( gradop, M1, options );
    % prox step
    u = u - perform_l1ball_projection( u, tv0*mu );
    % display
    if disp
        imageplot(M1); drawnow;
    end
    % monitor error
        err(i) = 1/2*norm(M1,'fro')^2;  % + tv0 * max(max( sqrt(sum(u.^2,3)) ));
        tv(i) = compute_total_variation(M1);
end