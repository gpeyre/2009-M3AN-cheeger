function [M,M0] = load_3d_shape(name, n, options)

% load_3d_shape - load an implicit function in 3D
%
%   [M,M0] = load_3d_shape(name, n, options)
% 
%   Copyright (c) 2007 Gabriel Peyre

x = linspace(-1,1,n);
[X,Y,Z] = ndgrid(x,x,x);

r = getoptions(options, 'radius', .8);
c = getoptions(options, 'center', [0 0 0]);
        
M0 = [];
switch name
    case 'cube'
        if length(r)==1
            r = [r r r];
        end
        M = (abs(X-c(1))<r(1))&(abs(Y-c(2))<r(2))&(abs(Z-c(3))<r(3));
    case 'sphere'
        M = (X-c(1)).^2+(Y-c(2)).^2+(Z-c(3)).^2 <=r^2;   
        M0 = smooth3(M,'box',3);     
    case 'twospheres'
        delta = getoptions(options, 'delta',.35);
        options.radius = .95-delta;
        options.center = [0 0 delta];
        M1 = load_3d_shape('sphere', n, options);
        options.center = [0 0 -delta];
        M2 = load_3d_shape('sphere', n, options);
        M = double( M1 | M2);
        M0 = smooth3(M,'box',3);
    case 'cylinder'
        M = sqrt(X.^2+Y.^2)<r.^2;
        h = getoptions(options, 'height', .7);
        M = double(M.*(abs(Z)<=h));
        M0 = smooth3(M,'box',3);
    case 'spheres-tube'
        options.radius = .55;
        options.delta = .55;
        M = load_3d_shape('twospheres', n, options);
        options.radius = .4;
        M = double( M | load_3d_shape('cylinder', n, options) );
        M0 = smooth3(M,'box',3);
    case 'cubes-tube'
        delta = .6;        
        options.radius = .98-delta;
        options.center = [0 0 delta];
        M1 = load_3d_shape('cube', n, options);
        options.center = [0 0 -delta];
        M2 = load_3d_shape('cube', n, options);
        options.center = [0 0 0];
        options.radius = [.16 .16 .7];
        M3 = load_3d_shape('cube', n, options);
        M = double(M1|M2|M3);
    case 'cone'
        ax = getoptions(options, 'axis', [1 1 1]);
        A = rand(3); A(:,1) = ax;
        [A,R] = qr(A); 
        A = A * sign( sum(A(:,1).*ax(:)) );            
        pos = [X(:),Y(:),Z(:)]';
        pos = A'*pos;
        X = reshape(pos(1,:),n,n,n);
        Y = reshape(pos(2,:),n,n,n);
        Z = reshape(pos(3,:),n,n,n);
        D = sqrt(Y.^2+Z.^2);
        eta = getoptions(options, 'eta', .95);
        M = (X>=0) .* ( D <= r*(1-X/eta));
        M0 = smooth3(M,'box',3);
    case 'multi-cones'
        options.radius = .3;
        cones_layout = getoptions(options, 'cones_layout', 'ico');
        vertex = compute_base_mesh(cones_layout);
        M = zeros(n,n,n);
        for i=1:size(vertex,2)
            progressbar(i, size(vertex,2));
            options.axis = vertex(:,i);
            M = M | load_3d_shape('cone', n, options);
        end
        M0 = smooth3(M,'box',3);      
    otherwise
        error('Unknown');
end

if isempty(M0)
    M0 = M;
end