function [c,M] = compute_exact_cheeger(name,n, options)

if nargin<2
    n = 200;
end

options.null = 0;

switch name
    case 'square'
        R = 1/(2+sqrt(pi));
        
        %% curve
        t = linspace( R,1-R,n ); tr = t(end:-1:1);
        a = ones(1,n); z = a*0;
        theta = -pi/2 + linspace(0,2*pi,4*n);
        w = [ cos(theta); sin(theta) ]*R;
        w(1,:) = w(1,:)+1-R; w(2,:) = w(2,:)+R;
        c = [];
        c = [c [t;z]];
        c = [c w(:,1:end/4) ];
        c = [c [a;t]];
        
        w(2,:) = w(2,:)+1-2*R;
        c = [c w(:,end/4+1:end/2) ];
        c = [c [tr;a]];
        
        w(1,:) = w(1,:)-1+2*R;
        c = [c w(:,end/2+1:3*end/4) ];        
        c = [c [z;tr]];

        w(2,:) = w(2,:)-1+2*R;
        c = [c w(:,3*end/4:end) ];
        
        
        r = getoptions(options, 'radius', .9);
        c = c*r + (1-r)/2;

        %% image
        x = linspace(0,1,n)/r - (1-r)/(2*r);
        [Y,X] = meshgrid(x,x);
        C = { [R R], [R 1-R], [1-R R], [1-R 1-R] };
        M = zeros(n);
        for i=1:length(C)
            d = (X-C{i}(1)).^2 + (Y-C{i}(2)).^2;
            M(d<R^2) = 1;
        end
        M( X>0 & X<1 & Y>R & Y<1-R )=1;
        M( X>R & X<1-R & Y>0 & Y<1 )=1;
end
