function se = strel_sphere( r )
% create 3D structuring element : sphere of radius r
r = r-1;
[x y z] = meshgrid( -r:r, -r:r, -r:r );
se = x.^2 + y.^2 + z.^2 <= r.^2;
end