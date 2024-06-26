function [x,y,u,v] = get_circ_vectors(ocean)



% get grid coords and inverse of transport matrix
lon = ocean.lon(ocean.i(ocean.Ib));
lat = ocean.lat(ocean.j(ocean.Ib));
TM  = ocean.TMB{1};
TMi = TM';
vol = ocean.M(ocean.Ib);

% convert lat and lon to cartesian coordinates to avoid edge effects
az = deg2rad(lon+180); % convert to radians
el = deg2rad(lat); % convert to radians
[X,Y,Z] = sph2cart(az,el,1);

% Apply transport matrix
% concentration-weighted coordinate transform
X2 = sum(TM.*(vol.*X)',2)./(TM*vol); 
Y2 = sum(TM.*(vol.*Y)',2)./(TM*vol); 
Z2 = sum(TM.*(vol.*Z)',2)./(TM*vol); 

% convert back to lat and lon for plotting
[az2,el2,r2] = cart2sph(X2,Y2,Z2);

lon2 = rad2deg(az2); % convert to degrees
lat2 = rad2deg(el2); % convert to degrees

% recentre on 0 meridian
lon2(lon2<=0)=lon2(lon2<=0)+360;

% calculate velocity vectors
v = (lat-lat2).*vol;
u = (lon-lon2).*vol;

% saturate at 90th percentile
v = sign(v).*min(abs(v),prctile(abs(v),90));
u = sign(u).*min(abs(u),prctile(abs(u),90));

y = lat;
x = lon;

end