%% GPS signal plot

M = readmatrix("recordings/GPS/sensorLog_20250701_03.txt");

time= M(3:end-1,1);
lat = M(3:end-1,3);
lon = M(3:end-1,4);
alt = M(3:end-1,5);

figure
geoplot(lat,lon,"-.")
geobasemap topographic
title("Car route")


lat0 = lat(1);
lon0 = lon(1);
alt0 = alt(1);

ref = [lat0 lon0 alt0];

wgs84 = wgs84Ellipsoid('meters');
[x, y, zUp] = geodetic2enu(lat, lon, alt, lat0, lon0, alt0, wgs84);

% figure;
% plot(x, y, 'b.-');
% axis equal;
% xlabel('East [m]');
% ylabel('North [m]');

dx = diff(x);
dy = diff(y);

%%
time_sec = time/1000;
datetime_array = datetime(time_sec, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
sample_times = seconds(datetime_array - datetime_array(1));
disp(sample_times);

distances = (sqrt(dx.^2 + dy .^2));

vel_est = distances./diff(sample_times);
disp(vel_est)

tot_dist = sum(distances);
disp(tot_dist)