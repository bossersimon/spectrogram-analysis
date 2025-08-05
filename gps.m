%% GPS signal plot

M = readmatrix("recordings/GPS/sensorLog_20250701_02.txt");

time= M(:,1);
lat = M(:,3);
lon = M(:,4);
alt = M(:,5);

figure
geoplot(lat,lon,"-.",'LineWidth',1,'Color','r')
geobasemap streets
title("Car route")

lat0 = lat(1);
lon0 = lon(1);
alt0 = alt(1);

ref = [lat0 lon0 alt0];

wgs84 = wgs84Ellipsoid('meters');
[x, y, zUp] = geodetic2enu(lat, lon, alt, lat0, lon0, alt0, wgs84);

% Plot the full route
% plot3(x, y, zUp, 'b-'); % plot route in blue line
% hold on;
% 
% % Marker for start point (first coordinate)
% plot3(x(1), y(1), zUp(1), 'go', 'MarkerSize', 10, 'LineWidth', 2); % green circle
% 
% % Marker for end point (last coordinate)
% plot3(x(end), y(end), zUp(end), 'ro', 'MarkerSize', 10, 'LineWidth', 2); % red circle
% axis equal;

 % figure;
 % plot(x, y, 'b.-');
 % hold on
 % plot(x(1),y(1),'r*')
 % plot(x(end),y(end),'g*')
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
disp(vel_est);

tot_dist = sum(distances);
disp(tot_dist)