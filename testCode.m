numAngles = 32;
theta = linspace(0,2*pi,numAngles+1);%for circular calculation
theta = theta(1:end-1);

theta_degrees = theta/pi*180;
theta_degrees180 = mod(theta_degrees, 180);

theta2 = linspace(0,pi,numAngles+1);%for circular calculation
theta2 = theta2(1:end-1);
theta_degrees2 = theta2/pi*180;
theta_degrees2180 = mod(theta_degrees2, 180);

theta = [30, 120];
theta_rad = theta*pi/180;
theta2 = theta*2;
theta_rad2 = theta2*pi/180;

weight = ones(1,2);
theta_weigthed = circ_mean(theta_rad',weight');
theta_weigthed %2.37 %1.57
theta_weigthed/pi*180 % 136 %90
theta_weigthed2 = circ_mean(theta_rad2',weight');
theta_weigthed2 %1.52 %
theta_weigthed2/pi*180 % 87  %4

%original code
theta_weigthedOld = mod(theta_weigthed,2*pi);%from [-pi, pi] to [0 2pi]
theta_weigthedOld = theta_weigthedOld./2;%range 0 to pi.
theta_weigthedOld2 = mod(theta_weigthed2,2*pi);%from [-pi, pi] to [0 2pi]
theta_weigthedOld2 = theta_weigthedOld2./2;%range 0 to pi.

%fixed code
theta_weigthedNew = mod(theta_weigthed, pi);
theta_weigthedNew2 = mod(theta_weigthed2, pi);

theta_weigthedOld %1.18 %0.78
theta_weigthedOld/pi*180 %68 %45
theta_weigthedNew %2.37 %1.57
theta_weigthedNew/pi*180 %136 %90
theta_weigthedOld2 %0.76 %3.64
theta_weigthedOld2/pi*180 %43 %2
theta_weigthedNew2 %1.52 
theta_weigthedNew2/pi*180 %87

mean(theta) %3.04
mean(theta_degrees) %174
mean(theta_degrees180) %84
mean(theta2) %1.52
mean(theta_degrees2) %87
mean(theta_degrees2180) %87

orig = [1, 60, 120, 150, -150, -60];
orig180 = mod(orig, 180);
orig_rad = orig*pi/180;

%original code
orig_radOld = mod(orig_rad,2*pi);%from [-pi, pi] to [0 2pi]
orig_radOld = orig_radOld./2;%range 0 to pi.
orig_radOld = orig_radOld/pi*180;

%fixed code
orig_radNew = mod(orig_rad, pi);
orig_radNew = orig_radNew/pi*180;


orig_radOld
orig_radNew
