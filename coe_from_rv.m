%% this program is a function calculates the orbital elements from the initial position and velocity vectors 

function coe_from_rv()
%% calculate all six orbital emements from the initial position and velocity
    %
    % Noha Faty 
    % 21 august 2024
    %           
    % function orbitElements(R,V,mu)
    %
    % Purpose:  This function calculates the classic orbital elements.
    % 
    % Inputs:   o R  - A 1x3 vector describing the initial position of the
    %                  satellite.
    %           o V  - A 1x3 vector describing the initial velocity of the
    %                  satellite.
    %           o mu - Standard gravitationl parameter of the central body
    %                  [OPTIONAL]. Defaults to Earth (398600 [km^3/s^2])
    %
    % Output:   o h     - Specific angular momentum
    %           o e     - eccentricity
    %           o i     - orbital inclination
    %           o omega - right ascension of the ascending node
    %           o w     - argument of perigee
    %           o theta - true anomaly
    %
    
    clear r v vr H h i k N n E e omega w theta; clc;
    

%% User Input for Position Vector R
R = input ('Enter the initial position vector [R1, R2, R3] in km :');

%% User Input for Velocity  Vector V
V = input ('Enter the initial velocity vector [V1, V2, V3] IN km/s:');

%% User Input for Gravitational Parameter mu (optional)
mu = input ('Enter the gravitational parameter mu in km^3 / s^2 (or press Enter for Earth default )');

% Set default value for Earth's gravitational parameter if not provided
if isempty (mu)
    mu = 398600;
end



%% calculate the magnitudes of position and velocity vectors 
r = norm(R);
v = norm (V);

%% calculate  the specific angular momentum (h)

H = cross (R,V);
h = norm(H);

%% calculate  the eccentricity 

E = 1/mu *((v^2 - mu/r)*R - dot(R, V)*V);
e = norm(E);

%% calculate  the inclination

K = [0,0,1];  %% k unit vector
i = acos (dot(K, H)/h);

%% calculate the right ascention of ascending node (RAAN)
N = cross (K, H);  %% node line vector 
n = norm(N);
I = [1,0,0];   %% i unit vector


if N >=0
    omega = acos(dot(I, N)/ n);
else 
    omega = 360 -acos(dot(I, N)/ n);
end

%% calculate the argument of perigee (w) 
if E >=0 
    w= acos (dot (N, E)/ n*e);
else 
    w= 360-  acos (dot (N, E)/ n*e);
end 




%% calculate the true anomaly (theta)
theta = acos(dot(E, R)/ e*r);



  %% Display Results


fprintf('Specific Angular Momentum (h): %.2f km^2/s\n', h);
fprintf(' Eccentricity (e) : %.2f \n ', e);
fprintf('the inclination (i) : %.2f degree\n', i);
fprintf('the right ascention of ascending node (RAAN or omeaga ) : %.2f degree \n', omega);
fprintf('the argument of perigee (w) : %.2f degree \n ', w);
fprintf('the true anomaly  (theta) : %.2f degree ', theta);
end 
