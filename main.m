%% SAMUEL R PILON
% AE 426 - Project-2
% Date - 03/29/2023
% For personal use of SAMUEL PILON only. Not to be distributed
%% QUESTION STATEMENT  
%{
    Derive the equation of motion using Euler's rotational equation in terms of
    precession (φ), nutation (θ), and spin (ψ) angles. The center-of-mass of the top
    is located at ⃗r_G = lˆe3 in Fig. 1.
    
    Consider that the friction with the oor is adding resistance in such a way
    that there is a linear decrease of the spin with a torque applied in the ˆe3 axis
    leading to the spin stop in 1 minute. Graph the solution to the non-linear
    equations using MATLAB as a function of precession, nutation and spin and
    stop the motion when the surface of the top (5º with respect to the axis of
    symmetry) touches the oor. Explain the physical motion of the top with the
    simulation.
%}
%% MAIN FUNCTION (SPINNING TOP)
%{
    This program uses numerical methods to solve the equations of motion 
    known as Euler's equations, which describe the rotational motion of a 
    spinning top. By integrating these equations over time, the program 
    can predict the future orientation and angular velocity of the top based 
    on its initial conditions and the forces acting upon it.

    The quaternion allows for efficient and accurate updates to the top's 
    orientation as it undergoes rotational motion, providing a time history of 
    its orientation throughout the simulation.
 
    User subfunction required: Spinningtop.m

%}
 clear all; close all; clc
 
%% SET GIVEN VALUES
% This section of code initializes variables and constants related to a
% Spining top problem 

A = 12.e-4; % moment of inertia about one axis of rotation (kg*m^2)
C = 4.5e-4; % moment of inertia about another axis of rotation (kg*m^2)
m = 0.5; % mass (kg)
g = 9.79; % acceleration due to gravity (m/s^2)
d = 0.05; % distance from origin to center of mass (m)
initial_nutations = [1e-8]; % initial nutations in degrees
tspan = linspace(0,60,100); %time span of 60 seconds incremented in steps of 100 
syms  precession spin nutation precessiondot spindot nutationdot thetaddot psiddot phiddot mu 



%% EULERS EQUATIONS OF MOTION 
% Define left-hand and right-hand sides of Euler equations of motion
Momentright_body = [A*thetaddot, A*(phiddot*sin(nutation) + precessiondot*nutationdot*cos(nutation)), ...
    C*(psiddot + phiddot*cos(nutation) - precessiondot*nutationdot*sin(nutation))] + ...
    cross([nutationdot, precessiondot*sin(nutation), precessiondot*cos(nutation)], ...
    [A*nutationdot, A*precessiondot*sin(nutation), C*(spindot + precessiondot*cos(nutation))]); % torques acting on the spinning top
Momentleft_body = [m*g*d*sin(nutation), 0, -mu]; % forces acting on the spinning top


%% DIRECT COSINE MATRIX 
%Relate the coordinates of a vector in one reference 
%frame to the coordinates of the same vector in another reference frame.

%Represents a rotation around the z-axis by an angle nutation.
DCM_1_nutation = [cos(nutation),-sin(nutation),0;sin(nutation),cos(nutation),0;0,0,1];

%Represents a rotation around the x-axis by an angle precession.
DCM_2_precession = [1,0,0;0,cos(precession),-sin(precession);0,sin(precession),cos(precession)];

%Represents a rotation around the z-axis by an angle spin.
DCM_3_spin = [cos(spin),-sin(spin),0;sin(spin),cos(spin),0;0,0,1];

%create a single 3x3 matrix DCM that represents the overall 
%rotation from the inertial frame to the body frame.
DCM = DCM_1_nutation*DCM_2_precession*DCM_3_spin;

%Represents the angular velocity of the body frame relative 
%to the inertial frame, expressed in the body frame
w_body = [nutationdot,precessiondot*sin(nutation),precessiondot*cos(nutation)]; 

%Computes the angular velocity vector by transforming from the body frame to the inertial frame 
%using the DCM matrix.
w_inertial = w_body*DCM;

% Solve for the angular accelerations
thetaddot = solve(Momentleft_body(1) == Momentright_body(1), thetaddot)
phiddot = solve(Momentleft_body(2) == Momentright_body(2), phiddot)
psiddot = solve(Momentleft_body(3) == Momentright_body(3), psiddot)


%% ODE 45 (ORDINARY DIFFERENTIAL EQUATION SOLVER)

figure('Name', 'Euler Angles and Their Rates', 'color', [1 1 1])

opt = odeset('RelTol', 1e-9,'AbsTol',1e-9); % Sets the relative and absolute error tolerances

for p = 1
    initial_nutations = initial_nutations(p); % Setting the initial nutation angle for this iteration
    Mu_initial = 0; % Setting the initial friction factor to zero

    % Initial conditions in order: precession, nutation, spin, precessiondot, nutationdot, spindot
    ini = [0, initial_nutations*pi/180, 0, 0, 0, 100*pi/180]; % Converting initial nutation angle from degrees to radians

    % Increasing friction factor until the spin stops
    defining_angle = [9,9;9,9]; % Defining an angle matrix to fall in loop
    while defining_angle(end,end) > 0 % While the spin angular velocity is greater than zero
        Mu_initial = Mu_initial + 0.001; % Increase the friction factor
        [t,defining_angle] = ode45(@(t,k) Spiningtop(t,k,Mu_initial), tspan, ini, opt); % Solve the ODE using the current friction factor
    end

    % This section of the code is used to identify when the nutation angle reaches 85 degrees
    % and then to extract the corresponding time and angle values up to that point, 
    % effectively stopping the motion.
    t_index = find(defining_angle(:,2) > 85*pi/180); % Find the index where the nutation angle reaches 85 degrees
    t_index = t_index(1); % Get the first index (i.e., the earliest time) where the nutation angle exceeds 85 degrees
    t = t(1:t_index); % Keep time values up to the point where nutation angle exceeds 85 degrees
    defining_angle = defining_angle(1:t_index,:); % Keep angle values up to the point where nutation angle exceeds 85 degrees
%% PLOTS 

    % Plots
    figure(p); % Create a new figure for each initial nutation angle
    titles = {'Precession (φ)', 'Nutation (θ)', 'Spin (ψ)', 'Angular Velocity (Precession)', 'Angular Velocity (Nutation)', 'Angular Velocity (Spin)'};
    for i = 1:6
    % Create subplot for current angle/velocity
    subplot_index = subplot(2,3,i);
    subplot_handle = subplot(subplot_index);
    
    % Plot angle/velocity vs time with red color and linewidth of 1.5
    plot_handle = plot(t, defining_angle(:,i), 'r-', 'LineWidth', 1.5);
    
    % Customize subplot properties
    set(subplot_handle, 'GridLineStyle', ':', 'GridColor', 'k', 'XMinorGrid', 'on');
    xlabel('Time (s)');
    if i < 4
        ylabel('Angle (rad)');
    else
        ylabel('Angular Velocity (rad/s)');
    end
    title(titles(i));
    
    end

    if p == 1
        sgtitle('Spinning Top Simulation for Small Initial Nutation Angle (0.01e-5)'); % Add a super title to the first figure
    end
end
