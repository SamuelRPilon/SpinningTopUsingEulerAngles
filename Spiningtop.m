
%% MOTION OF THE SPINING TOP 
% This function defines the system of ordinary differential equations (ODEs)
% that describe the motion of the spinning top.
%
% INPUTS:
% - t: current time
% - k: vector of current values for the six state variables [phi, theta, psi, phi_dot, theta_dot, psi_dot]
% - friction: the friction factor of the spinning top
%
% OUTPUT:
% - difft: vector of the time derivative of the six state variables
%
% The function uses the input values to compute the time derivatives of
% the six state variables in order to solve the ODEs. The equations for
% these time derivatives were derived from the equations of motion for
% the spinning top, taking into account the effects of precession, nutation,
% and spin. The friction factor is also included in the equations.

function difft = Spiningtop(t,k,friction)

difft(1) = k(4);
difft(2) = k(5);
difft(3) = k(6);
difft(4) = (3*k(6)*k(5) - 13*k(4)*k(5)*cos(k(2)))/(8*sin(k(2)));
difft(5) = (5*cos(k(2))*sin(k(2))*k(5)^2)/8 - (3*k(4)*sin(k(2))*((k(6))/8 ...
    + (4895*k(4)*cos(k(2)))/24));
difft(6) = k(4)*k(5)*sin(k(2)) - difft(4)*cos(k(2)) - (20000*friction)/9;
difft = difft';

end