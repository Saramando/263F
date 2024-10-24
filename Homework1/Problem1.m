%% Related to Assignment 1
% Number of nodes
N = 3;
ndof = N * 2; % number of degrees of freedom
dt = 0.01; % second - Time step size
RodLength = 0.1; % l = 1 m
deltaL = RodLength / (N-1);
% Radii of spheres
R1 = 0.005; % meter
R2 = 0.025; % meter
R3 = 0.005; % meter
% Density
rho_metal = 7000; % kg/m^3
rho_f = 1000; % fluid
rho = rho_metal - rho_f;
r0 = 0.001; % meter - rod radius
Y = 1e9; % Young's modulus (Y instead of E for clarity)
g = 9.8; % m/s^2 - gravity
visc = 1000; % pa-s
totalTime = 10; % second - total simulation time
% Utility parameter
ne = N - 1; % number of edges
EI = Y * pi * r0^4 / 4; % Nm^2 - bending stiffness
EA = Y * pi * r0^2; % Newton
% Geometry - initial configuration
nodes = zeros(N,2);
for c=1:N % Loop over all the nodes
 nodes(c,1) = (c-1) * deltaL; % x coordinates
 nodes(c,2) = 0;
end
% Mass, M
M = zeros(ndof, ndof);
M(1,1) = 4/3*pi*R1^3*rho_metal;
M(2,2) = 4/3*pi*R1^3*rho_metal;
M(3,3) = 4/3*pi*R2^3*rho_metal;
M(4,4) = 4/3*pi*R2^3*rho_metal;
M(5,5) = 4/3*pi*R3^3*rho_metal;
M(6,6) = 4/3*pi*R3^3*rho_metal;
% Viscous damping matrix, C
C = zeros(6,6);
C1 = 6 * pi * visc * R1;
C2 = 6 * pi * visc * R2;
C3 = 6 * pi * visc * R3;
C(1,1) = C1;
C(2,2) = C1;
C(3,3) = C2;
C(4,4) = C2;
C(5,5) = C3;
C(6,6) = C3;
% Weight vector, W
W = zeros(ndof, 1);
W(2) = -4/3*pi*R1^3*rho*g;
W(4) = -4/3*pi*R2^3*rho*g;
W(6) = -4/3*pi*R3^3*rho*g;
% Initial DOF
q0 = zeros(ndof, 1);
for c=1:N % loop over nodes
 q0( 2*c-1 ) = nodes(c,1); % x1, x2, x3
 q0( 2*c ) = nodes(c,2); % y1, y2, y3
end
u0 = zeros(ndof, 1); % old velocity (initial velocity)
% tolerance
tol = EI/RodLength^2 * 1e-3; % small enouch force that can be neglected
% Time marching scheme
Nsteps = round(totalTime/dt);
% Storage for y-velocity of the middle node
all_mid_v = zeros(Nsteps, 1);
for c = 2:Nsteps

 fprintf('Time = %f\n', (c-1) * dt);

 % Guess
 q = q0; % New DOFs are initialized to be equal to old DOFs
 % Newton Raphson
 err = 10 * tol;
 while err > tol
 f = M / dt * ( (q-q0)/dt - u0 );
 J = M / dt^2;

 %
 % Elastic forces
 %
 % Linear spring between nodes 1 and 2
 xk = q(1);
 yk = q(2);
 xkp1 = q(3);
 ykp1 = q(4);
 l_k = deltaL;
 dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
 dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
 f(1:4) = f(1:4) + dF;
 J(1:4,1:4) = J(1:4,1:4) + dJ;

 % Linear spring between nodes 2 and 3
 xk = q(3);
 yk = q(4);
 xkp1 = q(5);
 ykp1 = q(6);
 l_k = deltaL;
 dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
 dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
 f(3:6) = f(3:6) + dF;
 J(3:6, 3:6) = J(3:6, 3:6) + dJ;

 % Bending spring at node 2
 xkm1 = q(1);
 ykm1 = q(2);
 xk = q(3);
 yk = q(4);
 xkp1 = q(5);
 ykp1 = q(6);
 curvature0 = 0;
 l_k = deltaL;
 dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
 dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
 f(1:6) = f(1:6) + dF;
 J(1:6, 1:6) = J(1:6, 1:6) + dJ;

 % Viscous force
 f = f + C * (q-q0) / dt;
 J = J + C / dt;

 % Weight
 f = f - W;

 % Update
 q = q - J \ f;

 err = sum ( abs(f) );

end
 % New velocity
 u = (q - q0) / dt;

 % Store some information
 all_mid_v(c) = u(4);

 % Plot
 figure(1);
 plot( q(1:2:end), q(2:2:end), 'ro-');
 axis equal
 xlabel('x [meter]');
 ylabel('y [meter]');
 title(['Structure Configuration at Time t = ', num2str(t), ' s']);
 drawnow

 % Update (new becomes old)
 q0 = q;
 u0 = u;
end
% Plot middle node downward velocity
figure(2);
timeArray = (1:Nsteps) * dt;
plot(timeArray, all_mid_v, 'k-');
xlabel('Time, t [sec]');
ylabel('Velocity (vertical) of middle node, v [m/s]');

function F = gradEs(xk, yk, xkp1, ykp1, l_k, EA)
%
% This function returns the derivative of stretching energy E_k^s with 
% respect to x_{k-1}, y_{k-1}, x_k, and y_k.
%
F = zeros(4,1);

F(1) = -(0.1e1 - sqrt((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k) * ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1) / l_k * (-0.2e1 * xkp1 + 0.2e1 * xk);
F(2) = -(0.1e1 - sqrt((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k) * ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1) / l_k * (-0.2e1 * ykp1 + 0.2e1 * yk);
F(3) = -(0.1e1 - sqrt((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k) * ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1) / l_k * (0.2e1 * xkp1 - 0.2e1 * xk);
F(4) = -(0.1e1 - sqrt((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k) * ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1) / l_k * (0.2e1 * ykp1 - 0.2e1 * yk);

F = 0.5 * EA * l_k * F;

end
function J = hessEs(xk, yk, xkp1, ykp1, l_k, EA)
%
% This function returns the 4x4 hessian of the stretching energy E_k^s with
% respect to x_k, y_k, x_{k+1}, and y_{k+1}.
%
J11 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (-2 * xkp1 + 2 * xk) ^ 2) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * ((-2 * xkp1 + 2 * xk) ^ 2) / 0.2e1 - 0.2e1 * (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1)) / l_k;
J12 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (-2 * ykp1 + 2 * yk) * (-2 * xkp1 + 2 * xk)) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * (-2 * xkp1 + 2 * xk) * (-2 * ykp1 + 2 * yk) / 0.2e1;
J13 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * xkp1 - 2 * xk) * (-2 * xkp1 + 2 * xk)) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * (-2 * xkp1 + 2 * xk) * (2 * xkp1 - 2 * xk) / 0.2e1 + 0.2e1 * (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1)) / l_k;
J14 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * ykp1 - 2 * yk) * (-2 * xkp1 + 2 * xk)) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * (-2 * xkp1 + 2 * xk) * (2 * ykp1 - 2 * yk) / 0.2e1;
J22 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (-2 * ykp1 + 2 * yk) ^ 2) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * ((-2 * ykp1 + 2 * yk) ^ 2) / 0.2e1 - 0.2e1 * (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1)) / l_k;
J23 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * xkp1 - 2 * xk) * (-2 * ykp1 + 2 * yk)) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * (-2 * ykp1 + 2 * yk) * (2 * xkp1 - 2 * xk) / 0.2e1;
J24 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * ykp1 - 2 * yk) * (-2 * ykp1 + 2 * yk)) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * (-2 * ykp1 + 2 * yk) * (2 * ykp1 - 2 * yk) / 0.2e1 + 0.2e1 * (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1)) / l_k;
J33 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * xkp1 - 2 * xk) ^ 2) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * ((2 * xkp1 - 2 * xk) ^ 2) / 0.2e1 - 0.2e1 * (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1)) / l_k;
J34 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * ykp1 - 2 * yk) * (2 * xkp1 - 2 * xk)) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * (2 * xkp1 - 2 * xk) * (2 * ykp1 - 2 * yk) / 0.2e1;
J44 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * ykp1 - 2 * yk) ^ 2) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * ((2 * ykp1 - 2 * yk) ^ 2) / 0.2e1 - 0.2e1 * (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1)) / l_k;

J = [J11 J12 J13 J14;
     J12 J22 J23 J24;
     J13 J23 J33 J34;
     J14 J24 J34 J44];

J = 0.5 * EA * l_k * J;

end

%% Copyright M. Khalid Jawed (khalidjm@seas.ucla.edu)
% You should use this code at your own risk. Copy and redistribution is not
% permitted. Written permission is required.

function dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI)

%
% This function returns the derivative of bending energy E_k^b with respect
% to x_{k-1}, y_{k-1}, x_k, y_k, x_{k+1}, and y_{k+1}.
%
% curvature0 is the "discrete" natural curvature [dimensionless] at node (xk, yk).
% l_k is the voronoi length of node (xk, yk).
% EI is the bending stiffness.
%

node0 = [xkm1, ykm1, 0];
node1 = [xk, yk, 0];
node2 = [xkp1, ykp1, 0];
%     m1e, 
m2e = [0 0 1];
%     m1f,
m2f = [0 0 1];
kappaBar = curvature0;

%% Computation of gradient of the two curvatures
gradKappa = zeros(6,1);

ee = node1 - node0;
ef = node2 - node1;

norm_e = norm(ee);
norm_f = norm(ef);

te = ee / norm_e;
tf = ef / norm_f;

% Curvature binormal
kb = 2.0 * cross(te, tf) / (1.0 + dot(te, tf));

chi = 1.0 + dot(te, tf);
tilde_t = (te + tf) / chi;
tilde_d2 = (m2e + m2f) / chi;

% Curvatures
kappa1 = kb(3); % 0.5 * dot( kb, m2e + m2f); % CHECKED

Dkappa1De = 1.0 / norm_e * (-kappa1 * tilde_t + cross(tf,tilde_d2));
Dkappa1Df = 1.0 / norm_f * (-kappa1 * tilde_t - cross(te,tilde_d2));

gradKappa(1:2, 1) = -Dkappa1De(1:2);
gradKappa(3:4, 1) = Dkappa1De(1:2) - Dkappa1Df(1:2);
gradKappa(5:6, 1) = Dkappa1Df(1:2);

%% Gradient of Eb
dkappa = kappa1 - kappaBar;
dF = gradKappa * EI * dkappa / l_k;
end

%% Copyright M. Khalid Jawed (khalidjm@seas.ucla.edu)
% You should use this code at your own risk. Copy and redistribution is not
% permitted. Written permission is required.

function dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI)

%
% This function returns the derivative of bending energy E_k^b with respect
% to x_{k-1}, y_{k-1}, x_k, y_k, x_{k+1}, and y_{k+1}.
%
% curvature0 is the "discrete" natural curvature [dimensionless] at node (xk, yk).
% l_k is the voronoi length of node (xk, yk).
% EI is the bending stiffness.
%

node0 = [xkm1, ykm1, 0];
node1 = [xk, yk, 0];
node2 = [xkp1, ykp1, 0];
%     m1e, 
m2e = [0 0 1];
%     m1f,
m2f = [0 0 1];
kappaBar = curvature0;

%% Computation of gradient of the two curvatures
gradKappa = zeros(6,1);

ee = node1 - node0;
ef = node2 - node1;

norm_e = norm(ee);
norm_f = norm(ef);

te = ee / norm_e;
tf = ef / norm_f;

% Curvature binormal
kb = 2.0 * cross(te, tf) / (1.0 + dot(te, tf));

chi = 1.0 + dot(te, tf);
tilde_t = (te + tf) / chi;
tilde_d2 = (m2e + m2f) / chi;

% Curvatures
kappa1 = kb(3); % 0.5 * dot( kb, m2e + m2f); % CHECKED

Dkappa1De = 1.0 / norm_e * (-kappa1 * tilde_t + cross(tf,tilde_d2));
Dkappa1Df = 1.0 / norm_f * (-kappa1 * tilde_t - cross(te,tilde_d2));

gradKappa(1:2, 1) = -Dkappa1De(1:2);
gradKappa(3:4, 1) = Dkappa1De(1:2) - Dkappa1Df(1:2);
gradKappa(5:6, 1) = Dkappa1Df(1:2);

%% Computation of hessian of the two curvatures
DDkappa1 = zeros(6, 6);
% DDkappa2 = zeros(11, 11);

norm2_e = norm_e^2;
norm2_f = norm_f^2;

tt_o_tt = tilde_t' * tilde_t; % must be 3x3. tilde_t is 1x3
tmp = cross(tf, tilde_d2);
tf_c_d2t_o_tt = tmp' * tilde_t; % must be 3x3
tt_o_tf_c_d2t = tf_c_d2t_o_tt'; % must be 3x3
kb_o_d2e = kb' * m2e; % must be 3x3
d2e_o_kb = kb_o_d2e'; % must be 3x3

Id3 = eye(3);
D2kappa1De2 ...
    = 1.0 / norm2_e * (2 * kappa1 * tt_o_tt - tf_c_d2t_o_tt - tt_o_tf_c_d2t) ...
    - kappa1 / (chi * norm2_e) * (Id3 - te'*te) ...
    + 1.0 / (4.0 * norm2_e) * (kb_o_d2e + d2e_o_kb);

tmp = cross(te, tilde_d2);
te_c_d2t_o_tt = tmp' * tilde_t;
tt_o_te_c_d2t = te_c_d2t_o_tt';
kb_o_d2f = kb' * m2f;
d2f_o_kb = kb_o_d2f';

D2kappa1Df2 ...
    = 1.0 / norm2_f * (2 * kappa1 * tt_o_tt + te_c_d2t_o_tt + tt_o_te_c_d2t) ...
    - kappa1 / (chi * norm2_f) * (Id3 - tf'*tf) ...
    + 1.0 / (4.0 * norm2_f) * (kb_o_d2f + d2f_o_kb);

D2kappa1DeDf ...
    = -kappa1/(chi * norm_e * norm_f) * (Id3 + te'*tf) ...
    + 1.0 / (norm_e*norm_f) * (2 * kappa1 * tt_o_tt - tf_c_d2t_o_tt + ...
    tt_o_te_c_d2t - crossMat(tilde_d2));
D2kappa1DfDe = D2kappa1DeDf';

% Curvature terms
DDkappa1(1:2, 1:2)  =   D2kappa1De2(1:2, 1:2);
DDkappa1(1:2, 3:4)  = - D2kappa1De2(1:2, 1:2) + D2kappa1DeDf(1:2, 1:2);
DDkappa1(1:2, 5:6) =               - D2kappa1DeDf(1:2, 1:2);
DDkappa1(3:4, 1:2)  = - D2kappa1De2(1:2, 1:2)                + D2kappa1DfDe(1:2, 1:2);
DDkappa1(3:4, 3:4)  =   D2kappa1De2(1:2, 1:2) - D2kappa1DeDf(1:2, 1:2) - ...
    D2kappa1DfDe(1:2, 1:2) + D2kappa1Df2(1:2, 1:2);
DDkappa1(3:4, 5:6) =                 D2kappa1DeDf(1:2, 1:2)                - D2kappa1Df2(1:2, 1:2);
DDkappa1(5:6, 1:2)  =                              - D2kappa1DfDe(1:2, 1:2);
DDkappa1(5:6, 3:4)  =                                D2kappa1DfDe(1:2, 1:2) - D2kappa1Df2(1:2, 1:2);
DDkappa1(5:6, 5:6) =                                               D2kappa1Df2(1:2, 1:2);

%% Hessian of Eb
dkappa = kappa1 - kappaBar;
dJ = 1.0 / l_k * EI * gradKappa * transpose(gradKappa);
temp = 1.0 / l_k * dkappa * EI;
dJ = dJ + temp * DDkappa1;

end

function A = crossMat(a)
A = [0, -a(3), a(2); ...
    a(3), 0, -a(1); ...
    -a(2), a(1), 0];
end