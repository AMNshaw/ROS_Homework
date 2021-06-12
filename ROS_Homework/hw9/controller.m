function out = controller(u,P)

% input(25*1):desired trajectory and full state feedback, x v R Omega time
% output(4*1): force and moment control input

%parameters
m = P.mass;
g = P.gravity;
k_x = P.kx;
k_v = P.kv;
k_R = P.kR;
k_Omega = P.kOmega;
J = [P.Jxx 0 0;...
    0 P.Jyy 0;...
    0 0 P.Jzz];
e_3 = [0; 0; 1];

% process inputs
R_d=eye(3);
Omega_d=[0;0.1;1];
x=Omega_d(1);
y=Omega_d(2);
z=Omega_d(3);
omega_hat=[0 -z y;...
           z 0 -x;...
           -y x 0];
delta_t = 0.05; 

x_d = u(1:3);
b1_d = u(4:6);
v_d = [0; 0; 0];
a_d = [0; 0; 0];
% current state
x  = u(7:9);
v  = u(10:12);
R  = reshape(u(13:21),3,3);
Omega = u(22:24);
t = u(end);
%error
e_x = x(3) - x_d(3);
e_v = v(3) - v_d(3);
e_R = 1/2 * vee(R_d'*R - R'*R_d);
e_Omega = Omega - R'*R_d*Omega_d;
% Desired total thrust & Moment
f = (k_x*e_x + k_v*e_v + m*g - m*a_d(3))/dot(e_3, R*e_3)  ; 
M = -k_R*e_R - k_Omega*e_Omega + cross(Omega, J*Omega);

out = [f;M;e_R];
end