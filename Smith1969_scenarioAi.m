clc;
clear;
close all;

% For scenario A(i) constatn eddy viscosity and no-slip BC
Pc = 94000;   %(940 hpa-central pressure)   
Pg = 100000;  %(1000 hpa- geostrophic pressure)   
Rg = 1000000; %(1000 km- geostrophic radius)
Rm = 40000;   %(40 km- radius of max wind)
rho = 1.2;    %density of air
f = 5e-5;      
Km = 50;  % for scenario (i)

Zg=sqrt(Km/f);             % Ekman length scale

%Gradient wind profile (generlally like Holland 1980's profile)
b=Rm/Rg;
m=((Pg-Pc)/rho)/(Rg^2 * f^2);
% we need to solve for constant x that satisfies condition 1<x<2 from
% equation (31) of Smith 1968 and then use it on equation (29)
f_solve= @(x) m*x*((x-1)^2)*exp(x*(b-1))-((2-x)*(b^2));
x_val=fzero(f_solve,[1,2]); %nonlinear solver
fprintf('Calculated x value: %.3f\n', x_val);

%Radial profile
R=linspace(50,1000000,10000);
r=R/Rg;
v_tilde=(-0.5*r)+sqrt((0.25*(r.^2))+(m*x_val*b./r).*exp(x_val*b*(1-(1./r))));
Ro=-0.5+sqrt(0.25+(m*x_val*b)); % is 0.078 which is descirbed in page 481 of Smith (1968)
v_gr=v_tilde/Ro;
V_g=Ro*Rg*f; %page 476 after equatiuon (10)
V_gr=V_g*v_gr;

% we begin at r=1 (the geostrophic radius, 1000 km) and integrate inward
% initial Conditions: At r=1, the flow is purely geostrophic, so E=1 and delta=sqrt(2)​
% for A(i) constant eddy viscosity K_m and no-slip boundry condition
I1=1/8; I2=5/8; I3=1/2; I4=-3/8; I5=-1/2; fp0=-1; gp0=1;
A=I2/I1;
B=(I5-(2*I4))/(I5-I4);
C=fp0/I1;
D=gp0/(I5-I4);
X=I3/I1;
Y=I5/(I4-I5);
k=1;

dr_step=1e-3;
v_t_1=(-0.5)+sqrt((0.25)+(m*x_val*b/1));  %same v_tilde
r_plus = 1 + dr_step;
v_t_p = (-0.5*r_plus) + sqrt((0.25*(r_plus^2)) + (m*x_val*b/r_plus)*exp(x_val*b*(1-(1/r_plus))));
r_minus = 1 - dr_step;
v_t_m = (-0.5*r_minus) + sqrt((0.25*(r_minus^2)) + (m*x_val*b/r_minus)*exp(x_val*b*(1-(1/r_minus))));
dvgr_dr=(v_t_p-v_t_m)/(2*dr_step*Ro);

% Find E and delta such that dE/dr = 0 and d_delta/dr = 0 using Eq 21 & 22
% Initial guess: E=1, delta=sqrt(2)
% Setting dE/dr = 0 and d_delta/dr = 0 
% y(2) = E*delta^2 (Eq 22)
y2_start = -(k*gp0) / (Ro*((2 + 2*dvgr_dr)*I4 - (1 + dvgr_dr)*I5) + I5);
% y(1) = E^2 (Eq 21)
y1_start = -(y2_start * (Ro*I2 + I3)) / (Ro*y2_start*(1 +(2*dvgr_dr))*I1 + k*fp0);

%define the sysmte of ODEs (equations 24 and 25)
% y(1)=E^2; y(2)=E*(delta^2)
%we integrate from r=1 to r=0.03 (near the center)

r_span=[1 0.03]; 
initial_conds=[y1_start;y2_start]; %E=1, delta=sqrt(2) >> E^2=1, E*(delta^2)=2

ode_system=@(r,y) smith_odes(r,y,m,x_val,b,A,B,C,D,X,Y,k);

options=odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[r_out,y_out]=ode45(ode_system,r_span,initial_conds, options);

E_sq=y_out(:,1);
E_delta_sq=y_out(:,2);
E=sqrt(E_sq);
delta=sqrt(E_delta_sq./E);


% figures
figure (1);
plot(R/1000,V_gr,'k-','LineWidth',2);
xlabel('R(km)');
ylabel('V_{gr} (m/s)');
title('Figure 2: Zonal velocty at the top of BL');

figure (2);
plot(r_out * (Rg/1000), delta, 'k', 'LineWidth', 1.5);
set(gca, 'XDir');
xlabel('R(km)');
ylabel('\delta');
title('Figure 3: Non-dimensional scale BL thickness (\delta)');


v_t_out = (-0.5 * r_out) + sqrt((0.25 * (r_out.^2)) + (m * x_val * b ./ r_out) .* exp(x_val * b * (1 - (1 ./ r_out))));
v_gr_out=v_t_out / Ro;
V_gr_out=V_g*v_gr_out;
eta=pi/4;
fmax= exp(-eta)* sin(eta);
EV_gr=E.*V_gr_out*fmax; 

figure(3);
plot(r_out * (Rg/1000), EV_gr, 'k-', 'LineWidth', 2);
set(gca, 'XDir');
xlabel('R (km)');
ylabel('EV_{gr} (m/s)');
title('Figure 4: scale inflow velocity');

dr_val = 0.1;
r_eval = r_out; 
r_p = r_eval + dr_val;
vt_r = (-0.5*r_eval) + sqrt((0.25*(r_eval.^2)) + (m*x_val*b./r_eval).*exp(x_val*b*(1-(1./r_eval))));
vt_p = (-0.5*r_p) + sqrt((0.25*(r_p.^2)) + (m*x_val*b./r_p).*exp(x_val*b*(1-(1./r_p))));
term2_deriv = ((r_p.*vt_p) - (r_eval.*vt_r)) / dr_val;
w_non = -(I5 ./ (Ro .* delta)) .* ( D + E .* (delta.^2) .* ( ((1-B)./r_out) .* term2_deriv - Y ) );
W_cms=(V_g*Zg/Rg)*w_non*100;

figure(4);
plot(r_out * (Rg/1000), W_cms, 'k-', 'LineWidth', 2);
set(gca, 'XDir');
xlabel('R (km)');
ylabel('W_{gr} (cm/sec)');
title('Figure 5: Vertical velocity');

function dydr = smith_odes(r, y, m, x_val, b, A, B, C, D, X, Y, k)
    E_sq=y(1);
    E_delta_sq=y(2);
    E=sqrt(E_sq);

    exp_term = exp(x_val * b * (1 - 1/r));
    G = 0.25 * r^2 + (m * x_val * b / r) * exp_term;
    v_t = -0.5 * r + sqrt(G);
    dGdr = 0.5 * r + (m * x_val * b / r^2) * exp_term * (x_val * b / r - 1);
    dv_tdr = -0.5 + dGdr / (2 * sqrt(G));
    term1_deriv = v_t^2 + 2 * r * v_t * dv_tdr;
    term2_deriv = v_t + r * dv_tdr;

    % equation (24)
    dE2dr=E_sq*((-2/(r*(v_t^2)))*(term1_deriv - B*v_t*term2_deriv)-((1/E_sq)*((2*A)/r+(2*X)/v_t))+((2*Y)/v_t)-(2*(C + D)/(v_t*E_delta_sq)));
    % equation (25)
    dEdelta2dr = E_delta_sq * ((1/(r*v_t))*(term1_deriv - 3*B*v_t*term2_deriv) + (A/r + X/v_t)*(1/E_sq) - (3*Y)/v_t + (C + 3*D)/(v_t*E_delta_sq) );
    dydr=[dE2dr; dEdelta2dr];
end












