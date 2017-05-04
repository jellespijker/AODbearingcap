%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%IHC MTI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculation script Dredge Crawler Archimedes screw torque prediction

% This script is made to predict the required torque to initiate movement 
% through a certain type of soil. The script is inspired by the excelsheet 
% "Propulsion" created by Rick Lotman in 2009.

% See for an thorough explaination:

% 64115R03 Rev B Applicable theory for the modelling of the propulsion system of the SMT

% Bear in mind that it is NOT an exact copy! Some methodologies used in the
% model presented underneath are thought to be more precise. 

% Script made by R.G. van de Wetering, tel: +31 (6) 24 38 06 37
% Created at 10/04/2017

tic 
clc
clear all
close all

%% User Input 
disp('----------------------------------------------------------------------------')
disp('Choose a certain soil type if the properties are not known, else type "user"')
disp('----------------------------------------------------------------------------')

% User input 
soil = input('Soil Type? (silt/clay/sand):','s')

    if strcmpi(soil,'user') == 1 
        c       = input('Cohesion? [Pa] :')
        rho_ins = input('In-situ bottom density? [kg/m^3] :')
        phi_deg = input('Internal friction angle? [deg] :')
        K0      = input('Coefficient of lateral earth pressure at rest? [-] :')
    elseif strcmpi(soil,'silt') == 1
        c       = 3e3;                      % Cohesion                      [Pa]
        rho_ins = 1300;                     % In-situ bottom density        [kg/m^3]
        phi_deg = 0;                        % Internal friction angle       [deg]
        K0      = 0.47;                     % Coeff of lateral earth press  [-]
    elseif strcmpi(soil,'clay') == 1  
        c       = 5e3;                      % Cohesion                      [Pa]
        rho_ins = 1400;                     % In-situ bottom density        [kg/m^3]
        phi_deg = 0;                        % Internal friction angle       [deg]
        K0      = 0.84;                     % Coeff of lateral earth press  [-]                     
    elseif strcmpi(soil,'sand') == 1
        c       = 0;                        % Cohesion                      [Pa]
        rho_ins = 2000;                     % In-situ bottom density        [kg/m^3]
        phi_deg = 35;                       % Internal friction angle       [deg]
        phi     = phi_deg*(pi/180);         % Internal friction angle       [rad]
        K0      = 1-sin(phi);               % Coeff of lateral earth press  [-]
    else
        disp('--------------------------------------------------------------------------')
        disp('Current input cannot be used! Please re-run script and enter a valid input')
        disp('--------------------------------------------------------------------------')
        return
    end 
    
%% Soil & water parameters

%c           = 3e3;                          % Cohesion                      [Pa]
%phi_deg     = 35;                           % Internal friction angle       [deg]
%rho_ins     = 1300;                         % In-situ bottom density        [kg/m^3]
%K0          = 0.5;                          % Coeff of lateral earth press  [-]    
rho_s       = 2650;                         % Soil density                  [kg/m^3]
rho_air     = 1.225;                        % Air density at sealevel T=15  [kg/m^3]
phi         = phi_deg*(pi/180);             % Internal friction angle       [rad]
delta       = tan(phi);                     % Steel-soil friction factor    [-]
T_f         = 15;                           % Water temperature             [C]
g           = 9.81;                         % Gravitational constant        [m^2/s]
                  
% Water density dependency on temperature (validity 5<T_f<100C) (Matousek, 2004) [kg/m^3]
rho_f = 999.7-0.10512*(T_f-10)-0.005121*(T_f-10)^2+0.00001329*(T_f-10)^3;

% Water dynamic & kinematic viscosity as function of temperature (Matousek, 2004) [Pa s]
mu_f  = 0.10/(2.1482*((T_f-8.435)+sqrt(8078.4+(T_f-8.435)^2))-120);
nu_f  = mu_f/rho_f;

% Porosity                                                                  [-]
n           = (rho_s - rho_ins)/(rho_s-rho_f);

S_s         = rho_s/rho_f;                  % Specific gravity/ rel. density[-]
gamma       = g*(rho_ins-rho_f);            % Submerged soil weight         [N/m^3] 

%% Archimedes screw parameters

% Screw geometry 
no_screw    = 2;                            % Number of screws              [-]
D_screw     = 0.6;                          % Cylinder diameter             [m]
r_screw     = D_screw/2;                    % Cylinder radius               [m]
L_screw     = 1.92;                         % Cylinder lenght               [m]
O_screw     = pi*D_screw;                   % Cylinder circumference        [m]     
% Archimedes screw volume                                                   [m^3]
V_screw     = (1/4)*pi*(D_screw^2)*L_screw;  

% Total buoyancy obtained due to both archimedes screws                     [N]
F_buoy_f    = V_screw*2*(rho_f-rho_air)*g;

% Helix geometry 
no_helix    = 1;                           % Number of parallel helices     [-]
p_helix     = 0.478;                       % Pitch helix                    [m]
h_vane      = 0.1;                         % Vane height                    [m]
D_helix     = D_screw+(2*h_vane);          % Outer diameter helix           [m]
O_helix     = pi*D_helix;                  % Outer circumference helix      [m]
D_helix_mid = D_screw+(2*0.5*h_vane);      % Mid helix diameter             [m]
r_helix_mid = D_helix_mid/2;               % Mid helix radius               [m]
O_helix_mid = pi*D_helix_mid;              % Mid helix circumference        [m]
alpha       = atan(p_helix/O_screw);       % Lead angle cylinder surface    [rad]
alpha_mid   = atan(p_helix/O_helix_mid);   % Lead angle mid helix           [rad] 
no_thread   = no_helix*(L_screw/p_helix);  % Number of threads              [-]
L_helix_mid = no_thread*sqrt((p_helix^2)...
    +(D_helix_mid*pi)^2);                  % Mid total helix length         [m]
A_helix     = L_helix_mid*h_vane;          % Total surface area             [m^2]

%% Crawler parameters 

W_dry       = 4100;                        % Dry weight                     [kg]
F_dry       = W_dry*g;                     % Dry weight grav. force         [N]
F_wet       = (F_dry - F_buoy_f)/2;        % Submerged weight per screw     [N]
V           = 0.4;                         % Velocity                       [m/s]
V_rot       = V/p_helix;                   % Rotational velocity helix      [1/s]
V_rot2      = V_rot*2*pi;                  % Rotational velocity helix      [rad/s]

%% Bearing Capacity 

% Dimensionless constants in the Brinch-Hansen model (verruijt, 2009)
N_q         = (1+sin(phi))/(1-sin(phi))*exp(pi*tan(phi)); 
N_gamma     = 2*(N_q-1)*tan(phi);

    if phi == 0
        N_c = 2+pi;
    else
        N_c = (N_q-1)*cot(phi);   
    end 

% Initial values 
D(1)            = 0.49*D_screw;                             % Initial condition depth           [m]
q(1)            = gamma*D(1);                               % Submerged Soil weight             [N/m^2]

% Width "foundation" as function of sinkage depth Brinch Hansen model                           [-] 
B_screw(1)      = sqrt((r_screw^2)-(r_screw-D(1))^2);      

theta(1)        = asin(B_screw(1)/r_screw);                 % Sink-angle                        [rad]

V_soil(1)       = 0.5*r_screw^2*(2*theta(1)...              % Volume of screw submerged in bottom [m^3]
    -sin(2*theta(1)))*L_screw;
F_buoy_soil(1)  = V_soil(1)*(rho_ins-rho_f)*g;              % Extra soil buoyancy force         [N]

S_c(1)          = 1+0.2*((2*B_screw(1))/L_screw);           % Shape factor Brinch Hansen        [-]
S_q(1)          = 1+((2*B_screw(1))/L_screw)*sin(phi);      % Shape factor Brinch Hansen        [-]
S_gamma(1)      = 1-0.3*((2*B_screw(1))/L_screw);           % Shape factor Brinch Hansen        [-]

% Allowed load according to Brinch Hansen                                                       [N/m^2]
p(1)            = (S_c(1)*c*N_c)+(S_q(1)*q(1)*N_q)+(S_gamma(1)*0.5*gamma*2*B_screw(1)*N_gamma);

A_req(1)        = (F_wet - F_buoy_soil(1))/p(1);            % required surface area for equilibrium [m^2]

% Required foundation width "B" as function of constant length L_screw                          [m]
B_screw_req(1)  = A_req(1)/L_screw;

% Convergence criterium                                                                         [-]
epsilon         = 1E-6;                                     

% Iteratively solving the dependency between p,q and D using the Brinch Hansen model
for i = 1:1000
    % Calculate next "sink-angle"                                                               [rad]
    theta(i+1) = asin((0.5*B_screw_req(i))/r_screw);
    
        % Model is only valid if the sink angle remains between 0-90 degrees (D = r_screw)
        if rad2deg(theta(i+1)) >= 90 
            disp('----------------------------------------------------')
            disp('WARNING, sink angle too large! Calculation has ended')
            disp('----------------------------------------------------')
                return
        elseif rad2deg(theta(i+1)) <= 0
            disp('---------------------------------------------------------')
            disp('WARNING, no sinkage seems to occur! Calculation has ended')
            disp('---------------------------------------------------------')
                return
        end 
    
    % Calculate next depth as function of required foundation width B_screw_req: [m]
    D(i+1) = r_screw - r_screw*cos(theta(i+1));
    
    % Calculate next submerged soil weight overburden                       [N/m^2]
    q(i+1) = gamma*D(i+1);
    
    % Calculate next shape factors 
    S_c(i+1)     = 1+0.2*(B_screw_req(i)/L_screw);                          %[-]
    S_q(i+1)     = 1+(B_screw_req(i)/L_screw)*sin(phi);                     %[-]
    S_gamma(i+1) = 1-0.3*(B_screw_req(i)/L_screw);                          %[-]
    
    % Calculate next maximum allowed load according to Brinch Hansen 
    p(i+1)       = (S_c(i+1)*c*N_c) + (S_q(i+1)*q(i+1)*N_q) + (S_gamma(i+1)*0.5*gamma*B_screw_req(i)*N_gamma);
    
    % Calculate next screw volume in soil                                   [m^3]
    V_soil(i+1)  = 0.5*r_screw^2*(2*theta(i+1)-sin(2*theta(i+1)))*L_screw;
    
    % Calculate next extra buoyancy force due to sinkage                    [N]
    F_buoy_soil(i+1) = V_soil(i+1)*(rho_ins-rho_f)*g;
    
    % Calculate next required surface area for equilibrium                  [m^2]
    A_req(i+1)   = (F_wet - F_buoy_soil(i+1))/p(i+1); 
    
    % Calculate next required foundation width "B_screw_req"                [m]
    B_screw_req(i+1) = A_req(i+1)/L_screw; 
    
        % Convergence criterium                                             [-]
        if abs(D(i+1) - D(i)) < epsilon
            break
        end
end

%% Soil Friction Forces on Archimedes screw

% Maximum vertical TOTAL soil stress using the bearing capacity equilibrium [N/m^2]
sigma_v       = p(end);

% Maximum horizontal TOTAL soil stress assuming soil is at rest             [N/m^2]
sigma_h       = sigma_v*K0;

% Total skin friction force on cylinder                                     [N]
F_fric_cyl    = (F_dry - F_)*delta;  

% Skin friction force cylinder in rotational direction                      [N]
F_fric_cyl_x  = F_fric_cyl*cos(alpha);  

% Resulting torque to overcome cylinder friction                            [Nm]
T_fric_cyl    = F_fric_cyl*r_screw; 

% Skin friction force cylinder in thrust direction                          [N]
F_fric_cyl_y  = F_fric_cyl*sin(alpha);

% Determine resulting angle in sinkage triangle with theta:                 [rad]
beta          = pi-0.5*pi-theta;

% Determine length side B                                                   [m]
B             = h_vane*tan(beta(end));

% Frictional area vane in soil as function of depth (D)                     [m^2]
A_helix_soil  = (2*no_thread*(theta(end)*(((r_screw+h_vane)^2)-...
        (r_screw^2)) + (B*h_vane)))/cos(alpha_mid); 
    
% Skin friction force  
F_fric_helix  = sigma_h*A_helix_soil*delta;

% Skin friction force helix in thrust direction                             [N]
F_fric_helix_x = F_fric_helix*cos(alpha_mid);

% Resulting torque to overcome helix friction                               [Nm]
T_fric_helix   = F_fric_helix_x*(r_screw+0.5*h_vane);

% Skin friction force helix in rotational direction                         [N]
F_fric_helix_y = F_fric_helix*sin(alpha_mid); 

% Total torque required to overcome friction in rotational direction        [Nm]
T_fric_tot     = T_fric_cyl + T_fric_helix;  

%% Propulsion forces due to torque measurement

% From the performed torque measurements the maximum deliverd torque is
% known. From this torque the propulsion forces can be calculated, knowing
% the geometry of the Archimedes screw

% Delivered torque per screw                                                [Nm]                              
T_del   = 1770;                                                              

    if T_del <= T_fric_tot || T_fric_tot == 0
        disp('--------------------------------------------------------------')
        disp('Not enough torque available to overcome frictional resistance!')
        disp('--------------------------------------------------------------')
    else 
        disp('-------------------------------------------------------')
        disp('Torque is sufficient to overcome frictional resistance!')
        disp('-------------------------------------------------------')
    end 
    
% Calculate the force in rotational direction due to delivered torque       [N]
F_del_x     = T_del/(r_screw+0.5*h_vane);

% Calculate force in thrust direction due to delivered torque               [N]
F_del_y     = F_del_x/tan(alpha_mid);

% Calculate normal force due to delivered torque                            [N]
F_del_n     = F_del_x/sin(alpha_mid);

%% Mohr Coulomb failure analysis - Rankine method

%     % Coefficient of passive soil pressure in Rankine's method (verruijt, 2009) [-]
%     if phi == 0
%         Kp = 1;
%     else 
%         Kp = (1+sin(phi))/(1-sin(phi));
%     end 
% 
% % Maximum allowable horizontal soil stress (Rankine's method)               [N/m^2]
% sigma_h_allow = (Kp*sigma_v)+(2*c*sqrt(Kp));  
% 
%     if sigma_h_allow >= sigma_h 
%         disp('Shear strengh soil is sufficient to initiate relative motion!')
%     else 
%         disp('Soil failure occurs, no relative motion is found!')
%     end 

toc 