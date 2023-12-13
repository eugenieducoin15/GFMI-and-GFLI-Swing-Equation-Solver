tc = 10;%2; %time of step change
t_end = 20;%2.1; %end time
w_nom = 2*pi*50; %nominal frequency
n_tot = 14; %total number of buses

Pset = [1.0081 0.40 0.30 0 0 0.30 0 0.30 0 0 0.100 0.100 0.100 0]; %Power setpoint (0 if no generation)
Pl = [0 0.217 0.942 0.478 0.076 0.112 0 0 0.295 0.09 0.035 0.061 0.135 0.149]; %load
DPl = -[0 0 0 0 0 0 0 0 0 0.1 0 0 0 0]; %Pnom - Pl: change in load power (negative for additional load) 

%line reactance (inf in diagonal, 0 in top right triangle, inf if buses are
%not connected
Xl = [inf,             0,         0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0;
      0.05917,       inf,         0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0;
      inf,       0.19797,       inf,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0;
      inf,       0.17632,   0.17103,      inf,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0;
      0.22304,   0.17388,       inf,  0.04211,      inf,        0,        0,        0,        0,        0,        0,        0,        0,        0;
      inf,           inf,       inf,      inf,  0.25202,      inf,        0,        0,        0,        0,        0,        0,        0,        0;
      inf,           inf,       inf,  0.20912,      inf,      inf,      inf,        0,        0,        0,        0,        0,        0,        0;
      inf,           inf,       inf,      inf,      inf,      inf,  0.17615,      inf,        0,        0,        0,        0,        0,        0;
      inf,           inf,       inf,  0.55618,      inf,      inf,  0.11001,      inf,      inf,        0,        0,        0,        0,        0;
      inf,           inf,       inf,      inf,      inf,      inf,      inf,      inf,  0.08450,      inf,        0,        0,        0,        0;
      inf,           inf,       inf,      inf,      inf,  0.19890,      inf,      inf,      inf,  0.19207,      inf,        0,        0,        0;
      inf,           inf,       inf,      inf,      inf,  0.25581,      inf,      inf,      inf,      inf,      inf,      inf,        0,        0;
      inf,           inf,       inf,      inf,      inf,  0.13027,      inf,      inf,      inf,      inf,      inf,  0.19988,      inf,        0;
      inf,           inf,       inf,      inf,      inf,      inf,      inf,      inf,  0.27038,      inf,      inf,      inf,  0.34802,      inf];
Xl = Xl+Xl.'; % make Xl matrix symmetric

V = [1.06 1.045  1.01 1.029 1.031 1.07 1.065 1.09 1.058 1.053 1.062 1.065 1.056 1.039]; %voltage magnitude

n_gfm = 5; %number of GFM with droop
gfm_buses = [1 2 3 6 8]; %buses at which gfm with droop are connected
H = [2.66 4.99 1.52 1.2 1.2]; %in (pu of power).s or GVA.s
Dgfm = [100 20 15 15 15]; %in pu of power

n_gfl = 3; %number of GFL with droop
gfl_buses = [11 12 13]; %buses at which gfl with droop are connected
w_pll = 2*pi*[0.25 0.25 0.25]; %PLL bandwidth in rad/s. Order of values should match gfl_buses
zeta = [1 1 1]; %PLL damping in rad/s. Order of values should match gfl_buses
Dgfl = [2 2 2]; %in pu of power