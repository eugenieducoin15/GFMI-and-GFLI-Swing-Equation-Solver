%function [w_out, wcoi_out, Pg_out] = swing_equation_14bus(tc, t_end, w_nom, n_tot, Pn, Pl, DPl, Xl, V, n_gfm, gfm_buses, Hin, Dgfm, n_gfl, gfl_buses, Dgfl, w_pll, zeta)
% user_main;
%%
Hin = H;
H = zeros(1,n_tot);
D = zeros(1,n_tot);
mp = zeros(1,n_tot);
H(gfm_buses) = 2*Hin/w_nom; %in this script, Hscript = 2H
D(gfm_buses) = Dgfm/w_nom;
mp(gfm_buses) = 1./D(gfm_buses);

kp = zeros(1,n_tot);
kp(gfl_buses) = 2*zeta.*w_pll;
ki = zeros(1,n_tot);
ki(gfl_buses) = w_pll.^2;
D(gfl_buses) = Dgfl/w_nom;
mp(gfl_buses) = D(gfl_buses); % in this script, for GFL, we have DP = mp Dw whereas in thesis, DP = 1/mp Dw.
H(gfl_buses) = mp(gfl_buses).*kp(gfl_buses)./ki(gfl_buses);

gfm_buses = sort(gfm_buses);
gfl_buses = sort(gfl_buses);
connected2GFM_buses = isConnected2GFM(Xl, n_tot, gfm_buses, gfl_buses);
n_connected2GFM = length(connected2GFM_buses);
all_buses = 1:n_tot;
other_buses = nonzeros(not(ismember(all_buses,[gfm_buses gfl_buses connected2GFM_buses])).*all_buses)';

Xl = rearrange_buses_mat(Xl, gfm_buses, gfl_buses, connected2GFM_buses, other_buses);
Pn = rearrange_buses_vec(Pset, gfm_buses, gfl_buses, connected2GFM_buses, other_buses);
Pl = rearrange_buses_vec(Pl, gfm_buses, gfl_buses, connected2GFM_buses, other_buses);
DPl = rearrange_buses_vec(DPl, gfm_buses, gfl_buses, connected2GFM_buses, other_buses);
V = rearrange_buses_vec(V, gfm_buses, gfl_buses, connected2GFM_buses, other_buses);
kp = rearrange_buses_vec(kp, gfm_buses, gfl_buses, connected2GFM_buses, other_buses);
ki = rearrange_buses_vec(ki, gfm_buses, gfl_buses, connected2GFM_buses, other_buses);
H = rearrange_buses_vec(H, gfm_buses, gfl_buses, connected2GFM_buses, other_buses);
D = rearrange_buses_vec(D, gfm_buses, gfl_buses, connected2GFM_buses, other_buses);
mp = rearrange_buses_vec(mp, gfm_buses, gfl_buses, connected2GFM_buses, other_buses);

kp = kp(1:n_gfm+n_gfl);
ki = ki(1:n_gfm+n_gfl);
H = H(1:n_gfm+n_gfl);
D = D(1:n_gfm+n_gfl);
mp = mp(1:n_gfm+n_gfl);

DtGFM = sum(D(1:n_gfm));
HtGFM = sum(H(1:n_gfm));

Dt = sum(D);
Ht = sum(H);
Pnt = sum(Pn);
Plt = sum(Pl);
DPlt = sum(DPl);

%initial values
[d0, S0] = powerflow_PV(Xl, Pn, Pl, V, n_tot, 1e3);
w0 = ones(1,n_tot)*w_nom;
si0 = zeros(1,n_tot);    
dz = angle(S0);
i0 = real(S0)./(V.*cos(dz));
initCond = [w0(1:n_gfm), si0(n_gfm+1:n_gfm+n_gfl), d0(1:n_gfm), i0(n_gfm+1:n_gfm+n_gfl), w_nom]; %wgfm sigfl dgfm igfl wgfm

%get str and sym formulas for d w and d2i
[d_str_pre_fault, w_sym_pre_fault, d_sym_pre_fault, num2_denom2_pre_fault, d_str_post_fault, w_sym_post_fault, d_sym_post_fault, num2_denom2_post_fault, B] = GFLdiffeq_ieee14_fct(n_tot, n_gfm, n_gfl, n_connected2GFM, Xl, Pl, DPl, dz, ki, kp, mp, Pn, w_nom, V);

%%
%solve differential equations
disp('solving equations')

varPassedOut = zeros(n_gfl+1,1);
opts = odeset('MaxStep',1e-3); 
%     opts = odeset('MaxStep',1e-1);
out = ode15s(@(t, x) sys_step_noCOI(t, x, n_gfm, n_gfl, n_connected2GFM, n_tot, V, Xl, w_nom, H, D, Pn, mp, ki, kp, Ht, Pnt, Dt, dz, Pl, DPl, tc, Plt, DPlt, d_str_pre_fault, w_sym_pre_fault, d_sym_pre_fault, num2_denom2_pre_fault, d_str_post_fault, w_sym_post_fault, d_sym_post_fault, num2_denom2_post_fault, B), [0 t_end], initCond, opts);
disp('equations solved')
%%
disp('processing data')
t_graph = out.x;
w = zeros(n_tot, length(t_graph));
d = zeros(n_tot, length(t_graph));
Pg = zeros(n_tot, length(t_graph));
ii = zeros(n_tot, length(t_graph));
dPg = zeros(n_tot, length(t_graph));
dii = zeros(n_tot, length(t_graph));

for i = 1:n_gfm
    w(i,:) = out.y(i,:);
    d(i,:) = out.y(n_gfl+n_gfm+i,:);
end
for i = n_gfm+1:n_gfm+n_gfl
    dii(i,:) = out.y(i,:);
    dPg(i,:) = V(i)*dii(i,:)*cos(dz(i));
    ii(i,:) = out.y(n_gfl+n_gfm+i,:);
    Pg(i,:) = V(i)*ii(i,:)*cos(dz(i)); 
end

Rl = zeros(n_tot, length(t_graph));
for i = n_gfm+1:n_gfm+n_gfl
    Rl(i,:) = V(i)^2./getPl(t_graph, tc, Pl(i), DPl(i));
end

[~,tc_idx] = min(abs(t_graph-tc));
if t_graph(tc_idx)-tc < 0
    [~,temp] = min(abs(t_graph(tc_idx:end)-tc));
    tc_idx = temp + tc_idx;
end

for i = n_gfm+1:n_gfm+n_gfl+n_connected2GFM
    d(i,1:tc_idx-1) = getd_ieee14(i, n_gfm, n_gfl, [zeros(n_gfm+n_gfl,tc_idx-1); d(1:n_gfm,1:tc_idx-1); ii(n_gfm+1:n_gfm+n_gfl,1:tc_idx-1)], d_str_pre_fault, d_sym_pre_fault); 
    d(i,tc_idx:end) = getd_ieee14(i, n_gfm, n_gfl, [zeros(n_gfm+n_gfl,length(t_graph)-(tc_idx-1)); d(1:n_gfm,tc_idx:end); ii(n_gfm+1:n_gfm+n_gfl,tc_idx:end)], d_str_post_fault, d_sym_post_fault); 
    for j = 1:length(t_graph)
        while abs(d(i,j)-d(1,j))>pi/2
            d(i,j)=d(i,j)+pi;
        end
    end
end

done = false;
i = 1;
while ~done
    if i < length(varPassedOut) && i <= length(t_graph)
        if varPassedOut(1,i) ~= t_graph(i)
            varPassedOut = [varPassedOut(:,1:i-1) varPassedOut(:,i+1:end)];
        else
            i = i + 1;
        end
    elseif i <= length(varPassedOut) && i > length(t_graph)
        if varPassedOut(1,1:i-1) == t_graph
            varPassedOut = varPassedOut(:,1:i-1);
            done = true;
        end
    else
        done = true;
    end
end
w(n_gfm+1:n_gfm+n_gfl,:) = varPassedOut(2:end,:);

for i=1:n_gfm
    Ptemp = 0;
    for j=1:n_tot
        if Xl(i,j) ~= inf
            Ptemp = Ptemp + V(i)*V(j)/Xl(i,j).*sin(d(i,:)-d(j,:));
        end
    end
    Pg(i,:) = getPl(t_graph, tc, Pl(i), DPl(i))+Ptemp;
end

wcoi_out = sum((w(1:8,:).*H'),1)/Ht;

disp('data processed')

w_out = zeros(n_tot, length(t_graph));
Pg_out = zeros(n_tot, length(t_graph));

w_out(gfm_buses,:) = w(1:n_gfm,:);
w_out(gfl_buses,:) = w(1+n_gfm:n_gfm+n_gfl,:);  
Pg_out(gfm_buses,:) = Pg(1:n_gfm,:);
Pg_out(gfl_buses,:) = Pg(1+n_gfm:n_gfm+n_gfl,:); 


%%
function  dx = sys_step_noCOI(t, x, n_gfm, n_gfl, n_connected2GFM, n_tot, V, Xl, w_nom, H, D, Pn, mp, ki, kp, Ht, Pnt, Dt, dz, Pl, DPl, tc, Plt, DPlt, d_str_pre_fault, w_sym_pre_fault, d_sym_pre_fault, num2_denom2_pre_fault, d_str_post_fault, w_sym_post_fault, d_sym_post_fault, num2_denom2_post_fault, B)
    
    dx = zeros(2*(n_gfm+n_gfl)+1,1); %w dii delta ii

%     if t>15
%         disp('t>15s')
%     end
    if t < tc
        d_str = d_str_pre_fault;
        d_sym = d_sym_pre_fault;
        w_sym = w_sym_pre_fault;
        num2_denom2 = num2_denom2_pre_fault;
    else
        d_str = d_str_post_fault;
        d_sym = d_sym_post_fault;
        w_sym = w_sym_post_fault;
        num2_denom2 = num2_denom2_post_fault;
    end
    
    %voltage angle delta and frequency
    delta = zeros(1,n_gfm+n_gfl+n_connected2GFM); 
    w = zeros(1,n_gfm+n_gfl); 
    %GFM delta
    Rl = zeros(1, n_tot);
    for i = 1:n_gfm
        Rl(i) = V(i)^2./getPl(t, tc, Pl(i), DPl(i));
        delta(i) = x(n_gfm+n_gfl+i);
        w(i) = x(i);
    end
    %GFL & connected2GFM delta    
    for i = n_gfm+1:n_tot
        Rl(i) = V(i)^2./getPl(t, tc, Pl(i), DPl(i));
    end
    for i = n_gfm+1:n_gfm+n_gfl+n_connected2GFM
        delta(i) = getd_ieee14(i, n_gfm, n_gfl, x, d_str, d_sym); 
        while abs(delta(i)-delta(1))>pi/2
            delta(i)=delta(i)+pi;
        end
    end
    
    %generated real power Pg
    C = zeros(n_tot);
    for i = 1:n_tot
        for j = 1:n_tot
            C(i,j) = V(i)*V(j)/Xl(i,j);
        end
    end
    Pg = zeros(1, n_gfm+n_gfl);
    for i = n_gfm+1:n_gfm+n_gfl
        Pg(i) = V(i)*x(n_gfm+n_gfl+i)*cos(dz(i));
    end
    Psum = zeros(1,n_gfm);
    for i = 1:n_gfm
        for j = 1:n_gfm+n_gfl+n_connected2GFM
            if C(i,j) ~= 0
                Psum(i) = Psum(i) + C(i,j)*sin(delta(i)-delta(j));
            end
        end
        Pg(i) = getPl(t, tc, Pl(i), DPl(i)) + Psum(i);
    end
%     if n_gfm ~= 0        
%         Pg(n_gfm) = getPl(t, tc, Plt, DPlt) - sum(Pg);
%     end
    
    %GFM d/dt w = swing
    for i = 1:n_gfm     
        dx(i) = -D(i)/H(i)*(x(i)-w_nom)+(Pn(i)-Pg(i))/H(i);
    end    
    
    % GFM d/dt delta = w 
    for i = n_gfm+n_gfl+1:n_gfm+n_gfl+n_gfm
        dx(i) = x(i-(n_gfm+n_gfl));
    end 
    
    %GFL d/dt i = si
    for i = n_gfm+n_gfl+n_gfm+1:2*(n_gfm+n_gfl)
        dx(i) = x(i-(n_gfm+n_gfl));
    end  
    
    %GFL w
    varToPassOut = [t; zeros(n_gfl,1)];
    for i = n_gfm+1:n_gfm+n_gfl    
        w(i) = getw_ieee14(i, n_gfm, n_gfl, x, w_sym); 
        varToPassOut(i+1-n_gfm) = w(i);
    end      
    assignin('base','varInBase',varToPassOut);
    evalin('base','varPassedOut(:,end+1)=[varInBase];');
    
    %GFL d/dt si = swing
    if n_gfl ~= 0
        dx(n_gfm+1:n_gfm+n_gfl) = getd2ii_ieee14(n_gfm, n_gfl, x, dx, B, num2_denom2, mp, kp, ki, dz, w, V);
    end         
    
    %COI d/dt wcoi
    dx(2*(n_gfm+n_gfl)+1) = 1/Ht*(-sum(D.*w) + Dt*w_nom + Pnt - getPl(t, tc, Plt, DPlt) - sum(V(n_gfm+1:(n_gfm+n_gfl)).*cos(dz(n_gfm+1:(n_gfm+n_gfl))) ...
        ./ki(n_gfm+1:(n_gfm+n_gfl)).*(dx(n_gfm+1:(n_gfm+n_gfl))' + kp(n_gfm+1:(n_gfm+n_gfl)).*x(n_gfm+1:(n_gfm+n_gfl))'))); 
    
    if t==-1
        disp(['t = ', num2str(t)])
        for i = 1:2*(n_gfm+n_gfl)+1
            disp(['x', num2str(i), ' = ', num2str(x(i))])
        end
        for i = 1:2*(n_gfm+n_gfl)+1
            disp(['dx', num2str(i), ' = ', num2str(dx(i))])
        end
    end

end


function y = getPl(t, tc, Pl_init, Pl_change)
    y = Pl_init - heaviside(t-tc).*Pl_change;
end

function X_new= rearrange_buses_mat(X_old, gfm_buses, gfl_buses, connected2GFM_buses, other_buses)
    idx_old = [gfm_buses gfl_buses connected2GFM_buses other_buses];
    X_new = X_old(idx_old,idx_old);
end

function R_new= rearrange_buses_vec(R_old, gfm_buses, gfl_buses, connected2GFM_buses, other_buses)
    idx_old = [gfm_buses gfl_buses connected2GFM_buses other_buses];
    R_new = R_old(idx_old);
end

function idx_delta = isConnected2GFM(X, n_tot, gfm_buses, gfl_buses)
        % connected to GFM
    idx_delta = [];
    all_buses = 1:n_tot;
    for i = gfm_buses   
        idx_delta = [idx_delta nonzeros(not(isinf(X(i,:))).*all_buses)'];
    end
    idx_delta = unique(idx_delta);
    % remove GFM & GFL buses
    idx_delta = nonzeros(not(ismember(idx_delta,[gfm_buses gfl_buses])).*idx_delta)';
end