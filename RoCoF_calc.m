function [RoCoF_out,RoCoF_coi] = RoCoF_calc(n_tot,n_gfm,n_gfl,Xl,DPl,gfm_buses,gfl_buses,Hgfm_in,Dgfl,V,w_pll_in,zeta,w_nom,Pn,Pl)

H = zeros(1,n_tot);
H(gfm_buses) = 2*Hgfm_in/w_nom;

w_pll = zeros(1,n_tot);
w_pll(gfl_buses) = w_pll_in;
kp = zeros(1,n_tot);
kp(gfl_buses) = 2*zeta.*w_pll(gfl_buses);
ki = zeros(1,n_tot);
ki(gfl_buses) = w_pll(gfl_buses).^2;
H(gfl_buses) = Dgfl/w_nom.*kp(gfl_buses)./ki(gfl_buses);

gfm_buses = sort(gfm_buses);
gfl_buses = sort(gfl_buses);
connected2GFM_buses = isConnected2GFM(Xl, n_tot, gfm_buses, gfl_buses);
n_connected2GFM = length(connected2GFM_buses);
all_buses = 1:n_tot;
other_buses = nonzeros(not(ismember(all_buses,[gfm_buses gfl_buses connected2GFM_buses])).*all_buses)';

Xl = rearrange_buses_mat(Xl, gfm_buses, gfl_buses, connected2GFM_buses, other_buses);
DPl = rearrange_buses_vec(DPl, gfm_buses, gfl_buses, connected2GFM_buses, other_buses);
V = rearrange_buses_vec(V, gfm_buses, gfl_buses, connected2GFM_buses, other_buses);
H = rearrange_buses_vec(H, gfm_buses, gfl_buses, connected2GFM_buses, other_buses);
Pn = rearrange_buses_vec(Pn, gfm_buses, gfl_buses, connected2GFM_buses, other_buses);
Pl = rearrange_buses_vec(Pl, gfm_buses, gfl_buses, connected2GFM_buses, other_buses);
w_pll = rearrange_buses_vec(w_pll, gfm_buses, gfl_buses, connected2GFM_buses, other_buses);

% H = H(1:n_gfm+n_gfl);
HtGFM = sum(H(1:n_gfm));

[d0, ~] = powerflow_PV(Xl, Pn, Pl, V, n_tot, 1e3);
  
iterations = 1000;

type = [ones(1,n_gfm) 2*ones(1,n_tot-n_gfm)]; % 1: slack, 2: PV

[P_step, delta_step] = power_flow_for_RoCoF(Xl, -DPl, V, n_tot, iterations, type, d0);
RoCoF = zeros(1,n_tot);
RoCoF(1:n_gfm) = -P_step(1:n_gfm)./H(1:n_gfm)/(2*pi);

A_matx = diag(-ones(1,n_tot));
B_matx = zeros(n_tot,1);
for i = n_gfm+1:n_tot
    if i > n_gfm && i <= n_gfm+length(gfl_buses)
        denom = H(i)*w_pll(i)^2 + V(i)*sum(V./Xl(i,:));
    else
        denom = 0 + V(i)*sum(V./Xl(i,:));
    end
    
    for j = 1:n_gfm
        B_matx(i) = B_matx(i) - V(i)*V(j)/Xl(i,j)/denom*RoCoF(j);
    end
    for j = n_gfm+1:n_tot
        if j ~= i
            A_matx(i, j) = V(i)*V(j)/Xl(i,j)/denom;
        end
    end
end
sol = linsolve(A_matx,B_matx); %Ax=B

RoCoF = RoCoF + sol';
% RoCoF_coi = -sum(DPl)/HtGFM/(2*pi);
RoCoF_coi = sum(H.*RoCoF)/sum(H);
RoCoF_out([gfm_buses,gfl_buses,connected2GFM_buses,other_buses]) = RoCoF;

end
%%
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