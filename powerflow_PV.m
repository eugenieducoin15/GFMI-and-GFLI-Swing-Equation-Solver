%Power Flow with PV buses
function [delta_final, S_final] = powerflow_PV(Xl, Pn, Pl, V_mag, n_tot, iterations)
y = -1i./Xl;
Y = -y;
for i_pf = 1:n_tot
    Y(i_pf,i_pf) = sum([y(i_pf,1:i_pf-1),y(i_pf,i_pf+1:end)]);
end
DP = Pn-Pl;

iterations = iterations+1;
V = zeros(iterations,n_tot);
V(1,:) = V_mag;
Vt = zeros(iterations,n_tot);
S = zeros(iterations,n_tot);
for iteration_pf = 2:iterations
    %bus 1
    V(iteration_pf,1) = V_mag(1)+0*1i;
    temp = Y(1,1)*V(iteration_pf,1);
    for j_pf = 2:n_tot
        temp = temp + Y(1,j_pf)*V(iteration_pf-1,j_pf);
    end      
    S(iteration_pf,1) = conj(conj(V(iteration_pf,1))*temp);
   
    %bus 2 to n_tot
    for i_pf = 2:n_tot
        temp = 0;
        for j_pf = 1:i_pf-1
            temp = temp + Y(i_pf,j_pf)*V(iteration_pf,j_pf);
        end
        for j_pf = i_pf:n_tot
            temp = temp + Y(i_pf,j_pf)*V(iteration_pf-1,j_pf);
        end
        S(iteration_pf,i_pf) = DP(i_pf) + 1i*(-imag(conj(V(iteration_pf-1,i_pf))*temp));

        temp = 0;
        for j_pf = 1:i_pf-1
            temp = temp + Y(i_pf,j_pf)*V(iteration_pf,j_pf);
        end
        for j_pf = i_pf+1:n_tot
            temp = temp + Y(i_pf,j_pf)*V(iteration_pf-1,j_pf);
        end
        Vt(iteration_pf,i_pf) = (conj(S(iteration_pf,i_pf))/conj(V(iteration_pf-1,i_pf))-temp)/Y(i_pf,i_pf);

        [tempx,tempy] = pol2cart(angle(Vt(iteration_pf,i_pf)),V_mag(i_pf));
        V(iteration_pf,i_pf) = complex(tempx,tempy);
    end    
end

delta_final = angle(V(iterations,:));   
S_final = S(iterations,:)+Pl;

%%
% results = zeros(iterations,7);
% results(:,1)=1:iterations;
% results(:,2)=abs(V(iteration_pf,1));
% results(:,3)=abs(V(iteration_pf,2));
% results(:,4)=abs(V(iteration_pf,3));
% results(:,5)=abs(V(iteration_pf,4));
% results(:,6)=angle(V(iteration_pf,1));
% results(:,7)=angle(V(iteration_pf,2));
% results(:,8)=angle(V(iteration_pf,3));
% results(:,9)=angle(V(iteration_pf,4));
% 
% P12 = results(length(results),2)*results(length(results),3)/Xl(1,2)*sin(results(length(results),6)-results(length(results),7));
% P13 = results(length(results),2)*results(length(results),4)/Xl(1,3)*sin(results(length(results),6)-results(length(results),8));
% P14 = results(length(results),2)*results(length(results),5)/Xl(1,4)*sin(results(length(results),6)-results(length(results),9));
% P23 = results(length(results),3)*results(length(results),4)/Xl(2,3)*sin(results(length(results),7)-results(length(results),8));
% P24 = results(length(results),3)*results(length(results),5)/Xl(2,4)*sin(results(length(results),7)-results(length(results),9));
% P34 = results(length(results),4)*results(length(results),5)/Xl(3,4)*sin(results(length(results),8)-results(length(results),9));
% 
% DP(2,1) = P12+P13+P14;
% DP(2,2) = -P12+P23+P24;
% DP(2,3) = -P13-P23+P34;
% DP(2,4) = -P14-P24-P34;



