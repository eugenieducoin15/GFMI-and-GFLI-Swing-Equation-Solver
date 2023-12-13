function [d_gfl_str_pre_fault, w_gfl_pre_fault, d_gfl_pre_fault, num2_denom2_pre_fault, d_gfl_str_post_fault, w_gfl_post_fault, d_gfl_post_fault, num2_denom2_post_fault, B] = GFLdiffeq_ieee14_fct(n_tot, n_gfm, n_gfl, n_connected2GFM, X, Pl, DPl, dz, ki, kp, mp, Pn, wn, V)
if n_gfl ~= 0
    
    disp('looking for symbolic functions for d, w and d2ii')
    
    t = sym('t'); 
    iic = sym('iic',[1 n_tot]);
    dc = sym('dc',[1 n_tot]);
    
    assume(t, 'real');
    assumeAlso(t,'positive');
    assumeAlso(dc, 'real');
    assumeAlso(iic, 'real');
    
    for i = n_gfm+n_gfl+1:n_tot
        iic(i) = subs(iic(i),iic(i),0);
    end
    
    Y_pre_fault = 1i./X;
    Y_post_fault = Y_pre_fault;
    for i_de = 1:n_tot
        Y_pre_fault(i_de,i_de) = sum([Pl(i_de)/V(i_de)^2,-Y_pre_fault(i_de,1:i_de-1),-Y_pre_fault(i_de,i_de+1:end)]);
        Y_post_fault(i_de,i_de) = sum([(Pl(i_de)-DPl(i_de))/V(i_de)^2,-Y_post_fault(i_de,1:i_de-1),-Y_post_fault(i_de,i_de+1:end)]);
    end
    
    Yiv_pre_fault = Y_pre_fault(n_gfm+1:n_tot,1:n_gfm);
    Yii_inv_pre_fault = Y_pre_fault(n_gfm+1:n_tot,n_gfm+1:n_tot)^-1;
    Yiv_post_fault = Y_post_fault(n_gfm+1:n_tot,1:n_gfm);
    Yii_inv_post_fault = Y_post_fault(n_gfm+1:n_tot,n_gfm+1:n_tot)^-1;
    
    e_jd_pre_fault = -(diag(V(n_gfm+1:n_tot))-Yii_inv_pre_fault*diag(iic(n_gfm+1:end))*diag(cos(-dz(n_gfm+1:end))+1i*sin(-dz(n_gfm+1:end))))^-1*Yii_inv_pre_fault*Yiv_pre_fault*diag(V(1:n_gfm))*(cos(dc(1:n_gfm))+1i*sin(dc(1:n_gfm))).';
    e_jd_post_fault = -(diag(V(n_gfm+1:n_tot))-Yii_inv_post_fault*diag(iic(n_gfm+1:end))*diag(cos(-dz(n_gfm+1:end))+1i*sin(-dz(n_gfm+1:end))))^-1*Yii_inv_post_fault*Yiv_post_fault*diag(V(1:n_gfm))*(cos(dc(1:n_gfm))+1i*sin(dc(1:n_gfm))).';
    
    %only buses connected to GFM and GFL buses
    e_jd_pre_fault = simplify(e_jd_pre_fault(1:n_gfl+n_connected2GFM));
    tan_d_pre_fault = simplify(imag(e_jd_pre_fault)./real(e_jd_pre_fault));
    e_jd_post_fault = simplify(e_jd_post_fault(1:n_gfl+n_connected2GFM));
    tan_d_post_fault = simplify(imag(e_jd_post_fault)./real(e_jd_post_fault));
    disp('function for tan(d) obtained')
    
    %only GFL buses
    [num_pre_fault, denom_pre_fault] = numden(tan_d_pre_fault(1:n_gfl));
    [num_post_fault, denom_post_fault] = numden(tan_d_post_fault(1:n_gfl));
    
    %%
    terms = [iic(n_gfm+1:n_gfm+n_gfl), cos(dc(1:n_gfm)), sin(dc(1:n_gfm))];
    
    num_pre_fault = collect(num_pre_fault,terms);
    denom_pre_fault = collect(denom_pre_fault,terms);
    Cn_d_pre_fault = cell(n_gfm+n_gfl, 1);
    Tn_d_pre_fault = cell(n_gfm+n_gfl, 1);
    Cd_d_pre_fault = cell(n_gfm+n_gfl, 1);
    Td_d_pre_fault = cell(n_gfm+n_gfl, 1);
    for i_de = n_gfm+1:n_gfm+n_gfl
        [temp1,temp2] = coeffs(num_pre_fault(i_de-n_gfm),terms);
        [temp3,temp4] = coeffs(denom_pre_fault(i_de-n_gfm),terms);
        Cn_d_pre_fault(i_de,:) = {temp1};
        Tn_d_pre_fault(i_de,:) = {temp2};
        Cd_d_pre_fault(i_de,:) = {temp3};
        Td_d_pre_fault(i_de,:) = {temp4};
    end
    
    num_post_fault = collect(num_post_fault,terms);
    denom_post_fault = collect(denom_post_fault,terms);
    Cn_d_post_fault = cell(n_gfm+n_gfl, 1);
    Tn_d_post_fault = cell(n_gfm+n_gfl, 1);
    Cd_d_post_fault = cell(n_gfm+n_gfl, 1);
    Td_d_post_fault = cell(n_gfm+n_gfl, 1);
    for i_de = n_gfm+1:n_gfm+n_gfl
        [temp1,temp2] = coeffs(num_post_fault(i_de-n_gfm),terms);
        [temp3,temp4] = coeffs(denom_post_fault(i_de-n_gfm),terms);
        Cn_d_post_fault(i_de,:) = {temp1};
        Tn_d_post_fault(i_de,:) = {temp2};
        Cd_d_post_fault(i_de,:) = {temp3};
        Td_d_post_fault(i_de,:) = {temp4};
    end
    
    %%
    if isequal(Tn_d_pre_fault, Tn_d_post_fault)
        post_fault_calculated = true;
        disp('calculating both pre and post fault')
    else
        post_fault_calculated = false;
        disp('unable to calculate post fault')
    end
    
    cn = cell(n_gfm+n_gfl, 1);
    cd = cell(n_gfm+n_gfl, 1);
    for i_de = n_gfm+1:n_gfm+n_gfl
        cn(i_de,:) = {sym(['cn',num2str(i_de),'_'],[1 length(Tn_d_pre_fault{i_de,:})])}; 
        cd(i_de,:) = {sym(['cd',num2str(i_de),'_'],[1 length(Td_d_pre_fault{i_de,:})])}; 
    end
    
    d = getsymfund(n_gfm);
    ii = getsymfunii(n_gfm+n_gfl);
    
    old_terms = [dc(1:n_gfm) iic(n_gfm+1:n_gfm+n_gfl)];
    new_terms = [d(1:n_gfm) ii(n_gfm+1:n_gfm+n_gfl)];
    
    for i_de = n_gfm+1:n_gfm+n_gfl   
        num_tan_d(i_de) = sum(cn{i_de,:}.*subs(Tn_d_pre_fault{i_de,:}, old_terms, new_terms));
        denom_tan_d(i_de) = sum(cd{i_de,:}.*subs(Tn_d_pre_fault{i_de,:}, old_terms, new_terms));
    end
    %%
    % d = atan(u/v)
    % d'=(u'v-uv')/(u^2+v^2)
    % d''=((u''v-uv'')(u^2+v^2) - 2(uv(u'^2-v'^2) + u'v'(v^2-u^2)))/(u^2+v^2)^2;
    
    for i_de = n_gfm+1:n_gfm+n_gfl   
        %first derivative d
        num1(i_de) = diff(num_tan_d(i_de),t)*denom_tan_d(i_de)-num_tan_d(i_de)*diff(denom_tan_d(i_de),t); %(u'v-uv')
        denom1(i_de) = num_tan_d(i_de)^2 + denom_tan_d(i_de)^2; %(u^2+v^2)
        %second derivative d
        num2(i_de) = (diff(num_tan_d(i_de),t,2)*denom_tan_d(i_de) - num_tan_d(i_de)*diff(denom_tan_d(i_de),t,2))*(num_tan_d(i_de)^2 + denom_tan_d(i_de)^2) - 2*(num_tan_d(i_de)*denom_tan_d(i_de)*(diff(num_tan_d(i_de),t)^2-diff(denom_tan_d(i_de),t)^2) + diff(num_tan_d(i_de),t)*diff(denom_tan_d(i_de),t)*(denom_tan_d(i_de)^2-num_tan_d(i_de)^2)); %((u''v-uv'')(u^2+v^2) - 2(uv(u'^2-v'^2) + u'v'(v^2-u^2)))
        denom2(i_de) = (num_tan_d(i_de)^2 + denom_tan_d(i_de)^2)^2; %(u^2+v^2)^2
    end
    disp('function for w obtained')
    
    %%
    d2d = sym('d2d',[1 n_tot]);
    d2ii = sym('d2ii',[1 n_tot]);
    dd = sym('dd',[1 n_tot]);
    dii = sym('dii',[1 n_tot]);
    
    diffd = getsymfundiffd(n_gfm);
    diffii = getsymfundiffii(n_gfm+n_gfl);
    
    old_terms = [d(1:n_gfm) ii(n_gfm+1:n_gfm+n_gfl) diffd(1:n_gfm) diffii(n_gfm+1:n_gfm+n_gfl) diffd(n_gfm+1:2*n_gfm) diffii(2*n_gfm+n_gfl+1:2*(n_gfm+n_gfl))];
    new_terms = [dc(1:n_gfm) iic(n_gfm+1:n_gfm+n_gfl) dd(1:n_gfm) dii(n_gfm+1:n_gfm+n_gfl) d2d(1:n_gfm) d2ii(n_gfm+1:n_gfm+n_gfl)];
    
    for i_de = n_gfm+1:n_gfm+n_gfl
        num1(i_de) = subs(num1(i_de), old_terms, new_terms);
        num2(i_de) = subs(num2(i_de), old_terms, new_terms);
        denom1(i_de) = subs(denom1(i_de), old_terms, new_terms);
        denom2(i_de) = subs(denom2(i_de), old_terms, new_terms);
    end
    
    %%
    old_terms = [];
    new_terms_pre_fault = [];
    new_terms_post_fault = [];
    for i_de = n_gfm+1:n_gfm+n_gfl
        old_terms = [old_terms cn{i_de} cd{i_de}];
        new_terms_pre_fault = [new_terms_pre_fault Cn_d_pre_fault{i_de} Cd_d_pre_fault{i_de}];
        if post_fault_calculated
            new_terms_post_fault = [new_terms_post_fault Cn_d_post_fault{i_de} Cd_d_post_fault{i_de}]; 
        end
    end
    for i_de = n_gfm+1:n_gfm+n_gfl
        num2_pre_fault(i_de) = subs(num2(i_de), old_terms, new_terms_pre_fault);
        denom2_pre_fault(i_de) = subs(denom2(i_de), old_terms, new_terms_pre_fault);
        num2_denom2_pre_fault(i_de) = vpa(num2_pre_fault(i_de)/denom2_pre_fault(i_de));
        if post_fault_calculated
            num2_post_fault(i_de) = subs(num2(i_de), old_terms, new_terms_post_fault);
            denom2_post_fault(i_de) = subs(denom2(i_de), old_terms, new_terms_post_fault); 
            num2_denom2_post_fault(i_de) = vpa(num2_post_fault(i_de)/denom2_post_fault(i_de));
        end
    end
    
    B = ki(n_gfm+1:n_gfm+n_gfl)'./(V(n_gfm+1:n_gfm+n_gfl)'.*cos(dz(n_gfm+1:n_gfm+n_gfl)')).*(Pn(n_gfm+1:n_gfm+n_gfl)'+mp(n_gfm+1:n_gfm+n_gfl)'*wn);
    disp('matrices for d2ii obtained')
    
    %% w
    for i_de = n_gfm+1:n_gfm+n_gfl
        num1_pre_fault(i_de) = subs(num1(i_de), old_terms, new_terms_pre_fault);
        denom1_pre_fault(i_de) = subs(denom1(i_de), old_terms, new_terms_pre_fault);
        w_gfl_pre_fault(i_de-n_gfm,1) = num1_pre_fault(i_de)/denom1_pre_fault(i_de);
        if post_fault_calculated
            num1_post_fault(i_de) = subs(num1(i_de), old_terms, new_terms_pre_fault);
            denom1_post_fault(i_de) = subs(denom1(i_de), old_terms, new_terms_pre_fault);
            w_gfl_post_fault(i_de-n_gfm,1) = num1_post_fault(i_de)/denom1_post_fault(i_de);
        end
    end
    %%
    d_gfl_pre_fault = atan(tan_d_pre_fault);
    d_gfl_str_pre_fault = replace_str(d_gfl_pre_fault, false, n_gfm, n_gfl);
    if post_fault_calculated
        d_gfl_post_fault = atan(tan_d_post_fault);
        d_gfl_str_post_fault = replace_str(d_gfl_post_fault, false, n_gfm, n_gfl);
    else
        d_gfl_str_post_fault = d_gfl_str_pre_fault;
        w_gfl_post_fault = w_gfl_pre_fault;
        d_gfl_post_fault = d_gfl_pre_fault;
        num2_denom2_post_fault = num2_denom2_pre_fault;
    end
    
    disp('all functions obtained')
else
    d_gfl_str_pre_fault = 0;
    w_gfl_pre_fault = 0;
    d_gfl_pre_fault = 0;
    num2_denom2_pre_fault = 0; 
    d_gfl_str_post_fault = 0; 
    w_gfl_post_fault = 0;
    d_gfl_post_fault = 0; 
    num2_denom2_post_fault = 0; 
    B = 0;
end
    
    %%
    function out = getsymfund(n)
        d1=str2sym('d1(t)');
        d2=str2sym('d2(t)');
        d3=str2sym('d3(t)');
        d4=str2sym('d4(t)');
        d5=str2sym('d5(t)');
        d6=str2sym('d6(t)');
        d7=str2sym('d7(t)');
        d8=str2sym('d8(t)');
        d9=str2sym('d9(t)');
        d10=str2sym('d10(t)');
        d11=str2sym('d11(t)');
        d12=str2sym('d12(t)');
        d13=str2sym('d13(t)');
        d14=str2sym('d14(t)');
        temp_getsymfund = [d1 d2 d3 d4 d5 d6 d7 d8 d9 d10 d11 d12 d13 d14];
        out = temp_getsymfund(1:n);
    end
    
    function out = getsymfunii(n)
        ii1=str2sym('ii1(t)');
        ii2=str2sym('ii2(t)');
        ii3=str2sym('ii3(t)');
        ii4=str2sym('ii4(t)');
        ii5=str2sym('ii5(t)');
        ii6=str2sym('ii6(t)');
        ii7=str2sym('ii7(t)');
        ii8=str2sym('ii8(t)');
        ii9=str2sym('ii9(t)');
        ii10=str2sym('ii10(t)');
        ii11=str2sym('ii11(t)');
        ii12=str2sym('ii12(t)');
        ii13=str2sym('ii13(t)');
        ii14=str2sym('ii14(t)');
        temp_getsymfunii = [ii1 ii2 ii3 ii4 ii5 ii6 ii7 ii8 ii9 ii10 ii11 ii12 ii13 ii14];
        out = temp_getsymfunii(1:n);
    end

    function out = getsymfundiffd(n)
        d1=str2sym('d1(t)');
        d2=str2sym('d2(t)');
        d3=str2sym('d3(t)');
        d4=str2sym('d4(t)');
        d5=str2sym('d5(t)');
        d6=str2sym('d6(t)');
        d7=str2sym('d7(t)');
        d8=str2sym('d8(t)');
        d9=str2sym('d9(t)');
        d10=str2sym('d10(t)');
        d11=str2sym('d11(t)');
        d12=str2sym('d12(t)');
        d13=str2sym('d13(t)');
        d14=str2sym('d14(t)');
        temp1_getsymfundiffd = [diff(d1, t) diff(d2, t) diff(d3, t) diff(d4, t) diff(d5, t) diff(d6, t) diff(d7, t) diff(d8, t) diff(d9, t) diff(d10, t) diff(d11, t) diff(d12, t) diff(d13, t) diff(d14, t)];
        temp2_getsymfundiffd = [diff(d1, t, t) diff(d2, t, t) diff(d3, t, t) diff(d4, t, t) diff(d5, t, t) diff(d6, t, t) diff(d7, t, t) diff(d8, t, t) diff(d9, t, t) diff(d10, t, t) diff(d11, t, t) diff(d12, t, t) diff(d13, t, t) diff(d14, t, t)];
        out = [temp1_getsymfundiffd(1:n),temp2_getsymfundiffd(1:n)];
    end

    function out = getsymfundiffii(n)
        ii1=str2sym('ii1(t)');
        ii2=str2sym('ii2(t)');
        ii3=str2sym('ii3(t)');
        ii4=str2sym('ii4(t)');
        ii5=str2sym('ii5(t)');
        ii6=str2sym('ii6(t)');
        ii7=str2sym('ii7(t)');
        ii8=str2sym('ii8(t)');
        ii9=str2sym('ii9(t)');
        ii10=str2sym('ii10(t)');
        ii11=str2sym('ii11(t)');
        ii12=str2sym('ii12(t)');
        ii13=str2sym('ii13(t)');
        ii14=str2sym('ii14(t)');
        temp1_getsymfundiffii = [diff(ii1, t) diff(ii2, t) diff(ii3, t) diff(ii4, t) diff(ii5, t) diff(ii6, t) diff(ii7, t) diff(ii8, t) diff(ii9, t) diff(ii10, t) diff(ii11, t) diff(ii12, t) diff(ii13, t) diff(ii14, t)];
        temp2_getsymfundiffii = [diff(ii1, t, t) diff(ii2, t, t) diff(ii3, t, t) diff(ii4, t, t) diff(ii5, t, t) diff(ii6, t, t) diff(ii7, t, t) diff(ii8, t, t) diff(ii9, t, t) diff(ii10, t, t) diff(ii11, t, t) diff(ii12, t, t) diff(ii13, t, t) diff(ii14, t, t)];
        out = [temp1_getsymfundiffii(1:n),temp2_getsymfundiffii(1:n)];
    end
    
    function out = replace_str(in, unique, n_gfm, n_gfl)

        in_str = string(in);
        
        old_str = {'dc1','dc2','dc3','dc4','dc5','dc6','dc7','dc8','dc9','dc10','dc11','dc12','dc13','dc14'};
        if ~unique 
            new_str = {'d(1,:)','d(2,:)','d(3,:)','d(4,:)','d(5,:)','d(6,:)','d(7,:)','d(8,:)','d(9,:)','d(10,:)','d(11,:)','d(12,:)','d(13,:)','d(14,:)'};
        else
            new_str = {'d(1)','d(2)','d(3)','d(4)','d(5)','d(6)','d(7)','d(8)','d(9)','d(10)','d(11)','d(12)','d(13)','d(14)'};
        end
        out_str = replace(in_str, old_str(1:n_gfm), new_str(1:n_gfm));
        
        old_str = {'iic1','iic2','iic3','iic4','iic5','iic6','iic7','iic8','iic9','iic10','iic11','iic12','iic13','iic14'};
        if ~unique
            new_str = {'ii(1,:)','ii(2,:)','ii(3,:)','ii(4,:)','ii(5,:)','ii(6,:)','ii(7,:)','ii(8,:)','ii(9,:)','ii(10,:)','ii(11,:)','ii(12,:)','ii(13,:)','ii(14,:)'};
        else
            new_str = {'ii(1)','ii(2)','ii(3)','ii(4)','ii(5)','ii(6)','ii(7)','ii(8)','ii(9)','ii(10)','ii(11)','ii(12)','ii(13)','ii(14)'};
        end
        out_str = replace(out_str, old_str(n_gfm+1:n_gfm+n_gfl), new_str(n_gfm+1:n_gfm+n_gfl));
        
        old_str = {'dd1','dd2','dd3','dd4','dd5','dd6','dd7','dd8','dd9','dd10','dd11','dd12','dd13','dd14'};
        if ~unique
            new_str = {'dd(1,:)','dd(2,:)','dd(3,:)','dd(4,:)','dd(5,:)','dd(6,:)','dd(7,:)','dd(8,:)','dd(9,:)','dd(10,:)','dd(11,:)','dd(12,:)','dd(13,:)','dd(14,:)'};
        else
            new_str = {'dd(1)','dd(2)','dd(3)','dd(4)','dd(5)','dd(6)','dd(7)','dd(8)','dd(9)','dd(10)','dd(11)','dd(12)','dd(13)','dd(14)'};
        end
        out_str = replace(out_str, old_str(1:n_gfm), new_str(1:n_gfm));
        
        old_str = {'dii1','dii2','dii3','dii4','dii5','dii6','dii7','dii8','dii9','dii10','dii11','dii12','dii13','dii14'};
        if ~unique
            new_str = {'dii(1,:)','dii(2,:)','dii(3,:)','dii(4,:)','dii(5,:)','dii(6,:)','dii(7,:)','dii(8,:)','dii(9,:)','dii(10,:)''dii(11,:)','dii(12,:)','dii(13,:)','dii(14,:)'};
        else
            new_str = {'dii(1)','dii(2)','dii(3)','dii(4)','dii(5)','dii(6)','dii(7)','dii(8)','dii(9)','dii(10)','dii(11)','dii(12)','dii(13)','dii(14)'};
        end
        out_str = replace(out_str, old_str(n_gfm+1:n_gfm+n_gfl), new_str(n_gfm+1:n_gfm+n_gfl));
        
        old_str = {'d2d1','d2d2','d2d3','d2d4','d2d5','d2d6','d2d7','d2d8','d2d9','d2d10','d2d11','d2d12','d2d13','d2d14'};
        new_str = {'d2d(1)','d2d(2)','d2d(3)','d2d(4)','d2d(5)','d2d(6)','d2d(7)','d2d(8)','d2d(9)','d2d(10)','d2d(11)','d2d(12)','d2d(13)','d2d(14)'};
        out_str = replace(out_str, old_str(1:n_gfm), new_str(1:n_gfm));
        
        old_str = {'d2ii1','d2ii2','d2ii3','d2ii4','d2ii5','d2ii6','d2ii7','d2ii8','d2ii9','d2ii10','d2ii11','d2ii12','d2ii13','d2ii14'};
        new_str = {'d2ii(1)','d2ii(2)','d2ii(3)','d2ii(4)','d2ii(5)','d2ii(6)','d2ii(7)','d2ii(8)','d2ii(9)','d2ii(10)','d2ii(11)','d2ii(12)','d2ii(13)','d2ii(14)'};
        out_str = replace(out_str, old_str(n_gfm+1:n_gfm+n_gfl), new_str(n_gfm+1:n_gfm+n_gfl));
        
        if ~unique
            old_str = {'*','/','^'};
            new_str = {'.*','./','.^'};
            out_str = replace(out_str, old_str, new_str);
        end
        
        out = out_str;
        
    end
end