function d_gfl = getd_ieee14(bus, n_gfm, n_gfl, x, d_str, d_sym)
for i = 1:n_gfm
    d(i,:)=x(i+n_gfm+n_gfl,:);
end
for i = n_gfm+1:n_gfm+n_gfl
    ii(i,:)=x(i+n_gfm+n_gfl,:);
end
    d_gfl = eval(d_str(bus-n_gfm));
    if isnan(d_gfl)
        dc1 = sym('dc1');
        dc2 = sym('dc2');
        dc3 = sym('dc3');
        dc4 = sym('dc4');
        dc5 = sym('dc5');
        dc6 = sym('dc6');
        dc7 = sym('dc7');
        dc8 = sym('dc8');
        dc9 = sym('dc9');
        dc10 = sym('dc10');
        dc11 = sym('dc11');
        dc12 = sym('dc12');
        dc13 = sym('dc13');
        dc14 = sym('dc14');
        temp_dc = [dc1 dc2 dc3 dc4 dc5 dc6 dc7 dc8 dc9 dc10 dc11 dc12 dc13 dc14];
        iic1 = sym('iic1');
        iic2 = sym('iic2');
        iic3 = sym('iic3');
        iic4 = sym('iic4');
        iic5 = sym('iic5');
        iic6 = sym('iic6');
        iic7 = sym('iic7');
        iic8 = sym('iic8');
        iic9 = sym('iic9');
        iic10 = sym('iic10');
        iic11 = sym('iic11');
        iic12 = sym('iic12');
        iic13 = sym('iic13');
        iic14 = sym('iic14');
        temp_iic = [iic1 iic2 iic3 iic4 iic5 iic6 iic7 iic8 iic9 iic10 iic11 iic12 iic13 iic14];

        d_gfl = vpa(subs(d_sym(bus-n_gfm),[temp_dc(1:n_gfm) temp_iic(n_gfm+1:n_gfm+n_gfl)],x(n_gfm+n_gfl+1:end-1)'));
    end
end