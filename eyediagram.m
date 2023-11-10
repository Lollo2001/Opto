t0=5e-9;
T_bit=10;
T_w=linspace(t0-T_bit/2,t0+T_bit*3/2, 10^4);

figure()
for i=1:64
    itime=find(time>=t0-T_bit/2+2*T_bit*(i-1) & time<=t0+T_bit*3/2+2*T_bit*(i-1));
    Pf1=Py(itime);
    hold on
    plot(itime,Pf1)
end