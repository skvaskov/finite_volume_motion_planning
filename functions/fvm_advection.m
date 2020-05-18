function Q = fvm_advection(Q0,dt,h,ubar,UP,DQI,DQi,DF,T)
d = length(ubar);


Q = cell(1,floor(T/dt));
Q{1} = Q0;

for n = 1:(floor(T/dt)-1)
 
    Qstar = Q{n};
    
    %sweep through each dimension
for j = 1:d    
    %compute flux at each edge of dimension j
    dQi = DQi{j}*Qstar;
    dQI = DQI{j}*Qstar;
    theta = dQI./(dQi+3e-14);
    phi =(theta+abs(theta))./(1+abs(theta));
    delta = phi.*dQi;

    F = UP{j}*Qstar + 1/2*abs(ubar{j}).*(1-abs(ubar{j}*dt/h(j))).*delta;
    
   Qstar = Qstar-dt/h(j)*DF{j}*F;

end

Q{n+1} = Qstar;
end

end

