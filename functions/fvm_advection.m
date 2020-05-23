%solve the advection equation with the finite volume method, using the van
%leer flux limiter and dimensional splitting

%inputs: Q0: column vector containing the average value in each cell
%        dt: timestep
%         h: 1xd array with cell widths in each dimension
%      ubar: 1xd cell array. each cell contains an array with the average velocity at the
%            cell interfaces
%      UP,DQI,DQi,DF: matrices produced by build_matrices_for_fvm_advection
%         T: time horizon

%output: Q: 1 x floor(T/dt) cell array, each cell contains an array with the average
%function values at each timestep

%created by Sean Vaskov on May 11th, 2020


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

