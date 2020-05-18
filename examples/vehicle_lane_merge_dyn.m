function dzdt = vehicle_ss_dyn(~,z,p)

ksteer = [0.01,0.3];

y = z(2);
psi =  z(3);
vx  =  z(4);
 
vx_nom =  p(1);
y_nom =  p(2);

 m = 1558;
 lf = 1.462884;
 lr = 1.405516;
 l = lf + lr;
 
 Cf = 1.432394487827058e+05;
 Cr = 2.214094963126969e+05;
 g = 9.80655;
 
 
 a = -2*atan(vx-vx_nom);
 delta = -ksteer(1)*(y-y_nom)-ksteer(2)*(psi);
 
   
 kus = (m*g*lr/(l*Cf)-m*g*lf/(l*Cr));
 
 yr = delta*vx/(l+kus*vx^2/g);
 
 vy = yr*(lr-m*vx^2*lr/(Cr*l));

    dzdt = [vx*cos(psi)-vy*sin(psi);...
         vx*sin(psi)+vy*cos(psi);...
         yr;...
         a]; 

     
end