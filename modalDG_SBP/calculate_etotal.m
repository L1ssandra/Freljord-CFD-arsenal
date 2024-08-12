% calculate_etotal

v0 = vhpro(uh0);
v1 = vhpro(uh1mod);
v2 = vhpro(uh2mod);

etotal = 0;

for i = 1:Nx
    for d = 1:dimPk
    
        etotal = etotal + (v0(i,d)*du0(i,d) + 4*v1(i,d)*du1(i,d) + v2(i,d)*du2(i,d))*mm(d)*dt*hx/6;
        
    end
end