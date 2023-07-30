%draw_solution.m
uhG = zeros(Nx,NumGLP);
for i = 1:Nx
    for d = 1:dimPk
        uhG(i,:) = uhG(i,:) + uh(i,d)*phiG(:,d)';
    end
end
uhG = reshape(uhG',Nx*NumGLP,1)';

hold on
plot(Xq,uhG,'rx-');
%plot(Xc,uh(:,1),'r*'); 
axis([Xq(1),Xq(end),min(uhG) - 0.1,max(uhG) + 0.1])