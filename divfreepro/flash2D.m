%显示解随时间变化规律
h = figure();				% 创建图形窗口
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');	% 关闭相关的警告提示（因为调用了非公开接口）
jFrame = get(h,'JavaFrame');	% 获取底层 Java 结构相关句柄
pause(0.1);					% 在 Win 10，Matlab 2017b 环境下不加停顿会报 Java 底层错误。各人根据需要可以进行实验验证
set(jFrame,'Maximized',1);	%设置其最大化为真（0 为假）
pause(0.1);					% 个人实践中发现如果不停顿，窗口可能来不及变化，所获取的窗口大小还是原来的尺寸。各人根据需要可以进行实验验证
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');		% 打开相关警告设置

xlabel('X');
ylabel('Y');
ylabel('Q');

TT = 100; % 帧数
t0 = T(end)/TT; % 间隔

xc = zeros(Nx,Ny1);
y = zeros(Nx,Ny1);

for j = 1:Ny1
    xc(:,j) = Xc;
end

for i = 1:Nx
    y(i,:) = Y';
end

for i = 1:TT + 1
    tt = (i - 1)*t0;
    [~,j] = min(abs(T - tt));
    p = pcolor(xc,y,uflash(:,:,j)); colormap(jet); p.EdgeColor = 'none';
%    mesh(xc,y,uflash(:,:,j)); colormap(jet);axis([Xc(1),Xc(end),Y(1),Y(end),min(min(min(uflash(:,:,1)))) - 0.1,max(max(max(uflash(:,:,1)))) + 0.1]);
   pause(0.001);
end