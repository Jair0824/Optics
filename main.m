%% --------------------------------------------------------------------------%%
%% 清除空间
clear;
clc;
%% 夫琅禾夫衍射--圆孔
%%  三维--解析法
a=8.5*10^-2;%半径
lambda=632*10^(-9);%波长
z=a^(2)/lambda*100;%距离
N=16;%空间范围
L=512;%采样点
[x, y]=meshgrid(linspace(-N,N,L));%网格
p=sqrt(x.^2+y.^2);%每个点直接的距离
m=2*pi*a*p/z/lambda;%自变量
j1=besselj(1,m);%贝塞尔
i_0=(2*pi*a^(2)/(2*z*lambda))^2;%初始光强
i=i_0*(2*j1./m).^2;%透射光强
figure()%画图
mesh(x,y,i);%三维
zlabel('强度')
xlabel('X方向')
ylabel('Y方向')
title('圆孔-夫琅禾费衍射')
grid on;%上色
colormap;%显色
colorbar;%颜色棒
%%  三维--快速傅里叶变换
figure()
z1=zeros(L,L);%透射空间
for i =1:L%在透射空间中画出圆形，圆形内部透射率为1
    for j =1:L
        if x(i,j)^2+y(i,j)^2<=a%判断距离
            z1(i,j)=1;%透射率为1
        end 
    end
end
d=fftshift(abs((fft2(z1)).^2));%快速傅里叶变换，并移动坐标为原点
mesh(x,y,d);%三维
zlabel('强度')
xlabel('X方向')
ylabel('Y方向')
title('圆孔-夫琅禾费衍射')
grid on;%上色
colormap;%显色
colorbar;%颜色棒
%% 下面参数基本相同不做标记
%% 二维--解析法
figure()
p_1=-10:0.0001:10;
m_1=2*pi*a*p_1/z/lambda;
j1=besselj(1,m_1);
i_1=(2*j1./m_1).^2;
plot(m_1,i_1)
xlabel('离几何中心距离')
ylabel('强度')
title('圆孔-夫琅禾费衍射')
%% 二维--快速傅里叶变换
figure()
p_1=linspace(-8,8,L);
d=fftshift(abs((fft2(z1)).^2));
plot(p_1,d(L/2,:))
xlabel('离几何中心距离')
ylabel('强度')
title('圆孔-夫琅禾费衍射')
%% 菲涅尔衍射-圆孔
%%  三维--快速傅里叶变换
figure()
z=a^(2)*pi/4/lambda*100;
z1=zeros(L,L);
z2=zeros(L,L);
for i =1:L
    for j =1:L
        if x(i,j)^2+y(i,j)^2<=a
            z1(i,j)=1;
            z2(i,j)=x(i,j)^2+y(i,j)^2;
        end 
    end
end
d=fftshift(abs((fft2(z1.*exp(1i*2*pi/lambda/2/z*z2))).^2));
i_2=d;
mesh(x,y,i_2);grid on;
zlabel('强度')
xlabel('X方向')
ylabel('Y方向')
title('圆孔-菲涅尔衍射')
colormap;
colorbar;
%% 二维--快速傅里叶变换
figure()
p_1=linspace(-8,8,L);
d=fftshift(abs((fft2(z1.*exp(1i*2*pi/lambda/2/z*z2))).^2));
plot(p_1,d(L/2,:))
xlabel('离几何中心距离')
ylabel('强度')
title('圆孔-菲涅尔衍射')
%%
%% 
%解析法模拟计算圆形孔的菲涅耳衍射图样
clc;
clear;
close all;
warning off;
addpath(genpath(pwd));

step = 350;
lamda = 500e-6; %changed
k = 2*pi/lamda;
z = 12.5; %changed 
%确定衍射屏
N = 500; %圆屏采样点数
r = 0.25; %changed
I = zeros(N, N);
[m, n] = meshgrid(linspace(-N/step, N/step, N));
D = (m.^2+n.^2).^(1/2);
i = find(D <= r);
I(i) = 1;  %空半径范围内透射系数为1
q = exp(j*k*(m.^2+n.^2)/2/z);
subplot(2,2,1); %圆孔图像
imshow(I);
%imagesc(I) %衍射屏图像
%colormap([0 0 0;1 1 1]) %黑白区分
% 
% I = I.*q;
L = 500;
M = 500; %取相同点数用于矩阵运算
[x, y] = meshgrid(linspace(-L/step, L/step, M));
h = exp(j*k*z)*exp((j*k*(x.^2+y.^2))/(2*z))/(j*lamda*z); %接收屏
%H = fftshift(fft2(h));
B = fftshift(fft2(I.*q));
G = h.*B; %
% U = fftshift(ifft2(G));
%Br = (abs(G)/max(abs(G))); %归一化
C = abs(G);
subplot(2,2,2);imagesc(C);
axis image;
colormap(hot);
% %figure;
subplot(2,2,3);mesh(x,y,abs(G));
subplot(2,2,4);
axis image;
d = C(251,:);
d = d/max(d);
plot(d);

%% 夫琅禾夫衍射--方形孔
%%  三维--解析法
a=8.5*10^-2;
b=8.5*10^-2;
lambda=573.3*10^(-9);
N=32;L=1024;
z=a^(2)/lambda*100;
[x,y]=meshgrid(linspace(-N,N,L));
p1=x.*a/(lambda*z);p2=y.*b/(lambda*z);
i_0=(a*b/(lambda*z))^2;
i=i_0*(sinc(p1).^2).*(sinc(p2).^2);
figure()
mesh(x,y,i);grid on;
zlabel('强度')
xlabel('X方向')
ylabel('Y方向')
title('方形孔-夫琅禾费衍射')
colormap;
colorbar;
%% 菲涅尔衍射-方形孔
%%  三维--快速傅里叶变换
a=2.25*10^-1;
b=2.25*10^-1;
lambda=573*10^(-9);
N=32;L=1024;
z=a^(2)*pi/4/lambda*500;
[x,y]=meshgrid(linspace(-N,N,L));
z1=zeros(L,L);
z2=zeros(L,L);
for i =1:L
    for j =1:L
        if abs(x(i,j))<=a && abs(y(i,j))<=b
            z1(i,j)=1;
            z2(i,j)=x(i,j)^2+y(i,j)^2;
        end 
    end
end
d=fftshift(abs((fft2(z1.*exp(1i*2*pi/lambda/2/z*z2))).^2));
figure()
mesh(x,y,d);grid on;
zlabel('强度')
xlabel('X方向')
ylabel('Y方向')
title('方形孔-菲涅尔衍射')
colormap;
colorbar;
%% --------------------------------------------------------------------------%%
%%菲涅尔矩孔衍射--解析法
%解析法模拟计算矩形孔的菲涅耳衍射图样Rectangular_hole_Fresnel_diffraction  
clear;
clc;
close all;

%定义变量
%波长
lambda = 0.000632; %mm
%衍射屏到成像屏的距离d
d = 200; %mm
%取样数
N = 3000;
%波数
k = 2*pi/lambda;
%像屏宽度
wide = 8 ;%mm
%矩形孔半宽度
omegax = 3 ;%mm
omegay = 2 ;%mm
%衍射距离d
d = 500; %mm

%让用户输入衍射距离
d = input('请输入衍射距离d=?(mm)');


%利用循环进行求解,先对x进行循环
for p = 1:N
x = -wide/2 + wide*(p - 1)/N;%像屏面的x轴坐标的大小
alpha1 = sqrt(2/lambda/d)*(omegay + x);%这个为什么是oemgay而不是omegax呢，还没完全想明白。
alpha2 = sqrt(2/lambda/d)*(x - omegay);

S1 = S(alpha1);
S2 = S(alpha2);
C1 = C(alpha1);
C2 = C(alpha2);
Ix(p)=((S2-S1)^2+(C2-C1)^2)/2;

end


%利用循环进行求解,对y进行循环
for q = 1:N
y = -wide/2 + wide*(q - 1)/N;%像屏面的x轴坐标的大小
beta1 =  sqrt(2/lambda/d)*(omegax + y);
beta2 =  sqrt(2/lambda/d)*(y - omegax);

S1 = S(beta1);
S2 = S(beta2);
C1 = C(beta1);
C2 = C(beta2);
Iy(q)=((S2-S1)^2+(C2-C1)^2)/2;
end

%求得强度I
I=zeros(N,N);
for p=1:N
for q=1:N
I(p,q)=Ix(p)*Iy(q);
end
end

%x，y轴强度曲线显示
figurestrX=strcat('X方向宽度=',num2str(2*omegax),'mm'); 
x=1:N;
figure(1), plot(x,Ix);title('X方向轴向曲线');
xlabel(figurestrX);
figurestrY=strcat('Y方向宽度=',num2str(2*omegay),'mm'); 
y=1:N;
figure(2), plot(y,Iy);title('Y方向轴向曲线');
xlabel(figurestrY);

%IMAX限幅显示
% figure(3), imshow(I);
figstr=strcat('矩形孔X向宽=',num2str(2*omegax),'mm , Y向宽=',num2str(2*omegay),'mm，观测面宽=',num2str(wide),'mm，衍射距离',num2str(d),'mm');
Imax=max(max(I));
f=uint8(I./Imax*255);
figure, imshow(I,[0 Imax]);
title('矩形孔菲涅尔衍射图像');
xlabel(figstr)

%将图像缩小然后建立三维图。
%若不进行缩小则三维由于取样太密集变成黑色
g=0.02;
I1=imresize(I,g);
M=g*N;
n=1:M;
x=-wide/2+wide/M*(n-1);                     
y=x;
figure(4),surf(x,y,I1),colorbar;title('矩形孔衍射场强度分布形貌');

