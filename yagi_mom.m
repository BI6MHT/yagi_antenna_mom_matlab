clear
freq = 300*10^6%工作频率
w = 2*pi*freq
freq0 = 300*10^6 %中心频率
lambda0 = 3*10^8/freq0;
lambda = 3*10^8/freq;
u = 4*pi/lambda;
epsilong = 8.854*(10^(-12));%真空电容率，常量
k = 2*pi/lambda %波数
l1 = 0.39 %反射器长度
l2 = 0.356 %主振子长度
l3 = 0.3 %引向器1长度
l4 = 0.3;%引向器2长度

n1 = 10; %分段数
n2 = 10;
n3 = 10;
n4 = 10;

step = 3;
%所有分段
NN = n1+n2+n3+n4;

% 假如一个振子被分成9段，那么它需要记录的点数是8个
% 默认两端点电流为0，只需计算中间那些点的电流
% 故可得N=(n1-1)+(n2-1)+...=NN-7
% 7单元的八木
N = NN - 7;

d1 = 0.265;%振子位置坐标
d2 = d1+0.2575;
d3 = d2+0.225;

delta1 = l1/n1;
delta2 = l2/n2;
delta3 = l3/n3;
delta4 = l4/n4;

a = lambda0*0.005;%振子半径

% 此处的作用是什么
x = zeros(N);
y = zeros(N);
z = zeros(N);
%原点，中点，终点
y = struct('orig',[],'mid',[],'over',[]);
z = struct('orig',[],'mid',[],'over',[]);

% 反射器的坐标
% y(1)代表所有字段的第一个值
% y(1).mid代表所有字段中的第一个值中mid字段的值
y(1).mid = -d1;
y(1).orig = -d1;
y(1).over = -d1;
% 两端点不记录，拿去端点，还剩(l1-2*delta1)的长度
% 故中间点的中点为l1/2-delta1,原点为l1/2-delta1,终点为l1/2-3/2*delta1
z(1).orig = l1/2 - delta1/2;
z(1).mid = z(1).orig - delta1/2
z(1).over = z(1).orig -delta1;

% 循环结束时，ii+1=n1-1，即正好记录了段数-1个点
for ii=1:n1-2
    y(ii+1) = y(ii);
    z(ii+1).mid = z(ii).mid-delta1;
    z(ii+1).orig=z(ii).orig-delta1;
    z(ii+1).over=z(ii).over-delta1;
end;

%主振子的坐标
%已经有n1-1个点了
%对于主振子，需要记录的点只有n2-1个，其中心点为第n2/2个
y(n1-1+n2/2).mid=0;
y(n1-1+n2/2).orig=0;
y(n1-1+n2/2).over=0;
z(n1-1+n2/2).mid=0;
z(n1-1+n2/2).orig=delta2/2;
z(n1-1+n2/2).over=-delta2/2

for ii=n1-1+n2/2:-1:n1+1
    y(ii-1)=y(ii)
    z(ii-1).mid=z(ii).mid+delta2;
    z(ii-1).orig=z(ii).orig+delta2;
    z(ii-1).over=z(ii).over+delta2;
end
%结束时可得n1+n2-2=(n1-1)+(n2-1)
for ii=n1-1+n2/2:n1+n2-3
    y(ii+1)=y(ii);
    z(ii+1).mid=z(ii).mid-delta2;
    z(ii+1).orig=z(ii).orig-delta2;
    z(ii+1).over=z(ii).over-delta2;
end

%引向器1的坐标
y(n1+n2-1).mid=d2;
y(n1+n2-1).orig=d2;
y(n1+n2-1).over=d2;
z(n1+n2-1).orig=l3/2-delta3/2;
z(n1+n2-1).mid=z(n1+n2-1).orig-delta3/2;
z(n1+n2-1).over=z(n1+n2-1).orig-delta3;

for ii=n1+n2-1:n1+n2+n3-4
    y(ii+1)=y(ii);
    z(ii+1).mid=z(ii).mid-delta3;
    z(ii+1).orig=z(ii).orig-delta3;
    z(ii+1).over=z(ii).over-delta3;
end

%引向器2的坐标
y(n1+n2+n3-2).mid=d3;
y(n1+n2+n3-2).orig=d3;
y(n1+n2+n3-2).over=d3;
z(n1+n2+n3-2).orig=l4/2-delta4/2;
z(n1+n2+n3-2).mid=z(n1+n2+n3-2).orig-delta4/2;
z(n1+n2+n3-2).over=z(n1+n2+n3-2).orig-delta4;

for jj=n1+n2+n3-2:n1+n2+n3+n4-5
    y(jj+1)=y(jj);
    z(jj+1).mid=z(jj).mid-delta3;
    z(jj+1).orig=z(jj).orig-delta3;
    z(jj+1).over=z(jj).over-delta3;
end


% 遍历所有微元之间的组合，计算第m个微元和第n个微元之间的阻抗
% m可以等于或者不等于n
for m=1:N %计算阻抗矩阵元素
    for n=1:N
        % y(m).mid即第m个微元的中点的y坐标
        % z(m).mid即第m个微元的中点的z坐标
        rmn=mydistance(0,y(m).mid,z(m).mid,0,y(n).mid,z(n).mid);
        % m和n微元，末端之间的距离
        r1 = mydistance(0,y(m).over,z(m).over,0,y(n).over,z(n).over);
        % m和n微元，末端和始端的距离
        r2 = mydistance(0,y(m).over,z(m).over,0,y(n).orig,z(n).orig);
        % m和n微元，始端和末端的距离
        r3 = mydistance(0,y(m).orig,z(m).orig,0,y(n).over,z(n).over);
        % m和n微元，始端之间的距离
        r4 = mydistance(0,y(m).orig,z(m).orig,0,y(n).orig,z(n).orig);
        
        % 如果末端和末端重合，那m微元和n微元就是同一个微元了
        % 就套用《简明天线理论与应用设计》的(6-15)中m=n的情况
        % 即令q=q1
        if (r1==0)


            
            q1=1/(2*pi*r2)*log(r2/a)-j*k/(4*pi);
            q2=exp(-j*k*r2)/(4*pi*r2);
            q3=exp(-j*k*r3)/(4*pi*r3);
            q4=q1;
            q=q1;
            zz(m,n)=j*w*4*pi*(10^(-7))*(r2^2)*q+(1/(j*w*8.854*(10^(-12))))*(q1-q2-q3+q4);
        end

        % 末端和始端一致，而且末端不等于末端，始端和末端不一致
        if (r2==0)&(r1~=0)&(r3~=0)
            q1=exp(-j*k*r1)/(4*pi*r1);
            q2=1/(2*pi*(r3/2))*log((r3/2)/a)-j*k/(4*pi);
            q3=exp(-j*k*r3)/(4*pi*r3);
            q4=exp(-j*k*r4)/(4*pi*r4);
            q=exp(-j*k*rmn)/(1*pi*rmn);
            zz(m,n)=j*w*4*pi*(10^(-7))*((r3/2)^2)*q+(1/(j*w*8.854*(10^(-12))))*(q1-q2-q3+q4)
        end
        if (r3==0)&(r1~=0)&(r2~=0)
            q1=exp(-j*k*r1)/(4*pi*r1);
            
            q3=1/(2*pi*(r2/2))*log((r2/2)/a)-j*k/(4*pi);
            q2=exp(-j*k*r2)/(4*pi*r2);
            q4=exp(-j*k*r4)/(4*pi*r4);
            q=exp(-j*k*rmn)/(4*pi*rmn)
            zz(m,n) = j*w*4*pi*(10^(-7))*((r2/2)^2)*q+(1/(j*w*8.854*(10^(-12))))*(q1-q2-q3+q4);
        end
        % 两微元不在同一振子上
        if (r3~=0)&(r1~=0)&(r2~=0)
            q1=exp(-j*k*r1)/(4*pi*r1);
            q2=exp(-j*k*r2)/(4*pi*r2);
            q3=exp(-j*k*r3)/(4*pi*r3);
            q4=exp(-j*k*r4)/(4*pi*r4);
            q=exp(-j*k*rmn)/(4*pi*rmn)
            rm = mydistance(0,y(m).orig,z(m).orig,0,y(m).over,z(m).over);
            rn = mydistance(0,y(n).orig,z(n).orig,0,y(n).over,z(n).over);
            zz(m,n) = j*w*4*pi*(10^(-7))*(rm*rn)*q+(1/(j*w*8.854*(10^(-12))))*(q1-q2-q3+q4);
        end
    end
end

zz;
v = zeros(1,N)';
v(n1-1+n2/2)=1;
distribute = zz\v; %电流分布
zin = 1/(distribute(n1-1+n2/2))
rin = real(zin)
garma = (zin-50)/(zin+50)
vswr = (1+abs(garma))/(1-abs(garma))
dd = 180/step;
f1 = zeros(N,1);
f2 = zeros(N,1);
for step1 = 0:1:2*dd
    for step2 = 0:1:2*dd
        theta = step1*pi/dd;
        fai = step2*pi/dd;
        rx = sin(theta)*cos(fai);
        ry = sin(theta*sin(fai));
        rz = cos(theta);
        r = mydistance(rx,ry,rz,0,0,0);
        for n=1:N
            delta = mydistance(0,y(n).orig,z(n).orig,0,y(n).over,z(n).over)
            neiji = rx*x(n)+ry*y(n).mid+rz*z(n).mid;
            en(n) = (j*w*u*exp(-j*k*r)/(4*pi*r))*exp(j*k*neiji)*delta.*distribute(n)*sin(theta); % 场强
        end
        
        e1(step1+1,step2+1) = sum(en);
        y1(step1+1,step2+1) = abs(e1(step1+1,step2+1));
        g(step1+1,step2+1) = (4*pi*r*r*y1(step1+1,step2+1)^2)/(120*pi*rin*abs(distribute(n1+n2/2-1))^2); %增益
    end
end
gmax = max(g);
gmaxdB = 10*log10(max(gmax))
x = 0:1:2*dd;
y = 0:1:2*dd;

%以下为绘图命令
%[Th,Fi]=meshgrid(x,y);
%X=y1'.*sin(Th*pi/dd).*sin(Fi*pi/dd);
%Y=y1'.*sin(Th*pi/dd).*cos(Fi*pi/dd);
%Z=y1'.*cos(Th*pi/dd);
%figure;
%mesh(X,Y,Z)
%colormap([0 0 0]);
figure;
plot(abs(real(distribute))/abs(max(real(distribute))),'b:o')
xlabel('number');
ylabel('current');
grid on
figure;
polar(y*pi/dd,y1(90/step+1,:)/max(y1(90/step+1,:)),'b-');
hold on
polar(x'*pi/dd,y1(:,90/step+1)/max(y1(:,90/step+1)),'r--');
legend('H面方向图','E面方向图');