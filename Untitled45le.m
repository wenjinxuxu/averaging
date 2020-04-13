%equation: D^alpha X(t)=epsilon*X(t)*2*(sint)^2+sqrt(epsilon)*dB(t)/dt
alpha=0.85;%分数阶数
ga=gamma(alpha);
X0=0.1; %初值

epsilon=0.001;

bata=0.3;%levy
r=0.6;%levy测度常数
gd=0.5;%轨道
r1=1+r*gd^(4-bata)/( sqrt(epsilon)*(4-bata) );%平均后常数
jf=( gd^(4-bata)/(4-bata) );%积分


Dt=0.1; %时间步长
M=1000;  %离散点个数
t=0:Dt:Dt*M; %总时间步长

dB=sqrt(Dt)*randn(1,M+1);
B=[0 cumsum(dB)];

%dl=Dt^(1/alpha)*stblrnd(alpha,0,1,0,1,M+1);
%I=dl<0;
%dl(I)=[];
%dl=dl(1:M+1);

x=zeros(1,M+1); %X(t)向量
x(1)=X0;

for i=2:M+1
    k=i-1; %第k个时刻
    a=0;
    c=0;
    %校正系数
    for j=0:k-1;
        if j==0;
            %a=(Dt^alpha)/[alpha*(alpha+1)]*[(k-1)^(alpha+1)-(k-alpha-1)*(k^alpha)]*[ epsilon*x(1)*(cos(0))^2 + epsilon*2*(sin(0))^2 ];%(-x(1));%*f(0,X0);
            a=(Dt^alpha)/[alpha*(alpha+1)]*[(k-1)^(alpha+1)-(k-alpha-1)*(k^alpha)]*[ ( epsilon*2*(cos(0))^2+sqrt(epsilon)*r*2*(sin(0))^2*jf )*x(1)+sqrt(epsilon)*dB(1)/Dt ];
                       %(-x(1));%*f(0,X0);
           
        else
            %a=a0+(Dt^alpha)/[alpha*(alpha+1)]*[(k-j+1)^(alpha+1)-2*(k-j)^(alpha+1)+(k-j-1)^(alpha+1)]*[ epsilon*x(j+1)*(cos(Dt*j))^2 + epsilon*2*(sin(Dt*j))^2 ];%(-x(j+1));%f(Dt*j,x(j+1));
            a=a0+(Dt^alpha)/[alpha*(alpha+1)]*[(k-j+1)^(alpha+1)-2*(k-j)^(alpha+1)+(k-j-1)^(alpha+1)]*[ ( epsilon*2*(cos(Dt*j))^2+sqrt(epsilon)*r*2*(sin(Dt*j))^2*jf )*x(j+1)+sqrt(epsilon)*dB(j+1)/Dt];
                        %(-x(j+1));%f(Dt*j,x(j+1));
            
        end
        a0=a;
    end;
    
    %预估系数
    for j=0:k-1;   
        %c=c+(Dt^alpha)/alpha*[ (k-j)^alpha-(k-j-1)^alpha ]*[ epsilon*x(j+1)*(cos(Dt*j))^2 + epsilon*2*(sin(Dt*j))^2 ];%(-x(j+1));%*f(Dt*j,x(j+1));
        c=c+(Dt^alpha)/alpha*[ (k-j)^alpha-(k-j-1)^alpha ]*[ ( epsilon*2*(cos(Dt*j))^2+sqrt(epsilon)*r*2*(sin(Dt*j))^2*jf )*x(j+1)+sqrt(epsilon)*dB(j+1)/Dt];
                        %+sqrt(epsilon)*r*x(j+1)*2*(sin(Dt*j))^2*( gd^(4-bata)/(4-bata) )];%(-x(j+1));%*f(Dt*j,x(j+1));%(-x(j+1));%*f(Dt*j,x(j+1));
        
    end
    up=X0+1/ga*c;

    %x(i)=X0+1/ga* [ a0+(Dt^alpha)/[alpha*(alpha+1)]*[epsilon*up*(cos(Dt*k))^2 + epsilon*2*(sin(Dt*k))^2] ] ;%( -up );%f(Dt*k,u); 校正公式
    x(i)=X0+1/ga* [ a0+(Dt^alpha)/[alpha*(alpha+1)]*[( epsilon*2*(cos(Dt*k)^2+sqrt(epsilon)*r*2*(sin(0))^2*jf) )*up+sqrt(epsilon)*dB(i)/Dt] ] ;
        %( -up );%f(Dt*k,u); 校正公式
    
    
end
plot(t,x,'m.');
set(gca,'FontSize',13);%grid on;

%axis([0 Dt*M -1 14]);

xlabel ('$t$','Interpreter','latex')
%ylabel ('$X_t$','Interpreter','latex')   Times New Roman
set(gca,'FontSize',13);
hold on;

z=zeros(1,M+1);
z(1)=X0;
for i=2:M+1
    k=i-1; %第k个时刻
    a=0;
    c=0;
    %校正系数
    for j=0:k-1;
        if j==0;
            %a=(Dt^alpha)/[alpha*(alpha+1)]*[(k-1)^(alpha+1)-(k-alpha-1)*(k^alpha)]*[ epsilon*z(1)*0.5 + epsilon*1 ];%(-x(1));%*f(0,X0);
            a=(Dt^alpha)/[alpha*(alpha+1)]*[(k-1)^(alpha+1)-(k-alpha-1)*(k^alpha)]*[ epsilon*z(1)*r1+sqrt(epsilon)*dB(1)/Dt ];%(-x(1));%*f(0,X0); 
            
        else
            %a=a0+(Dt^alpha)/[alpha*(alpha+1)]*[(k-j+1)^(alpha+1)-2*(k-j)^(alpha+1)+(k-j-1)^(alpha+1)]*[ epsilon*z(j+1)*0.5 + epsilon*1 ];%(-x(j+1));%f(Dt*j,x(j+1));
            a=a0+(Dt^alpha)/[alpha*(alpha+1)]*[(k-j+1)^(alpha+1)-2*(k-j)^(alpha+1)+(k-j-1)^(alpha+1)]*[ epsilon*z(j+1)*r1+sqrt(epsilon)*dB(j+1)/Dt ];%(-x(j+1));%f(Dt*j,x(j+1));
           
        end
        a0=a;
    end;
    
    %预估系数
    for j=0:k-1;   
        %c=c+(Dt^alpha)/alpha*[ (k-j)^alpha-(k-j-1)^alpha ]*[ epsilon*z(j+1)*0.5 + epsilon*1 ];%(-x(j+1));%*f(Dt*j,x(j+1));
        c=c+(Dt^alpha)/alpha*[ (k-j)^alpha-(k-j-1)^alpha ]*[ epsilon*z(j+1)*r1+sqrt(epsilon)*dB(j+1)/Dt ];%(-x(j+1));%*f(Dt*j,x(j+1));
        
    end
    up=X0+1/ga*c;

    %z(i)=X0+1/ga* [ a0+(Dt^alpha)/[alpha*(alpha+1)]*[epsilon*up*0.5 + epsilon*1] ];%( -up );%f(Dt*k,u); 校正公式
    z(i)=X0+1/ga* [ a0+(Dt^alpha)/[alpha*(alpha+1)]*[epsilon*up*r1+sqrt(epsilon)*dB(i)/Dt] ];%( -up );%f(Dt*k,u); 校正公式
    
    
end
plot(t,z,'b-');

error=zeros(1,M+1);
for i=1:M+1
    error(i)=( (x(i)-z(i))^2 )^(1/2);%x(i)-z(i);
end
plot(t,error,'g-');
legend('X_\epsilon(t)','Z_\epsilon(t)','Er')
