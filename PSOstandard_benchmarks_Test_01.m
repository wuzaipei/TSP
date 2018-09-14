function PSOstandard_benchmarks_Test
%   感谢亲亲使用此代码，此代码解决您的问题了吗~(@^_^@)~
%   没解决的话告诉亲亲一个好消息，登录淘宝店铺“大成软件工作室”，可以下载(????)1分钱成品代码(′▽`〃)哦~
%   是的，亲亲真的没有看错，挠破头皮的问题真的1分钱就可以解决了(づ??????)づ
%   小的这就把传送门给您，记得要收藏好哦(づ￣3￣)づ╭?～
%   传送门：https://item.taobao.com/item.htm?spm=a1z10.1-c.w4004-15151018122.5.uwGoq5&id=538759553146
%   如果传送门失效，亲亲可以来店铺讨要，客服MM等亲亲来骚扰哦~(*/ω╲*)
web https://item.taobao.com/item.htm?spm=a1z10.1-c.w4004-15151018122.5.uwGoq5&id=538759553146 -browser
clear all;
close all;
c1=1.49445;c2=1.49445;%
global dimension  Size
dimension=40;Size=40;%种群维数 dimension、规模 Size
Tmax=1000;%%最大迭代次数 Tmax
%%选择不同测试函数的速度和位置限制范围%%
F_n=1;
switch F_n
  case 1   %%  f1_Sphere                    %%
           Vmax(1:dimension)= 30;  Vmin(1:dimension)=-30;
           Xmax(1:dimension)= 30;  Xmin(1:dimension)=-30;
  case 2   %%  f2_Quadric    [-100,100]     %% 
           Vmax(1:dimension)= 100;  Vmin(1:dimension)=-100;
           Xmax(1:dimension)= 100;  Xmin(1:dimension)=-100;
  case 3   %%  f3_Ackley     [-30,30]       %%
           Vmax(1:dimension)= 30;  Vmin(1:dimension)=-30;
           Xmax(1:dimension)= 30;  Xmin(1:dimension)=-30;
  case 4   %%  f4_griewank   [-600,600]     %%
           Vmax(1:dimension)= 600;  Vmin(1:dimension)=-600;
           Xmax(1:dimension)= 600;  Xmin(1:dimension)=-600;
  case 5   %%  f5_Rastrigin  [-5.12,5.12]   %%
           Vmax(1:dimension)= 5.12;  Vmin(1:dimension)=-5.12;
           Xmax(1:dimension)= 5.12;  Xmin(1:dimension)=-5.12;
  case 6   %%  f6_Rosenbrock [-2.408,2.408] %%
           Vmax(1:dimension)= 2.408;  Vmin(1:dimension)=-2.408;
           Xmax(1:dimension)= 2.408;  Xmin(1:dimension)=-2.408;
  case 7   %%  f7_Schaffer's f6 %%
           Vmax(1:dimension)= 2.408;  Vmin(1:dimension)=-2.408;
           Xmax(1:dimension)= 2.408;  Xmin(1:dimension)=-2.408;  
 end
%%三维显示粒子群运动变化%%
global Swarmscope; 
       Swarmscope = plot(0,0, '.');
       axis([Xmin(1) Xmax(1) Xmin(2) Xmax(2) Xmin(3) Xmax(3)]);   %初始轴的范围的设置
    %  axis square;
       grid on;
       set(Swarmscope,'EraseMode','xor','MarkerSize',12); %设置用来显示粒子.
%%initial Position Velocity%%
Position=zeros(dimension,Size);%以后位置Position统一为此种记法：行 dimension；列 Size；
Velocity=zeros(dimension,Size);%每个粒子的位置、速度对应于一列。
[Position,Velocity]=initial_Position_Velocity(dimension,Size,Xmax,Xmin,Vmax,Vmin);
%%个体最优 P_p 和全局最优 globe 初始赋值%%
P_p=Position;globe=zeros(dimension,1);
%%评价每个粒子适应值，寻找出 globle%%
for j=1:Size
    Pos=Position(:,j);
    fz(j)=Fitness_Function(Pos,F_n);
end
[P_g,I]=min(fz);%P_g  1*1 ?
globe=Position(:,I);
%%打散参数设置%%
 N_dismiss=51;%太小，不利于初始寻优
 N_dismissed=0;%记录被打散的次数
 deltaP_gg=0.001%种群过分收敛衡量标准值(适应度变化率)
%  reset = 1;  %设置reset = 1时指示粒子群过分收敛时将被打散，如果reset＝0则不打散
 reset_dismiss = 0;
%%迭代开始%%
for itrtn=1:Tmax
   time(itrtn)=itrtn;
%%过于集中时打散%%
   if reset_dismiss==1
       bit=1;
       if itrtn>N_dismiss
          bit=bit&((P_gg(itrtn-1)-P_gg(itrtn-N_dismiss))/P_gg(itrtn-1)< deltaP_gg);
          if bit==1
             [Position,Velocity]=initial_Position_Velocity(dimension,Size,Xmax,Xmin,Vmax,Vmin);%重新初始化位置和速度
             N_dismissed=N_dismissed+1;
             N_dismissed
             warning('粒子过分集中！重新初始化……');      %   给出信息
             itrtn
          end
       end 
    end
    
     Weight=0.4+0.5*(Tmax-itrtn)/Tmax;
%        Weight=1;
     r1=rand(1);r2=rand(1);
    for i=1:Size
        Velocity(:,i)=Weight*Velocity(:,i)+c1*r1*(P_p(:,i)-Position(:,i))+c2*r2*(globe-Position(:,i));%速度更新
    end
%%速度限制%%
    for i=1:Size
            %%引入速度边界变异%%
%         Vout_max=max(Velocity(:,i));
%         Vout_min=min(Velocity(:,i));
%         if Vout_max
        jj=1;
        K=ones(dimension,1);
        for row=1:dimension
            if Velocity(row,i)>Vmax(row) 
             K(jj)=Vmax(row)/Velocity(row,i);
             jj=jj+1; 
            elseif Velocity(row,i)<Vmin(row)
             K(jj)=Vmin(row)/Velocity(row,i);
             jj=jj+1;
            else
            end
        end
        Kmin=min(K);
        for row=1:dimension
            if Velocity(row,i)>Vmax(row)
                Velocity(row,i)=Velocity(row,i)*Kmin;
            elseif Velocity(row,i)<Vmin(row)
                Velocity(row,i)=Velocity(row,i)*Kmin;
            else
            end
        end
        
     end
%      K
    Position=Position+Velocity;%位置更新
%%位置限制%%
    for i=1:Size
        for row=1:dimension
            if Position(row,i)>Xmax(row)
                Position(row,i)=Xmax(row);
            elseif Position(row,i)<Xmin(row)
                Position(row,i)=Xmin(row);
            else
            end
        end
    end
%%重新评价每个粒子适应值，更新个体最优 P_p 和全局最优 globe%%
    for j=1:Size
        xx=Position(:,j)';
       fz1(j)=Fitness_Function(xx,F_n);
        if fz1(j)<fz(j)
            P_p(:,j)=Position(:,j);
            fz(j)=fz1(j);
        end
     % [P_g1,I]=min(fz1);%%%有改动
        if fz1(j)<P_g
            P_g=fz1(j);
     %      globe=Position(:,I);
        end
    end
     [P_g1 I]=min(fz);
     P_gg(itrtn)=P_g1;
     globe=P_p(:,I);
%     globe=Position(:,I);
%     itrtn
%     globe
%% draw 粒子群运动变化图%%    
       XX=Position(1,:);YY=Position(2,:);ZZ=Position(3,:);
        if dimension>= 3
             set(Swarmscope,'XData',XX,'YData', YY, 'ZData', ZZ);
        elseif dimension== 2
             set(Swarmscope,'XData',XX,'YData',YY );%设置
        end
        xlabel('粒子第一维');
        ylabel('粒子第二维');
        zlabel('粒子第三维');
        drawnow; 
end
%%画‘评价值’变化曲线%%
figure(1);
BestFitness_plot(time,P_gg);

%%画系统阶跃响应变化曲线%%
% figure(2);
% Step_2PID(globe)

function   BestFitness_plot(time,P_gg)
plot(time,P_gg);
xlabel('迭代的次数');ylabel('适应度值P_g');
function [Position,Velocity]=initial_Position_Velocity(dimension,Size,Xmax,Xmin,Vmax,Vmin)
   for i=1:dimension
      Position(i,:)=Xmin(i)+(Xmax(i)-Xmin(i))*rand(1,Size);
      Velocity(i,:)=Vmin(i)+(Vmax(i)-Vmin(i))*rand(1,Size);
   end

function Fitness=Fitness_Function(X,F_n)
  global dimension  Size
% F_n 标准测试函数选择,其中:
% n=1: f1_Sphere        测试   
% n=2: f2_Quadric       测试
% n=3: f3_Ackley        测试  
% n=4: f4_Griewank      测试
% n=5: f5_Rastrigin     测试
% n=6: f6_Rosenbrock    测试
% n=7: f7_Schaffer's f6 测试  注：此函数只接受两个变量，故dimension＝2。
switch F_n
  case 1
    %%  f1_Sphere        %%
       Func_Rastrigin=X(:)'*X(:);
       Fitness=Func_Rastrigin;
  case 2
    %%  f2_Quadric       %%
       res1=0;
       for row=1:dimension
           res1=res1+(sum(X(1:row)))^2;
       end
       Func_Quadric=res1;
       Fitness=Func_Quadric;
  case 3
       %%  f3_Ackley        %%
       Func_Ackley=-20*exp(-0.2*sqrt((1/dimension)*(X(:)'*X(:))))-exp((1/dimension)*((cos(2*pi*X(:)')*cos(2*pi*X(:)))))+20+exp(1);
       Fitness=Func_Ackley;
  case 4
      %%  f4_griewank      %%
      res1=X(:)'*X(:)/4000;
      res2=1;
      for row=1:dimension
          res2=res2*cos(X(row)/sqrt(row));
      end
      Func_Griewank=res1-res2+1;
      Fitness=Func_Griewank;
  case 5 
      %%  f5_Rastrigin     %%
      Func_Rastrigin=X(:)'*X(:)-10*sum(cos(X(:)*2*pi))+10*dimension;
      Fitness=Func_Rastrigin;
  case 6
      %%  f6_Rosenbrock    %%
      res1=0;
      for row=1:(dimension-1)
          res1=res1+100*(X(row+1)-X(row)^2)^2+(X(row)-1)^2;
      end
      Func_Rosenbrock=res1;
      Fitness=Func_Rosenbrock;
   case 7
      %%  f7_Schaffer's f6 %%
      Func_Schaffer=0.5-(sin(sqrt(X(1)^2+X(2)^2))^2-0.5)/(1+0.001*(X(1)^2+X(2)^2))^2;
      Fitness=Func_Schaffer;            
end