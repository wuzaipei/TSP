function PSOstandard_benchmarks_Test
%   ��л����ʹ�ô˴��룬�˴�����������������~(@^_^@)~
%   û����Ļ���������һ������Ϣ����¼�Ա����̡������������ҡ�����������(????)1��Ǯ��Ʒ����(�䨌`��)Ŷ~
%   �ǵģ��������û�п�������ͷƤ���������1��Ǯ�Ϳ��Խ����(��??????)��
%   С����ͰѴ����Ÿ������ǵ�Ҫ�ղغ�Ŷ(�ţ�3��)�Ũq?��
%   �����ţ�https://item.taobao.com/item.htm?spm=a1z10.1-c.w4004-15151018122.5.uwGoq5&id=538759553146
%   ���������ʧЧ�����׿�����������Ҫ���ͷ�MM��������ɧ��Ŷ~(*/�بv*)
web https://item.taobao.com/item.htm?spm=a1z10.1-c.w4004-15151018122.5.uwGoq5&id=538759553146 -browser
clear all;
close all;
c1=1.49445;c2=1.49445;%
global dimension  Size
dimension=40;Size=40;%��Ⱥά�� dimension����ģ Size
Tmax=1000;%%���������� Tmax
%%ѡ��ͬ���Ժ������ٶȺ�λ�����Ʒ�Χ%%
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
%%��ά��ʾ����Ⱥ�˶��仯%%
global Swarmscope; 
       Swarmscope = plot(0,0, '.');
       axis([Xmin(1) Xmax(1) Xmin(2) Xmax(2) Xmin(3) Xmax(3)]);   %��ʼ��ķ�Χ������
    %  axis square;
       grid on;
       set(Swarmscope,'EraseMode','xor','MarkerSize',12); %����������ʾ����.
%%initial Position Velocity%%
Position=zeros(dimension,Size);%�Ժ�λ��PositionͳһΪ���ּǷ����� dimension���� Size��
Velocity=zeros(dimension,Size);%ÿ�����ӵ�λ�á��ٶȶ�Ӧ��һ�С�
[Position,Velocity]=initial_Position_Velocity(dimension,Size,Xmax,Xmin,Vmax,Vmin);
%%�������� P_p ��ȫ������ globe ��ʼ��ֵ%%
P_p=Position;globe=zeros(dimension,1);
%%����ÿ��������Ӧֵ��Ѱ�ҳ� globle%%
for j=1:Size
    Pos=Position(:,j);
    fz(j)=Fitness_Function(Pos,F_n);
end
[P_g,I]=min(fz);%P_g  1*1 ?
globe=Position(:,I);
%%��ɢ��������%%
 N_dismiss=51;%̫С�������ڳ�ʼѰ��
 N_dismissed=0;%��¼����ɢ�Ĵ���
 deltaP_gg=0.001%��Ⱥ��������������׼ֵ(��Ӧ�ȱ仯��)
%  reset = 1;  %����reset = 1ʱָʾ����Ⱥ��������ʱ������ɢ�����reset��0�򲻴�ɢ
 reset_dismiss = 0;
%%������ʼ%%
for itrtn=1:Tmax
   time(itrtn)=itrtn;
%%���ڼ���ʱ��ɢ%%
   if reset_dismiss==1
       bit=1;
       if itrtn>N_dismiss
          bit=bit&((P_gg(itrtn-1)-P_gg(itrtn-N_dismiss))/P_gg(itrtn-1)< deltaP_gg);
          if bit==1
             [Position,Velocity]=initial_Position_Velocity(dimension,Size,Xmax,Xmin,Vmax,Vmin);%���³�ʼ��λ�ú��ٶ�
             N_dismissed=N_dismissed+1;
             N_dismissed
             warning('���ӹ��ּ��У����³�ʼ������');      %   ������Ϣ
             itrtn
          end
       end 
    end
    
     Weight=0.4+0.5*(Tmax-itrtn)/Tmax;
%        Weight=1;
     r1=rand(1);r2=rand(1);
    for i=1:Size
        Velocity(:,i)=Weight*Velocity(:,i)+c1*r1*(P_p(:,i)-Position(:,i))+c2*r2*(globe-Position(:,i));%�ٶȸ���
    end
%%�ٶ�����%%
    for i=1:Size
            %%�����ٶȱ߽����%%
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
    Position=Position+Velocity;%λ�ø���
%%λ������%%
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
%%��������ÿ��������Ӧֵ�����¸������� P_p ��ȫ������ globe%%
    for j=1:Size
        xx=Position(:,j)';
       fz1(j)=Fitness_Function(xx,F_n);
        if fz1(j)<fz(j)
            P_p(:,j)=Position(:,j);
            fz(j)=fz1(j);
        end
     % [P_g1,I]=min(fz1);%%%�иĶ�
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
%% draw ����Ⱥ�˶��仯ͼ%%    
       XX=Position(1,:);YY=Position(2,:);ZZ=Position(3,:);
        if dimension>= 3
             set(Swarmscope,'XData',XX,'YData', YY, 'ZData', ZZ);
        elseif dimension== 2
             set(Swarmscope,'XData',XX,'YData',YY );%����
        end
        xlabel('���ӵ�һά');
        ylabel('���ӵڶ�ά');
        zlabel('���ӵ���ά');
        drawnow; 
end
%%��������ֵ���仯����%%
figure(1);
BestFitness_plot(time,P_gg);

%%��ϵͳ��Ծ��Ӧ�仯����%%
% figure(2);
% Step_2PID(globe)

function   BestFitness_plot(time,P_gg)
plot(time,P_gg);
xlabel('�����Ĵ���');ylabel('��Ӧ��ֵP_g');
function [Position,Velocity]=initial_Position_Velocity(dimension,Size,Xmax,Xmin,Vmax,Vmin)
   for i=1:dimension
      Position(i,:)=Xmin(i)+(Xmax(i)-Xmin(i))*rand(1,Size);
      Velocity(i,:)=Vmin(i)+(Vmax(i)-Vmin(i))*rand(1,Size);
   end

function Fitness=Fitness_Function(X,F_n)
  global dimension  Size
% F_n ��׼���Ժ���ѡ��,����:
% n=1: f1_Sphere        ����   
% n=2: f2_Quadric       ����
% n=3: f3_Ackley        ����  
% n=4: f4_Griewank      ����
% n=5: f5_Rastrigin     ����
% n=6: f6_Rosenbrock    ����
% n=7: f7_Schaffer's f6 ����  ע���˺���ֻ����������������dimension��2��
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