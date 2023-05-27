%����l1�˵ĳ�ʼ����
l1 = 0.3;
alpha1 = 0:10:360;
omega1 = 5;
n = length(alpha1);

% ����l2�ĳ�ʼ����
l2 = 0.8;
M3 = 100;
%���ܳ���
m = 0.4;

% Ӧ��������1�ϵ�Mb
Mb = zeros(1,n);
% ֧��A��Լ����
Rax = zeros(1,n);
Ray = zeros(1,n);

% �ƶ���B��Լ����
F = zeros(1,n);

% ֧��C��Լ����
Rcx = zeros(1,n);
Rcy = zeros(1,n);

for iterTime = 1:n
    
%     �˶�ѧ����

syms l22 alpha2   % l22ָ����B��C֮��ľ���


% ʸ������ m+l1=l22
% ��x��
q1 = l1*cosd(alpha1(iterTime))-l22*cosd(alpha2);
q2 = l1*sind(alpha1(iterTime))+m-l22*sind(alpha2);
%���λ�Ʒ��̣��ó� l22 alpha2
%���õ�������ķ�������ǰһʱ�̵�λ����ȷ����һʱ�̵�λ��
if iterTime == 1
   x0 = [53,500];
end
T = vpasolve(q1,q2,x0);
l22 = T.l22;
alpha2 = T.alpha2;
   if iterTime > 1
   x0 = [alpha2,l22];
   end

   
%    ��ʸ��������һ�׵�

syms v omega2
% ��x��
q3 = -l1*sind(alpha1(iterTime))*omega1-v*cosd(alpha2)+l22*sind(alpha2)*omega2;
% ��y��
q4 = l1*cosd(alpha1(iterTime))*omega1-v*sind(alpha2)-l22*cosd(alpha2)*omega2;


P = vpasolve(q3,q4);
v = P.v;
omega2 = P.omega2;

% ��ʸ����������׵�

% �鹦�����Mb
syms Mbb;
eq1 = Mbb*omega1-M3*omega2;
Mb(iterTime) = vpasolve(eq1);

% ����l2����Լ������
syms Rcxx Rcyy FF
% ��C��ȡ��
eq1 = -M3+FF*l22;
% ˮƽ�������Ϊ0
eq2 = Rcxx-FF*sind(alpha2);
% ��ֱ�������Ϊ0
eq3 = Rcyy+FF*cosd(alpha2);


L = vpasolve(eq1,eq2,eq3);
Rcx(iterTime) = L.Rcxx;
Rcxx = L.Rcxx;
Rcy(iterTime) = L.Rcyy;
Rcyy = L.Rcyy;
F(iterTime) = L.FF;
FF = L.FF;


syms Rayy Raxx
% ��A��ȡ��
eq4 = Mbb-FF*l1*cosd(alpha1(iterTime)-alpha2);% �Զ����� �������
% ˮƽ�������Ϊ0
eq5 = Raxx+FF*sind(alpha2);
% ��ֱ�������Ϊ0
eq6 = Rayy-FF*cosd(alpha2);


N = vpasolve(eq5,eq6);
Rax(iterTime) = N.Raxx;
Ray(iterTime) = N.Rayy;

end


figure(1)
plot(alpha1, Rax,'-*');
title('���˻��� ֧��Aˮƽ�����Լ������');
xlabel('alpha1(��)');
ylabel('Rax(N)');


figure(2)
plot(alpha1, Ray,'-o');
title('���˻��� ֧��AǦ�������Լ������');
xlabel('alpha1(��)');
ylabel('Ray(N)');

figure(3)
plot(alpha1, F,'-^');
title('���˻��� �ƶ�����Լ������');
xlabel('alpha1(��)');
ylabel('F(N)');

figure(4)
plot(alpha1, Rcx,'-*');
title('���˻��� ֧��Cˮƽ�����Լ������');
xlabel('alpha1(��)');
ylabel('Rcx(N)');


figure(5)
plot(alpha1, Rcy,'-o');
title('���˻��� ֧��CǦ�������Լ������');
xlabel('alpha1(��)');
ylabel('Rcy(N)');

figure(6)
plot(alpha1, Mb,'-^');
title('���˻��� Ӧ��������1�ϵ�Mb');
xlabel('alpha1(��)');
ylabel('Mb(N*m)');
