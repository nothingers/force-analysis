% ����������������ĳ�ʼ����
l1 = 0.1;
l2 = 0.33;
omega1 = 50*pi;
m2 = 25/9.8;
m3 = 21/9.8;
J2 = 0.0425;
alpha1 = 45;
acceleration1 = 0;
% ʸ������ l1=S+l2
syms S alpha2
% ��x��
q1 = l1*cosd(alpha1)-S-l2*cosd(alpha2);
% ��y��
q2 = l1*sind(alpha1)-l2*sind(alpha2); 

x0=[0.393,168];
T = vpasolve(q1,q2,x0);
S = T.S;
alpha2 = T.alpha2;

% ʸ�����̶�ʱ����һ�׵�
syms v omega2
% ��x��
q3 = -l1*sind(alpha1)*omega1-v+l2*sind(alpha2)*omega2;
% ��y��
q4 = l1*cosd(alpha1)*omega1-l2*cosd(alpha2)*omega2;

P = vpasolve(q3,q4);
v = P.v;
omega2 = P.omega2;

% ʸ�����̶�ʱ������׵�
syms a acceleration2
% ��x��
q5 = -l1*cosd(alpha1)*omega1^2-l1*sind(alpha1)*acceleration1-a+l2*cosd(alpha2)*omega2^2+l2*sind(alpha2)*acceleration2;
% ��y��
q6 = -l1*sind(alpha1)*omega1^2+l1*cosd(alpha1)*acceleration1+l2*sind(alpha2)*omega2^2-l2*cosd(alpha2)*acceleration2;

Q = vpasolve(q5,q6);
a = Q.a;
acceleration2 = Q.acceleration2;

% l2�˵Ĺ��������
FIcx = double(+m2*(omega1^2*l1*cosd(alpha1)-omega2^2*(1/3*l2)*cosd(alpha2)-(1/3*l2)*acceleration2*sind(alpha2)));
FIcy = double(+m2*(omega1^2*l1*sind(alpha1)+omega2^2*(1/3*l2)*sind(alpha2)+(1/3*l2)*acceleration2*cosd(alpha2)));
MI = double(-J2*acceleration2);
% ����Ĺ���������
FI = double(-m3*a);
