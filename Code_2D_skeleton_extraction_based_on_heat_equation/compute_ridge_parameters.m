function [P,K,direct_max,direct_min ,plane_normal,other_dir] = compute_ridge_parameters( fx,fy,fxx,fxy,fyy,iidx )
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
% P ��ʾ�������������
% K ��ʾ���������ֵ
% direct_max direct_min�ֱ��ʾ������С������ ���ж��ǰ��д洢
u1=ones(size(fx));
u2=zeros(size(fx));
u3=fx;
u=[u1,u2,u3];

v1=zeros(size(fy));
v2=ones(size(fy));
v3=fy;
v=[v1,v2,v3];

A11=ones(size(fx))+fx.*fx;
A12=fx.*fy;
A21=fx.*fy;
A22=ones(size(fy))+fy.*fy;

B0=(ones(size(fx))+fx.*fx+fy.*fy).^(-1/2);
B11=B0.*fxx;
B12=B0.*fxy;
B21=B0.*fxy;
B22=B0.*fyy;
P=zeros(2,2,length(fx)); %�洢��������
K=zeros(length(fx),2);%�洢����ֵ
direct_max=zeros(length(fx),3);
direct_min=zeros(length(fx),3);
plane_normal=zeros(length(fx),3);
other_dir=zeros(length(fx),3);
disp('wait for eig');
for  i=1:length(fx)
   A=[A11(i),A12(i);A21(i),A22(i)];
   B=[B11(i),B12(i);B21(i),B22(i)];
   C=A\B;
   [Pi,Ki]=eig(C); %�������-����������� ���ﲻ�ӷ� ��V(:,3)ȡ��

   if (Ki(1,1)<Ki(2,2))
       tmp=Ki(1,1);
       Ki(1,1)=Ki(2,2);
       Ki(2,2)=tmp;
       tmp_v=Pi(:,1);
       Pi(:,1)=Pi(:,2);
       Pi(:,2)=tmp_v;
   end %�ô������ֵ��ǰ
    P(:,:,i)=Pi;
    K(i,1)=Ki(1,1);
    K(i,2)=Ki(2,2); % K��һ�б�ʾ�ϴ�����ֵ �ڶ����ǽ�С����ֵ
    direct_max(i,:)=Pi(1,1)*u(i,:) + Pi(2,1)*v(i,:);
    direct_min(i,:)=Pi(1,2)*u(i,:) + Pi(2,2)*v(i,:);
    plane_normal(i,:)=cross(u(i,:),v(i,:));
    other_dir(i,:)=cross(plane_normal(i,:),direct_max(i,:));
end

end

