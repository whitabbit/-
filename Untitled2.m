%��Ⱥ�㷨

clear all;%��չ����ռ�
close all;%�ر�����Figure���� 
clc;%��������

t1=clock;%��¼��������ʱ��

% ��ʼ����������
X1=[10 3 1 5 4 0 3 7 9 10 14 2 6 11 15 19 22 21 27 15];
Y1=[10 2 5 4 7 8 11 9 6 2 0 16 18 17 12 9 5 0 9 19]; 
position=[X1',Y1'];%������������ɵ�
number = length(position);%ʳƷ���۵������
dis=zeros(number);

%�����ڽӾ���
for n=1:number
    for m=1:number
        dis(n,m)=abs(position(n,1)-position(m,1))+abs(position(n,2)-position(m,2));
    end
end

%��Ⱥ�㷨��������ʼ��
alpha =1;%��������,��Ϣ�ص���Ҫ�̶�
beta = 3;%��������,ʳƷ���۵��������Ҫ�̶�
rou = 0.7;%��Ϣ�ز�������
Q = 2000;%��Ϣ����ǿϵ��
m = 40;%��������
Eta = 1./dis;%��ʾÿ���ߵ��ܼ���
Tau = ones(number,number);%Ϊ��Ϣ�ؾ���,�����ϻ�û�б�����ʳƷ���۵�ǰ��ÿ����·�Ͼ��Ѿ�������һ����������Ϣ��
Tauroute = zeros(m,number);%�洢����¼ÿ�ε���ʱÿֻ���Ͼ�����·������
NC_max = 50;%��������
R_best=zeros(NC_max,number);%�������·��
L_best=inf.*ones(NC_max,1);%�������·�ߵĳ���
NC = 1;%%��¼��ǰ���ϴ���

%��¼ÿһ������·��
X2 = [];
Y2 = [];

while NC <=NC_max%��ʾѭ����ֹ������������ֹ��
    %��mֻ���Ϸŵ�n��ʳƷ���۵���
    Randpos = [];%�漴��ȡ
    for i = 1:(ceil(m/number))%��mֻ��������ŵ�n��ʳƷ���۵���
        Randpos = [Randpos,randperm(number)];
    end
    Tauroute(:,1) = (Randpos(1,1:m))';%��ʼ��ÿֻ���ϵ�·��,1��Ϊÿֻ���ϵ���ʼ��
    %�������ϰ����ʺ���ѡ����һ��ʳƷ���۵�,��ɸ��Ե�����
    for j = 2:number
        for i = 1:m
            visited = Tauroute(i,1:(j-1));%�Ѿ��ҹ���ʳƷ���۵�
            prob = zeros(1,(number-j+1));%������ʳƷ���۵��ѡ����ʷֲ�,
            positions=1:number;
            visting=positions(ismember(positions,visited)==0);%�����ʵ�ʳƷ���۵�
            %visting=setdiff([1:number],visited);
            %��������ѡʳƷ���۵�ĸ��ʷֲ�
            for k = 1:length(visting)%��ÿֻ���ϻ�û�з��ʵ�ʳƷ���۵����μ������
                %Eta��ʾÿ���ߵ��ܼ���
                %TauΪ��Ϣ�ؾ���,ÿ�����Ѿ����ڵ���Ϣ��
                prob(k) = (Tau(visited(end),visting(k))^alpha)*(Eta(visited(end),visting(k))^beta);
            end
            prob = prob/(sum(prob));%ת��Ϊ����
            pcum = cumsum(prob);%��������
            Select = find(pcum>=rand);%���������ѡȡ��һ��Ҫѡ��ʳƷ���۵�
            tovisit = visting(Select(1));%���Ѿ�ѡ���ʳƷ���۵�������ѡһ��ʳƷ���۵�
            Tauroute(i,j) = tovisit;%��·���м�¼ѡ���ʳƷ���۵�
        end
    end
    if NC > 2
        Tauroute(1,:) = R_best(NC-1,:);%���������������2������һ�����·�����ǵ�һֻ���ϵ�·��
    end
    %��¼���ε������·��
    L = zeros(m,1);%��ʼ����Ϊ0,m*1��������
    for i = 1:m%����ÿһ�����ϵ�·��
        R = Tauroute(i,:);
        d = dis(R(end),R(1));
        for k = 1:length(R)-1%����·��
            d = d + dis(R(k),R(k+1));
        end
        L(i)=d;
    end
    L_best(NC) = min(L);%��¼���ε�����̵�·���ĳ���(L_best��һ��������)
	pos = find(L==L_best(NC));%�ҵ�����������Ӧ���±�
	R_best(NC,:) = Tauroute(pos(1),:);%ͨ�����·�����±��ҵ�����·������¼
	Rbest = Tauroute(pos(1),:);%��ͼ
	x=position(Rbest,1);
	y=position(Rbest,2);
    subplot(1,2,1);%���Ƶ�ǰ�����µ�·��ͼ
	plot(X1,Y1,'.k','LineWidth',5);%��עʳƷ���۵�
    if NC == 1
        figure(1)
        plot(X1(1),Y1(1),'p','markersize',10, ...
            'MarkerEdgeColor','r','MarkerFaceColor','g');%������ʼ��
        plot([x;x(1)],[y;y(1)],'b');%����ȱ��
        axis([-5,30,-2,20]);%�޶�����ϵ��Χ
        %��ÿ������б�ע
        for point_num=1:number
            text(X1(point_num),Y1(point_num),strcat(num2str(mod((point_num+18),20)+1),'(',num2str(X1(point_num)),',',num2str(Y1(point_num)),')'));
        end
        L_best(NC);
        xlabel(sprintf('��Ⱥ���� = %5i      ��ǰ��Ⱥ�����µ��ܾ��� = %6.1f',NC,L_best(NC)));%������ÿ����̾���
        drawnow;
    end
    figure(2);
    subplot(1,2,1);
	hold on;
    axis([-5,30,-2,20]);%�޶�����ϵ��Χ
	plot(X1(1),Y1(1),'p','markersize',10, ...
        'MarkerEdgeColor','r','MarkerFaceColor','g');%������ʼ��
	plot([x;x(1)],[y;y(1)],'b');%����ȱ��
    %��ÿ������б�ע
    for point_num=1:number
        text(X1(point_num),Y1(point_num),strcat(num2str(mod((point_num+18),20)+1),'(',num2str(X1(point_num)),',',num2str(Y1(point_num)),')'));
    end
    hold off;
	L_best(NC);
	xlabel(sprintf('��Ⱥ���� = %5i      ��ǰ��Ⱥ�����µ��ܾ��� = %6.1f',NC,L_best(NC)));%������ÿ����̾���
	drawnow;
    X2(end+1) = NC;
    Y2(end+1) = L_best(NC);
	NC = NC + 1;%����������
	%������Ϣ��
	Delta_Tau = zeros(number,number);%Ĭ����Ϣ�ؾ���Ϊȫ��
	for i = 1:m%��ÿֻ���ϵ�·��������Ϣ��
        for j = 1:(number-1)
            Delta_Tau(Tauroute(i,j),Tauroute(i,j+1)) = Delta_Tau(Tauroute(i,j),Tauroute(i,j+1))+Q/L(i);
        end
        %���һ��ʳƷ���۵�ص���ʼʳƷ���۵�
        Delta_Tau(Tauroute(i,number),Tauroute(i,1))=Delta_Tau(Tauroute(i,number),Tauroute(i,1))+Q/L(i);
    end
	Tau = (1-rou).*Tau + Delta_Tau;%�����������¼ӵ���Ϣ��
	Tauroute = zeros(m,number);%Tauroute����������һ��·���ļ�¼
end
figure(2);
subplot(1,2,2);
hold on;
axis([1,55,100,200]);
plot(X2,Y2,'b.-');%�������ϴ�������·�̵�ͼ
text(X2(1),Y2(1),strcat('(',num2str(X2(1)),',',num2str(Y2(1)),')'));
text(X2(50),Y2(50),strcat('(',num2str(X2(50)),',',num2str(Y2(50)),')'));
t2=clock;
etime(t2,t1)%��������ʱ��