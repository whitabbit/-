%蚁群算法

clear all;%清空工作空间
close all;%关闭所有Figure窗口 
clc;%清空命令窗口

t1=clock;%记录程序运行时间

% 初始化各点坐标
X1=[10 3 1 5 4 0 3 7 9 10 14 2 6 11 15 19 22 21 27 15];
Y1=[10 2 5 4 7 8 11 9 6 2 0 16 18 17 12 9 5 0 9 19]; 
position=[X1',Y1'];%坐标是随机生成的
number = length(position);%食品销售点的数量
dis=zeros(number);

%生成邻接矩阵
for n=1:number
    for m=1:number
        dis(n,m)=abs(position(n,1)-position(m,1))+abs(position(n,2)-position(m,2));
    end
end

%蚁群算法各参数初始化
alpha =1;%启发因子,信息素的重要程度
beta = 3;%期望因子,食品销售点间距离的重要程度
rou = 0.7;%信息素残留参数
Q = 2000;%信息素增强系数
m = 40;%蚂蚁数量
Eta = 1./dis;%表示每条边的能见度
Tau = ones(number,number);%为信息素矩阵,在蚂蚁还没有被放入食品销售点前，每条道路上就已经存在了一定含量的信息素
Tauroute = zeros(m,number);%存储并记录每次迭代时每只蚂蚁经历的路径生成
NC_max = 50;%迭代次数
R_best=zeros(NC_max,number);%各代最佳路线
L_best=inf.*ones(NC_max,1);%各代最佳路线的长度
NC = 1;%%记录当前蚂蚁代数

%记录每一代的总路程
X2 = [];
Y2 = [];

while NC <=NC_max%表示循环终止条件，迭代终止器
    %将m只蚂蚁放到n个食品销售点上
    Randpos = [];%随即存取
    for i = 1:(ceil(m/number))%将m只蚂蚁随机放到n个食品销售点上
        Randpos = [Randpos,randperm(number)];
    end
    Tauroute(:,1) = (Randpos(1,1:m))';%初始化每只蚂蚁的路径,1列为每只蚂蚁的起始点
    %所有蚂蚁按概率函数选择下一个食品销售点,完成各自的周游
    for j = 2:number
        for i = 1:m
            visited = Tauroute(i,1:(j-1));%已经找过的食品销售点
            prob = zeros(1,(number-j+1));%待访问食品销售点的选择概率分布,
            positions=1:number;
            visting=positions(ismember(positions,visited)==0);%待访问的食品销售点
            %visting=setdiff([1:number],visited);
            %下面计算待选食品销售点的概率分布
            for k = 1:length(visting)%对每只蚂蚁还没有访问的食品销售点依次计算概率
                %Eta表示每条边的能见度
                %Tau为信息素矩阵,每条边已经存在的信息素
                prob(k) = (Tau(visited(end),visting(k))^alpha)*(Eta(visited(end),visting(k))^beta);
            end
            prob = prob/(sum(prob));%转换为概率
            pcum = cumsum(prob);%生成轮盘
            Select = find(pcum>=rand);%按概率随机选取下一个要选的食品销售点
            tovisit = visting(Select(1));%从已经选择的食品销售点们中再选一个食品销售点
            Tauroute(i,j) = tovisit;%在路径中记录选择的食品销售点
        end
    end
    if NC > 2
        Tauroute(1,:) = R_best(NC-1,:);%如果迭代次数大于2，则将上一代最好路径覆盖第一只蚂蚁的路径
    end
    %记录本次迭代最佳路线
    L = zeros(m,1);%开始距离为0,m*1的列向量
    for i = 1:m%计算每一种蚂蚁的路径
        R = Tauroute(i,:);
        d = dis(R(end),R(1));
        for k = 1:length(R)-1%计算路径
            d = d + dis(R(k),R(k+1));
        end
        L(i)=d;
    end
    L_best(NC) = min(L);%记录本次迭代最短的路径的长度(L_best是一个列向量)
	pos = find(L==L_best(NC));%找到该列向量对应的下标
	R_best(NC,:) = Tauroute(pos(1),:);%通过最短路径的下标找到它的路径并记录
	Rbest = Tauroute(pos(1),:);%作图
	x=position(Rbest,1);
	y=position(Rbest,2);
    subplot(1,2,1);%绘制当前代数下的路径图
	plot(X1,Y1,'.k','LineWidth',5);%标注食品销售点
    if NC == 1
        figure(1)
        plot(X1(1),Y1(1),'p','markersize',10, ...
            'MarkerEdgeColor','r','MarkerFaceColor','g');%标明起始点
        plot([x;x(1)],[y;y(1)],'b');%补上缺口
        axis([-5,30,-2,20]);%限定坐标系范围
        %对每个点进行标注
        for point_num=1:number
            text(X1(point_num),Y1(point_num),strcat(num2str(mod((point_num+18),20)+1),'(',num2str(X1(point_num)),',',num2str(Y1(point_num)),')'));
        end
        L_best(NC);
        xlabel(sprintf('蚁群代数 = %5i      当前蚁群代数下的总距离 = %6.1f',NC,L_best(NC)));%代数和每代最短距离
        drawnow;
    end
    figure(2);
    subplot(1,2,1);
	hold on;
    axis([-5,30,-2,20]);%限定坐标系范围
	plot(X1(1),Y1(1),'p','markersize',10, ...
        'MarkerEdgeColor','r','MarkerFaceColor','g');%标明起始点
	plot([x;x(1)],[y;y(1)],'b');%补上缺口
    %对每个点进行标注
    for point_num=1:number
        text(X1(point_num),Y1(point_num),strcat(num2str(mod((point_num+18),20)+1),'(',num2str(X1(point_num)),',',num2str(Y1(point_num)),')'));
    end
    hold off;
	L_best(NC);
	xlabel(sprintf('蚁群代数 = %5i      当前蚁群代数下的总距离 = %6.1f',NC,L_best(NC)));%代数和每代最短距离
	drawnow;
    X2(end+1) = NC;
    Y2(end+1) = L_best(NC);
	NC = NC + 1;%代数自增长
	%更新信息素
	Delta_Tau = zeros(number,number);%默认信息素矩阵为全零
	for i = 1:m%对每只蚂蚁的路径更新信息素
        for j = 1:(number-1)
            Delta_Tau(Tauroute(i,j),Tauroute(i,j+1)) = Delta_Tau(Tauroute(i,j),Tauroute(i,j+1))+Q/L(i);
        end
        %最后一个食品销售点回到起始食品销售点
        Delta_Tau(Tauroute(i,number),Tauroute(i,1))=Delta_Tau(Tauroute(i,number),Tauroute(i,1))+Q/L(i);
    end
	Tau = (1-rou).*Tau + Delta_Tau;%蒸发量加上新加的信息素
	Tauroute = zeros(m,number);%Tauroute清零用于下一轮路径的记录
end
figure(2);
subplot(1,2,2);
hold on;
axis([1,55,100,200]);
plot(X2,Y2,'b.-');%绘制蚂蚁代数与总路程的图
text(X2(1),Y2(1),strcat('(',num2str(X2(1)),',',num2str(Y2(1)),')'));
text(X2(50),Y2(50),strcat('(',num2str(X2(50)),',',num2str(Y2(50)),')'));
t2=clock;
etime(t2,t1)%返回运行时间