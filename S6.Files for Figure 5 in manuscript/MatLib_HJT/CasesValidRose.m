function CasesValidRose(CasesValidTS,Figdir)
% ����������Ķ��ͳ��ָ�꣬��Bias��RMSE,R,ME,skill����õ��ͼһ��������
% ���뾶Խ������ɫԽ����ð���Խ��
% ע�⣺
%     Bias��RMSEԽСԽ�ã���˴�ʱ�뾶������ɫ��
%     R,ME,skillԽ��Խ�ã���˴�ʱ�뾶������ɫ��


% Ϊÿ����������Ƕ�
fcases = fieldnames(CasesValidTS);
ncases = numel(fcases); % number of cases
ocases = ncases;
switchAdd = 1;
while switchAdd
    if mod(360,ncases) == 0
        switchAdd = 0;
    else
        ncases = ncases + 1;
    end
end
casesAngleCenter = linspace(0,360,ncases+1);
casesAngleCenter = casesAngleCenter(1:end-1);
angleratio=0.9;
eachAngleRange = angleratio*360/ncases;

% ȷ��Ҫ�������򡢱�����ָ��
FigAreas = {'Damian','Lianxu'};
FigVarnames = {'Salt'};
FigIndicators = {'Bias','RMSE','R','ME','skill'};
nSections = numel(FigAreas)*numel(FigVarnames);
eachSectionAngle = eachAngleRange/nSections;

% ��ͼ����
figcolor = 'w';
minCircle = 0.12;
BiasCircle = 0.0471; BiasRange=[0 2]; 
BiasColor=flipud(feval('gray',256)); % ��ɫ�������� --> �С��
RMSECircle = 0.1983; RMSERange=[5 7]; 
RMSEColor=flipud([linspace(0,173/256,256)',linspace(0,235/256,256)',linspace(1,255/256,256)']); % ��ɫ�������� --> �С��
RCircle = 0.2065; RRange=[0.8 1]; 
RColor=flipud([linspace(0,204/256,256)',linspace(153/256,255/256,256)',linspace(51/256,204/256,256)']); % ��ɫ������С�� --> ���
MECircle = 0.2136; MERange=[0.7 0.9]; 
MEColor=flipud([linspace(222/256,255/256,256)',linspace(125/256,255/256,256)',linspace(0/256,204/256,256)']); % ��ɫ������С�� --> ���
skillCircle = 0.2145; skillRange=[0.92 0.97]; 
skillColor=flipud([linspace(255/256,255/256,256)',linspace(0/256,204/256,256)',linspace(0/256,204/256,256)']); % ��ɫ������С�� --> ���

% ��ͼ
figure
plot(0,0,'.','color',figcolor,'markeredgecolor',figcolor,'markerfacecolor',figcolor); % This will create an empty legend entry.
    hold on; axis([-1.01 1.3 -1.01 1.2 ]); axis equal; axis off;

for icase = 1:ocases
    CaseName = fcases{icase};
    theta = casesAngleCenter(icase);
    count = 1; % ��¼�ڸð������ǵڼ�������
    R = [minCircle;1];
    plot(R*sind(theta+(1-angleratio)*180/ncases),R*cosd(theta+(1-angleratio)*180/ncases),':','color','k'); % Draw radial axis
    text(sind(theta- eachAngleRange+180/ncases),cosd(theta- eachAngleRange+180/ncases),CaseName);
    for iarea = 1:numel(FigAreas)
        AreaName = FigAreas{iarea};
        for ivar = 1:numel(FigVarnames)
            varname = FigVarnames{ivar};
            theta_1 = theta - eachAngleRange + eachSectionAngle*(count-1);
            theta_2 = theta - eachAngleRange + eachSectionAngle*(count);
            count=count+1;
            for iInd=1:numel(FigIndicators)
                indicator = FigIndicators{iInd};
                str = sprintf('value = CasesValidTS.%s.%s%s_%s;',CaseName,AreaName,varname,indicator);
                eval(str);
                value=abs(value);
                str=sprintf('tempRange = linspace(min(%sRange),max(%sRange),256);',indicator,indicator);
                eval(str);
                
                % ����Bias��RMSE��ԽСԽ�ã���˴˴���Ϊ��ʹ��Bias��RMSE��R��ME��skill��һ��������Խ��Խ��
                switch indicator
                    case 'Bias'
                        tempRange = fliplr(tempRange);
                    case 'RMSE'
                        tempRange = fliplr(tempRange);
                end
                
                % �ҳ���valueֵ���Ӧ�������ɫ
                DiffValue = abs(value-tempRange);
                [VMin,VIndex] = min(DiffValue);
                
                value=(value-tempRange(1))/(tempRange(end)-tempRange(1));
                % ������ʵvalueֵ�������趨��Χʱ�Ĵ���ͼ
                if value<0
                    value = 0;
                end
                if value>1
                    value =1;
                end
                str=sprintf('value=value* %sCircle;',indicator);
                eval(str)
                str=sprintf('color=%sColor;',indicator);
                eval(str)
                               
                if iInd == 1
                    r1=minCircle; r2=r1+value;
                else
                    r1=r2; r2=r1+value;
                end
                thetaRange=linspace(theta_1,theta_2,100);
                x1 = r1*sind(fliplr(thetaRange)); x2 = r2*sind(thetaRange);
                y1 = r1*cosd(fliplr(thetaRange)); y2 = r2*cosd(thetaRange);
                x=[x1,x2];y=[y1,y2];
                fill(x,y,color(VIndex,:),'edgecolor','w'); 
            end
        end
    end
end

% ����Բ������
[x,y]=cylinder(1,1250); x = x(1,:); y = y(1,:); % Get x and y for a unit-radius circle
circles=[ 0.2 0.4 0.6 0.8 1];
plot(x'*circles,y'*circles,':','color','k'); hold on
plot(x'*circles(end),y'*circles(end),'color','k'); hold on

% ������ɫ��
% Bias
x=linspace(0.85,1.1,256);
y=[1.10 1.15];
for i=1:256
    plot([x(i) x(i)],y,'color',BiasColor(i,:)); hold on
end 
text('String',num2str(BiasRange(2)),'Position',[0.85 1.17],'HorizontalAlignment','left','FontSize',6);
text('String',num2str(BiasRange(1)),'Position',[1.10 1.17],'HorizontalAlignment','right','FontSize',6);
text('String','Bias','Position',[1.15,1.12],'HorizontalAlignment','left','FontSize',8)

% RMSE
x=linspace(0.85,1.1,256);
y=[1.00 1.05];
for i=1:256
    plot([x(i) x(i)],y,'color',RMSEColor(i,:)); hold on
end 
text('String',num2str(RMSERange(2)),'Position',[0.85 1.07],'HorizontalAlignment','left','FontSize',6);
text('String',num2str(RMSERange(1)),'Position',[1.10 1.07],'HorizontalAlignment','right','FontSize',6);
text('String','RMSE','Position',[1.15,1.02],'HorizontalAlignment','left','FontSize',8)

% R
x=linspace(0.85,1.1,256);
y=[0.90 0.95];
for i=1:256
    plot([x(i) x(i)],y,'color',RColor(i,:)); hold on
end 
text('String',num2str(RRange(1)),'Position',[0.85 0.97],'HorizontalAlignment','left','FontSize',6);
text('String',num2str(RRange(2)),'Position',[1.10 0.97],'HorizontalAlignment','right','FontSize',6);
text('String','R^2','Position',[1.15,0.92],'HorizontalAlignment','left','FontSize',8)

% ME
x=linspace(0.85,1.1,256);
y=[0.80 0.85];
for i=1:256
    plot([x(i) x(i)],y,'color',MEColor(i,:)); hold on
end 
text('String',num2str(MERange(1)),'Position',[0.85 0.87],'HorizontalAlignment','left','FontSize',6);
text('String',num2str(MERange(2)),'Position',[1.10 0.87],'HorizontalAlignment','right','FontSize',6);
text('String','ME','Position',[1.15,0.82],'HorizontalAlignment','left','FontSize',8)

% skill
x=linspace(0.85,1.1,256);
y=[0.70 0.75];
for i=1:256
    plot([x(i) x(i)],y,'color',skillColor(i,:)); hold on
end 
text('String',num2str(skillRange(1)),'Position',[0.85 0.77],'HorizontalAlignment','left','FontSize',6);
text('String',num2str(skillRange(2)),'Position',[1.10 0.77],'HorizontalAlignment','right','FontSize',6);
text('String','skill','Position',[1.15,0.72],'HorizontalAlignment','left','FontSize',8)

if ~exist(Figdir)
    mkdir(Figdir)
end
Figname=sprintf('%s/Fig',Figdir);
for ii=1:numel(FigAreas)
    Figname=sprintf('%s_%s',Figname,FigAreas{ii});
end
for ii=1:numel(FigVarnames)
    Figname=sprintf('%s_%s',Figname,FigVarnames{ii});
end
print('-dpng','-r400',Figname)
%%




    
    
    
    
    
    
    
    
    
    
    
    
    
    
