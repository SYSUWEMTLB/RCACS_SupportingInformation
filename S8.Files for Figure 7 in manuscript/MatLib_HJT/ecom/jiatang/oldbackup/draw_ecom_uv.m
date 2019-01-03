function [M,StopStep] = draw_ecom_uv(FileToRead,GridInfo,Coastal,PlotLevel, ...
    iniStep,endStep,OffsetStep,StartModelTime,varname,FigPropertyOpt)
% function draw_ecom_ts a function to draw scalar variable obtained from
% ECOM simulations
% Created by Jiatang, SYSU

for iStep = iniStep:endStep
    try
        m = memmapfile(...
            FileToRead,   ...
            'Format', {   ...
            'uint32'  [1]  'tag1';
            'single' [GridInfo.nx GridInfo.ny GridInfo.nz] 'u';
            'uint32'  [1]  'tag2';
            'uint32'  [1]  'tag3';
            'single' [GridInfo.nx GridInfo.ny GridInfo.nz] 'v';
            'uint32'  [1]  'tag4';
            %     'uint32'  [1]  'tag5';
            %     'single' [GridInfo.nx GridInfo.ny GridInfo.nz] 'w'; 没w
            %     'uint32'  [1]  'tag6';
            }, ...
            'Repeat',1,  ...
            'offset',(iStep-1)*8*(2+GridInfo.nx*GridInfo.ny*GridInfo.nz) ... %数据块的大小
            ); 
        mdata = m.data;    
    catch ME    
        StopStep = iStep-1;
        return;
    end
    StopStep = iStep;
    DateToPlot = (iStep+OffsetStep)/24.+StartModelTime;
    PlotTitle = sprintf('Level %02d: %s (%s)',PlotLevel,varname,datestr(DateToPlot,21));
    
    fprintf(' -- plotting %s: Level %02d,Step%04d \n',varname, PlotLevel,iStep);
        
    Uvel = mdata(1).u(:,:,:); 
    Vvel = mdata(1).v(:,:,:);
    uu = double(Uvel(:,:,PlotLevel));
    vv = double(Vvel(:,:,PlotLevel));
    % uuu = uu.*cos(GridInfo.angle)-vv.*sin(GridInfo.angle);
    % vvv = uu.*sin(GridInfo.angle)+vv.*cos(GridInfo.angle);
    %
    Vel = uu+sqrt(-1)*vv;
    % rotate into true east/north coordinates
    Vel0 = Vel.*exp(sqrt(-1)*GridInfo.angle);
    % Vel0 = uuu+sqrt(-1)*vvv;

    %
    figure('visible','off');    
    % FillLandMap(Coastal.lon,Coastal.lat);   % 画岸线及其陆地填充
    psliceuv(GridInfo.lon,GridInfo.lat,Vel0,1,0.02,'b');
    title(PlotTitle);
    xlabel('Latitude(degree)');
    ylabel('Longitude(degree)');
%     xlim([112.3,115.7]);
%     ylim([20.8,23]);
    xlim([113,114]);
    ylim([22,23]);    
    %
    amerc;    
    %
    M(iStep) = getframe(gcf); % movie
    
    FigDir = sprintf('Level%02d',PlotLevel);
    if ~exist(FigDir,'dir')
        mkdir(FigDir);
    end
    str = sprintf('%s/Small%s_level%02d_step%02d',FigDir,varname,PlotLevel,iStep+OffsetStep);
    print('-dpng','-r400',str);
    close
end

return
end