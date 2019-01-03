function [M,StopStep] = draw_ecom_ts(FileToRead,GridInfo,Coastal,PlotLevel, ...
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
            'single' [GridInfo.nx GridInfo.ny GridInfo.nz] 'TS';
            'uint32'  [1]  'tag2';
            }, ...
            'Repeat',1,  ...
            'offset',(iStep-1)*4*(2+GridInfo.nx*GridInfo.ny*GridInfo.nz) ... %���ݿ�Ĵ�С
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
    
    TS = mdata(1).TS(:,:,:);       
    Concs = double(TS(:,:,PlotLevel));
    %
    figure('visible','off');
    pcolorjw(GridInfo.lon,GridInfo.lat,Concs);
    % contourf(GridInfo.lon,GridInfo.lat,double(Concs),150);
    % shading interp;
    FillLandMap(Coastal.lon,Coastal.lat);   % �����߼���½�����
    colorbar;
    caxis(FigPropertyOpt.caxis);
    title(PlotTitle);
    xlabel('Latitude(degree)');
    ylabel('Longitude(degree)');
%     xlim([112.3,115.7]);
%     ylim([20.8,23]);
    %
    amerc;    
    %
    M(iStep) = getframe(gcf); % movie
    
    FigDir = sprintf('Level%02d',PlotLevel);
    if ~exist(FigDir,'dir')
        mkdir(FigDir);
    end
    str = sprintf('%s/%s_level%02d_step%02d',FigDir,varname,PlotLevel,iStep+OffsetStep);
    print('-dpng','-r400',str);
    close
end

return
end