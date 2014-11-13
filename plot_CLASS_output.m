function plot_CLASS_output(datafiles,varargin)
%  plot_CLASS_output(datafiles,selection,options)
%  plot_CLASS_output(datafiles,selection)
%  plot_CLASS_output(datafiles,options)
%  plot_CLASS_output(datafiles)
%  Thomas Tram, 12th of March 2014. thomas.tram@epfl.ch
%  Small plot utility for plotting background, thermodynamics and
%  perturbations files. (Compatibility with other files may be added later.)
%  Examples:
%
%    plot_CLASS_output('c:\class\test_background.dat')
%    plots every column of data in the background file
%    'c:\class\test_background.dat' and saves the plot in myplot.eps.
%
%
%    plot_CLASS_output('c:\class\test_background.dat',{'cdm','ur','crit'})
%    plots every mention of either 'cdm', 'ur' or 'crit' in the column titles.
%
%
%    plot_CLASS_output('c:\class\test_perturbations_k0_s.dat',{'delta'})
%    plots every mention of 'delta' in the column titles. Convenient for
%    plotting the density perturbations of all species.
%
%
%    plot_CLASS_output({'c:\class\model1_background.dat','c:\class\model2_background.dat'},{'rho'})
%    plots all mentions of rho in the files listed in the first cell array.
%
%
%    plot_CLASS_output('c:\class\test_perturbations_k0_s.dat',{'delta'},options)
%    options follows the usual MATLAB convention of
%    ...,paramname1,paramval1, paramname2, paramval2,...
%
%    Names:               Values:
%    ----------------------------------------------------------------------
%    EpsFilename          Filename for output
%    xvariable            'a' or 'z' (will convert one from the other)
%    xscale               'log' or 'linear'
%    yscale               'log' or 'linear'
%    xlim                 [xmin xmax]

if mod(nargin,2)==0
    opt = cell2struct(varargin(3:2:end),lower(varargin(2:2:end)),2);
    if ischar(varargin{1})
        opt.selection = varargin(1);
    else
        opt.selection = varargin{1};
    end
else
    opt = cell2struct(varargin(2:2:end),lower(varargin(1:2:end)),2);
    opt.selection = 'all';
end
%Filename for saving plot:
if isfield(opt,'epsfilename');
    epsfilename = opt.epsfilename;
else
    epsfilename = 'myplot';
end

%=================================================================
close all

if ischar(datafiles)
    datafiles = {datafiles};
end


linestyles = {'-','--',':','-.'};
legendcell = {};

xmin = inf;
xmax = -inf;
for fileidx = 1:length(datafiles)
    linestyle = linestyles{mod(fileidx-1,length(linestyles))+1};
    datafile = datafiles{fileidx};
    
    %Find column titles:
    fid = fopen(datafile);
    titleline = 0;
    tline = fgetl(fid);
    while tline(1)=='#'
        titleline = titleline+1;
        tline_old = tline;
        tline = fgetl(fid);
    end
    fclose(fid);
    
%     S=importdata(datafile);
%     data = S.data;
%    titlestring = S.textdata{titleline};
    titlestring = tline_old;
    S = importdata(datafile,' ',titleline);
    data = S.data;
    
    %remove leading #
    titlestring = titlestring(2:end);
    colonidx = [find(titlestring==':'),length(titlestring)+1];
    for j=1:(length(colonidx)-1)
        cellnames{j} = strtrim(titlestring(colonidx(j)+1:colonidx(j+1)-3));
    end
    
    %Determine independent variable:
    x = data(:,1);
    xlab = cellnames{1};
    xscale = 'linear';
    
    if (~isempty(strfind(cellnames{1},'z')))
        xscale = 'log';
        %First column is z
        if (~isfield(opt,'xvariable')) || (~strcmp(opt.xvariable,'z'))
            %if the field is empty or the field is not 'z' change to scale factor a
            x = 1./(x+1);
            xlab = 'a';
        end
    end
    
    if (~isempty(strfind(cellnames{1},'a')))
        xscale = 'log';
        %First column is a
        if (isfield(opt,'xvariable')) && (strcmp(opt.xvariable,'z'))
            %Change to z
            x = 1./x-1;
            xlab = 'z';
        end
    end
    
    if isfield(opt,'xscale')
        xscale = opt.xscale;
    end
    
    
    if strcmp(opt.selection,'all')
        indices = 2:length(cellnames);
    else
        tmp = [];
        for i=1:length(opt.selection)
            tmp2=strfind(cellnames,opt.selection{i});
            for j=1:length(tmp2)
                if ~isempty(tmp2{j})
                    tmp = [tmp,j];
                end
            end
        end
        indices = unique(tmp);
    end
    
    if isempty(indices)
        thenames = '';
        for i=1:length(opt.selection)
            thenames = [thenames,opt.selection{i},', '];
        end
        error(['No indices corresponding to the name(s) {',thenames,'} were found!'])
    end
    
    xmin = min(xmin,min(x));
    xmax = max(xmax,max(x));
    if isfield(opt,'xlim')
        %We need to restrict plotting:
        xl = opt.xlim;
        mask = (x>=xl(1))&(x<=xl(2));
    else
        mask = true(size(x));
    end  
    
    %semilogy(tau,data(:,indices))
    loglog(x(mask),abs(data(mask,indices)),'LineWidth',2,'LineStyle',linestyle)
    
    hold on
    legendcell = [legendcell,cellnames(indices)];
end

legend(legendcell,'Interpreter','none','Location','best')
yl = ylim;
%Fixes x label when y values go to 0 in log plot:
ylim([max(1e-100,yl(1)),yl(2)])
xlabel(xlab)
set(gca,'xscale',xscale);
if isfield(opt,'yscale')
    set(gca,'yscale',opt.yscale);
end
if isfield(opt,'xlim')
    xlim(opt.xlim)
else
    xlim([xmin,xmax])
end

saveas(gcf,[epsfilename,'.eps'],'epsc2')