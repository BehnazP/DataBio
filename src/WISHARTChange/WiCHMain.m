function WiCHMain(files,npols,pols_name,times_name,Nol,Pvalue,times,flag_size,flag_area,file_figure,file_table,shape_file)
% WICHMAIN or Wichart Change main function identify the changes in a sequence of SAR images
%
% WiCHMain(files,pols_name,times_name,Nol,Pvalue,times,file_figure,file_table,shape_file)
% Input:
% files         a cell structure containing all the images path and name
%               Example:
%               {'C:\Download\201712.VV_0'},{'C:\Download\201712.VV_1'}
%               {'C:\Download\201712.VH_0'},{'C:\Download\201712.VV_1'} and
%               so on, the program support geotiff, envi and emi images
% npols         a number represent the polarization
% pols_name     a char vector containing the name of polarization appear
%               after the name of image following a dot
%               Example: .VV .VH
% times_name    a vector containing either a sequence of number 0:10 or a
%               personalize numbers such as 062, 064, 067 at the end of the
%               name of the files
% Nol           number of looks
% Pvalue        a significant threshold for p-test
% times         number of times to be compared
% flag_size     an indicator; if 0 save the file in workspace if 1 read
%               line by line
%               default value 0
% flag_area     an indicator; if 1 read the flag_area area of the image, if 0 read the ROI
%               default value 1
% file_figure   name of the figure to be saved
% file_table    name of the table to be saved
% shape_file    the function accept either 0 for the case of choosing the ROI interactively or
%               [path and name] of the shape file to be read as ROI
%
% Return a table of average p-values in .csv format, a plot of first, last and
%        frequency of the changes both as pdf and 3 bands images with the same
%        format as original images.
%
% Behnaz Pirzamanbein
% bepi@dtu.dk
% Image Analysis and Computer Graphics section
% Department of Applied Mathematics and Computer Science
% Technical University of Denmark
% First version 17 July 2018
% last version 03-11-2018

p = gcp('nocreate');
if isempty(p)
    parpool;
end

disp('/////////////////////////////////////////////////');
disp('/////////////////////////////////////////////////');
disp('. . . changeDetection: FUNCTION IN PROGRESS . . .');
disp('/////////////////////////////////////////////////');

if ~iscell(files)
    tmp = textscan(files,'%s','Delimiter',',')';
    files = tmp{:}';
end

if ~iscell(pols_name)
    tmp = textscan(pols_name,'%s','Delimiter',',')';
    pols_name = tmp{:}';
end

if ~iscell(times_name)
    tmp = textscan(times_name,'%s','Delimiter',',')';
    times_name = tmp{:}';
end

if ~isnumeric(Nol)
    Nol = str2double(Nol);
end

if ~isnumeric(npols)
    npols = str2double(npols);
end
if ~isnumeric(Pvalue)
    Pvalue = str2double(Pvalue);
end

if ~isnumeric(times)
    times = str2double(times);
end

if nargin < 8
    flag_size = 0;
    flag_area = 1;
    path_name = fullfile(pwd,'/');
    file_figure = fullfile(path_name,'figure');
    file_table = fullfile(path_name,'table');
    shape_file = [];
else
    if ~isnumeric(flag_size)
        flag_size = str2double(flag_size);
    end
    if ~isnumeric(flag_area)
        flag_area = str2double(flag_area);
    end

    if ~ischar(file_figure)
        if isempty(file_figure)
            path = fullfile(pwd,'/');
            file_figure = fullfile(path,'figure');
        else
            tmp = textscan(file_figure,'%s')';
            file_figure = tmp{:}';
        end
    end

    if ~ischar(file_table)
        if isempty(file_table)
            path = fullfile(pwd,'/');
            file_table = fullfile(path,'table');
        else
            tmp = textscan(file_table,'%s')';
            file_table = tmp{:}';
        end
    end

    if ~isempty(shape_file)
        if ~str2double(shape_file)
            shape_file = [];
        else
            tmp = textscan(shape_file,'%s')';
            shape_file = tmp{:}';
        end
    end
end
switch flag_size
    case 0
        [images, ~] = read_data(files,pols_name,times_name);
        wc = WiCHParallel(images,Nol,npols);
        When = WiCHWhen(wc,Pvalue);
        % plot the figure
        fig = figure('Name','Results','NumberTitle','off');
        A = PlotWhen(fig,When,times);
        saveWhen(A,file_figure,files)

        switch flag_area
            case 1
                Prob = WiCHProb(wc);
                % print the table
                Ptable([],Prob,times,file_table);
            case 0
                ROI  = readROI(shape_file,fig);
                Prob_ROI = WiCHProb(wc,ROI.mask);
                % print the table
                file_table_ROI = [file_table,'_ROI'];
                Ptable([],Prob_ROI,times,file_table_ROI);
                plotROI(fig,ROI)
                file_figure = [file_figure,'_ROI'];
        end
        save_fig(fig,file_figure,'landscape');
    case 1
        [When, Prob] = WiCHLine(files,pols_name,times_name,Nol,Pvalue,npols);
        % plot the figure
        fig = figure('Name','Results','NumberTitle','off');
        A = PlotWhen(fig,When,times);
        saveWhen(A,file_figure,files)

        switch flag_area
            case 1
                % print the table
                Ptable([],Prob,times,file_table);
            case 0
                ROI  = readROI(shape_file,fig);
                [~, Prob_ROI] = WiCHLine(files,pols_name,times_name,Nol,Pvalue,npols,ROI.mask);

                % print the table
                file_table_ROI = [file_table,'_ROI'];
                Ptable([],Prob_ROI,times,file_table_ROI);
                plotROI(fig,ROI)
                file_figure = [file_figure,'_ROI'];
        end
        save_fig(fig,file_figure,'landscape');
end
disp('/////////////////////////////////////////////////////////');
disp('. . . changeDetection - WishartChange: PROCESS OVER . . .');
disp('/////////////////////////////////////////////////////////');
end
