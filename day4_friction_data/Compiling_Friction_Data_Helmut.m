clc
clear variables;
close all;


%% == opening normal force vs. distance number from average file
[avg_fileName,avg_pathName] = uigetfile('*.*','Select average file'); % the user selects the file
avg_filePath = horzcat(avg_pathName,avg_fileName); % total path name is created
[avg_data,avg_header,raw]=xlsread(avg_filePath); % reads data from selected cycle data file
avg_distance = avg_data(:,2);
avg_force = avg_data(:,5);

%% FIGURE 1: Fn vs. cycle
fig1 = figure('Name', 'Fn vs. cycle');
scatter(avg_distance, avg_force)
xlabel('Distance (m)');
ylabel('Fn (N)');


startLoop = input('Start cycle number: ');
endLoop = input('End cycle number: '); 
numFiles = endLoop-startLoop;

%c_info = getCursorInfo(dcm_obj); hold on; % gets x and y values selected x = Position(1) y = Position(2) 


%% == FRICTION_STRUCT: opening folder with all data files == %%
pathName = uigetdir('*.*','Select Data Folder'); % the user selects the folder with the desired files
allFiles = dir(fullfile(pathName,'*.xlsx')); % selecting all the files in the folder with the .xlsx format
allFiles = natsortfiles(allFiles); % sorting the files so they are in numerical order (uses a function that needed to be downloaded)
%endLoop = length(allFiles); % counting the number of files in the folder

%% initializing structurees
friction_struct = struct([]);
friction_crop_struct = struct([]);
mu_avg = zeros(numFiles,5);

%% == looping through all the data to store in a structure == %%
for k = startLoop:endLoop
    %% reading the data from the cycle files
    fileName = allFiles(k).name; % iterating through all the sorted filenames in the folder
    filePath = fullfile(pathName,fileName); % total path name is created
    [my_data,txt,raw]=xlsread(filePath);    % reads data from selected cycle data file


    %% === Read important columns into individual arrays === %% 
    nonnans = find(~isnan(my_data(:,1)));   % finds indices of all numbers in column 1
    startpoint = nonnans(1);            % index of first data point in column 5 that has a number
    clear nonnans

    cycle_time = my_data(startpoint:end,1);           % cycle time [seconds]
    Fn = my_data(startpoint:end,2);                   % normal force [µN]
    aveFn=mean(Fn)                                 %Gather average Fn value
    Ff = my_data(startpoint:end,3);                   % friction force [µN]
    x_stage = my_data(startpoint:end,5);              % x stage displacement [µm]
    %z_stage = my_data(startpoint:end,9);              % z stage displacement [µm]
    mu = Ff./aveFn;                         % friction coefficient (just dividing Ff by Fn)
    
    %% creating a structure to store all the data 
    friction_struct(k).name = fileName;
    friction_struct(k).cycle = k;
    friction_struct(k).x_stage = x_stage;
    friction_struct(k).Ff = Ff;
    friction_struct(k).Fn = Fn;
    friction_struct(k).mu = mu;
    friction_struct(k).data = [x_stage, Ff, Fn, mu];
end

%% == FIGURE 2: setting up Fn vs. x and u vs. x == %%
fig2 = tiledlayout(2,1);
% == Fn vs. x position == %
ax1 = nexttile;
for i = startLoop:endLoop
    plot(ax1, friction_struct(i).x_stage,friction_struct(i).Fn)
    %cycle = [1:numFiles];
    %legend(strcat('Cycle ',num2str(cycle')));
    title('Normal force vs. x position');
    xlabel('x position (mm)');
    ylabel('normal force (N)');
    hold on
end

ax2 = nexttile;
for j = startLoop:endLoop
    scatter(ax2, friction_struct(j).x_stage,friction_struct(j).mu)
    %ylim([-1 1])
    %cycle = [1:numFiles];
    %legend(strcat('Cycle ',num2str(cycle')));
    title('Friction coefficient vs. x position');
    xlabel('x position (mm)');
    ylabel('friction coefficient');
    hold on 
end
disp('Press ''Enter'' when ready to choose data range to analyze');
%pause;

%% == defining region of friction loop to analyze == %%
start_x = input('Start position (mm): '); 
x_range = input('Range (mm): ');
end_x = start_x+x_range; 

%% == FIGURE 3: iterating through mu vs. x positions == %% 
choice = input('Do you want to automatically analyze all cycles (a) or manually choose the cycles (m)?: ', 's');
switch choice
    case 'm'
        for j = startLoop:endLoop
            fig3 = figure('Name', "Cycle " + j);
            scatter(friction_struct(j).x_stage,friction_struct(j).mu)
            set(gcf, 'Position', get(0, 'Screensize')); %makes the graph full screen
            %ylim([-1 1])
            %cycle = [1:numFiles];
            %legend(strcat('Cycle ',num2str(cycle')));
            hold on;
            xline(start_x);
            xline(end_x);    
            title("Cycle " + j);
            xlabel('x position (um)');
            ylabel('friction coefficient'); 
            okay = input('Is this cycle okay (yes = 1 or no = 0)?: ');
            friction_struct(j).okay = okay;  
            close(fig3);
        end
    case 'a'
        for j = startLoop:endLoop
            friction_struct(j).okay = 1;
        end
end     

%% == FRICTION_CROP_STRUCT: calculating friction coefficients == %%
for m = startLoop:endLoop      
    %% finding the index values that are within the stated bounds
    index = find(friction_struct(m).x_stage>start_x & friction_struct(m).x_stage<end_x & friction_struct(m).okay == 1); 
    
    % creating a new structure with just the cropped data
    friction_crop_struct(m).name = friction_struct(m).name;
    friction_crop_struct(m).cycle = friction_struct(m).cycle;
    friction_crop_struct(m).x_stage = friction_struct(m).x_stage(index);
    friction_crop_struct(m).Ff = friction_struct(m).Ff(index);
    friction_crop_struct(m).Fn = friction_struct(m).Fn(index);
    friction_crop_struct(m).mu = friction_struct(m).mu(index);
    
    %% splitting the data into the forward and reverse direction
    [max_value, middle_index] = max(friction_crop_struct(m).x_stage); % finding the friction forces that are positive (forward direction)
    friction_crop_struct(m).Ff_forward = friction_crop_struct(m).Ff(1:middle_index);  
    friction_crop_struct(m).Ff_reverse = friction_crop_struct(m).Ff(middle_index+1:end);
 
    %friction_crop_struct(m).x_stage_forward = friction_crop_struct(m).x_stage(1:middle_index);
    %friction_crop_struct(m).Fn_forward = friction_crop_struct(m).Fn(1:middle_index);
    %friction_crop_struct(m).x_stage_reverse = friction_crop_struct(m).x_stage(middle_index+1:end);
    %friction_crop_struct(m).Fn_reverse = friction_crop_struct(m).Fn(middle_index+1:end);
        
    %% calculating the average and combined uncertainty 
    avg_Ff_forward = mean(friction_crop_struct(m).Ff_forward);
    std_Ff_forward = std(friction_crop_struct(m).Ff_forward);
    avg_Ff_reverse = mean(friction_crop_struct(m).Ff_reverse);
    std_Ff_reverse = std(friction_crop_struct(m).Ff_reverse);
    avg_Ff = (avg_Ff_forward - avg_Ff_reverse)/2;
    avg_Fn = mean(friction_crop_struct(m).Fn);
    std_Fn = std(friction_crop_struct(m).Fn);
    
    mu_avg(m,1) = avg_Ff/avg_Fn; %friction coefficient
    mu_avg(m,2) = ((avg_Ff_forward+std_Ff_forward)-avg_Ff_reverse)./(2*avg_Fn); % friction coefficient with Ff forward error 
    mu_avg(m,3) = (avg_Ff_forward-(avg_Ff_reverse+std_Ff_reverse))./(2*avg_Fn); % friction coefficient with Ff reverse error 
    mu_avg(m,4) = (avg_Ff_forward-avg_Ff_reverse)./(2*(avg_Fn+std_Fn)); % friction coefficient with Fn error
    mu_avg(m,5) = sqrt((mu_avg(m,2)-mu_avg(m,1)).^2+(mu_avg(m,3)-mu_avg(m,1)).^2+(mu_avg(m,4)-mu_avg(m,1)).^2); % combined uncertainty of friction coefficient
    
    %% saving the desired values
    friction_crop_struct(m).Ff_avg = avg_Ff;
    friction_crop_struct(m).Fn_avg = avg_Fn;
    friction_crop_struct(m).mu_avg = mu_avg(m,1); % saving average cycle friction coefficient to structure
    friction_crop_struct(m).mu_cu = mu_avg(m,5);  % saving average cycle friction coefficient combined uncertainty to structure
end

sample_number = 0;
fprintf('\nCycles with all NaN data:\n');
for q = startLoop:endLoop
    if friction_struct(q).okay == 0
        fprintf('%d\n', friction_crop_struct(q).cycle);
    end
    sample_number = sample_number + friction_struct(q).okay;
end

%% == FIGURE 4: plotting friction coefficient vs. x position in the desired region == %
% fig4 = figure('Name', 'Data in the desired range');
% for p = startLoop:numFiles
%     scatter(friction_crop_struct(p).x_stage, friction_crop_struct(p).mu) 
%     title('Friction coeffient vs. x stage')
%     xlim([start_x end_x])
%     %cycle = [1:numFiles];
%     %legend(strcat('Cycle ',num2str(cycle')));
%     xlabel('x stage (um)')
%     ylabel('friction coefficient')
%     hold on
% end

%% == FIGURE 5: plotting friction coefficient vs. cycle == %%
 fig5 = figure('Name', 'Friction coefficient vs. cycle');
 for n = startLoop:numFiles
     scatter(friction_crop_struct(n).cycle, friction_crop_struct(n).mu_avg)
     title('Friction coefficient vs. cycle number')
     xlabel('cycle number')
     ylabel('friction coefficient')
     hold on
 end

%% == calculating average friction coefficient for all cycles == %%
sample_mu_avg = mean([friction_crop_struct.mu_avg],'omitnan');
sample_mu_std = std([friction_crop_struct.mu_avg], 'omitnan');
sample_Fn_avg = mean([friction_crop_struct.Fn_avg], 'omitnan');
sample_Fn_std = std([friction_crop_struct.Fn_avg], 'omitnan');
sample_Ff_avg = mean([friction_crop_struct.Ff_avg], 'omitnan');
sample_Ff_std = std([friction_crop_struct.Ff_avg], 'omitnan');
friction_avg_matrix = [sample_Fn_avg sample_Fn_std sample_Ff_avg sample_Ff_std sample_mu_avg sample_mu_std sample_number];

fprintf('\nSample size: %d\n', sample_number);
fprintf('Sample friction coefficient average: %0.3f\n', sample_mu_avg);
fprintf('Sample friction coefficient std: %0.3f\n', sample_mu_std);
fprintf('Fn_avg: %0.0f +/- %0.0f\n', sample_Fn_avg, sample_Fn_std);


%% == saving structure files so I don't have to do it by hand == %%
splitName = split(avg_fileName,"_");
saveName = input('\nFilename: ', 's');
force = input('Force (N): ', 's');
saveFileName1 = strcat(saveName, '_friction_struct_', force, 'N.mat');
saveFileName2 = strcat(saveName, '_friction_crop_struct_', force, 'N.mat');
saveFileName3 = strcat(saveName, '_friction_avg_matrix_', force, 'N.mat');

save(fullfile(avg_pathName,saveFileName1), 'friction_struct'); 
save(fullfile(avg_pathName,saveFileName2), 'friction_crop_struct');
save(fullfile(avg_pathName,saveFileName3), 'friction_avg_matrix'); 

