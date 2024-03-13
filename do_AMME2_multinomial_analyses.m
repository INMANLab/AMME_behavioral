%Script to do multinomial regression on extended AMME dataset
%JRM 4/28/23
%Martina modified original code 05/09/2023
% The script incudes exploratory analyses with various predictor sets and various models
%ANALYSIS 1: quartile-based responder group, nominal multinomial regression
%ANALYSIS 2: quartile-based responder group, nominal multinomial regression no composite memory score, just RAVLT

%very difficult to understand the output from this model
%ANALYSIS 3: quartile-based responder group, ordinal multinomial regression
%ANALYSIS 4: quartile-based responder group, ordinal multinomial regression
%no composite memory score, just RAVLT


%LET'S BEGIN!

clear
close all

%Read in csv file as a data table
dt = readtable('foranalysis_responder_firstsession_quartilebased.csv');
numsubs = height(dt);%number of subjects

%Let's Z-score the two memory scores and average them for an overall memory score
dt.RAVLT_Z = normalize(dt.RAVLT_dprime);
dt.ReyO_Z = normalize(dt.ReyO_CF_delayedrecall);
%dt.WMSIV_VR_recog_Z = normalize(dt.WMSIV_VR_recog_score);
%WMSIV_VR_recog_Z = dt.WMSIV_VR_recog_Z;
%dt.Memory_Z = rmmissing(dt.WMSIV_VR_recog_Z);
dt.Memory_Z = nanmean([dt.RAVLT_Z dt.ReyO_Z], 2);
%dt.Memory = nanmean([dt.RAVLT_dprime dt.ReyO_CF_delayedrecall], 2);%nanmean because two patients are missing one memory score

%Let's also Z score the IQ scores so that the Beta scaling better aligns with other variable
dt.IQ_Z = normalize(dt.FSIQ);

%make sex be a scalar var: female = 1; male = 0;
dt.sex_num = strcmp(dt.sex, 'female');


%IED activity 
dt.IEDday1 = strcmp(dt.IED_duringimageday1, 'yes'); %1 for yes; 0 fr no
dt.IEDday2 = strcmp(dt.IED_duringimageday2, 'yes'); %1 for yes; 0 fr no
dt.IEDbilateral = strcmp(dt.IED_bilateral, 'yes'); %1 for yes; 0 fr no


%make Experiment names be numbers too
%A better approach, if it's possible, would be to make
% the experiment type be a categorical variable.
% % dt.Exp_nums = zeros(numsubs,1);
% % Exp_names = unique(dt.Experiment);
% % Exp_names = Exp_names([2,1,3]);%reorder so it's now orig, duration, timing
% % for ex = 1:length(Exp_names)%lenght is 3
% %   whichrows = strcmp(dt.Experiment, Exp_names{ex});%compares the exp name in the loop with the exp name in the raw data
% %   dt.Exp_nums(whichrows) = ex; %boolean var whichrows makes the ones (Trues) into the loop number (1, 2, or 3).
% % end%ex

%hemisphere of stimulation
%%dt.LeftStim = strcmp(dt.stim_hemisphere, 'L');%1 for left; 0 for right or bilateral

% % Run correlation between distance and dprime at 1-day test
% [r, p] = corrcoef(dt.dist_BLAtoHPC, dt.avg_stim_dprime_diff,'rows','complete');
% scatter(dt.dist_BLAtoHPC, dt.avg_stim_dprime_diff, 'filled')
% xlabel('Distance to HPC')
% ylabel('dprime at 1-day test')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BEGIN ANALYSIS 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%OUTCOME VARIABLE: QUARTILE-BASED NOMINAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     Predictors with composite memory score                                                      %
%                              You can change these two lines and the rest of the code will adapt                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FYI the first row of output variables from mnrfit will be the intercept and subsequent rows will be our predictor names
% Thus, we always have to include Intercept as the first output variable
% We'll combine the predictors into one matrix

%IQ, Memory, Sex, Dist, age
 my_predictors = [dt.IQ_Z dt.Memory_Z dt.sex_num dt.age dt.shortest_dist_MTL dt.IED_freq];
 my_output_var_names = {'Intercept'; 'IQ'; 'Baseline Memory'; 'Sex'; 'Age'; 'Shortest distance to MTL'; 'IED frequency'};%in order we determined above

%Hemisphere, Experiment, IQ, Memory, Sex, Dist
 %my_predictors = [dt.LeftStim dt.Exp_nums dt.IQ_Z dt.Memory_Z dt.sex_num dt.dist_BLAtoHPC];
 %my_output_var_names = {'Intercept'; 'Hemisphere'; 'Experiment'; 'IQ'; 'Memory'; 'Sex'; 'Distance to HPC'};%in order we determined above

%Hemisphere, IQ, Memory, Sex, Dist, age
 %my_predictors = [dt.LeftStim dt.IQ_Z dt.Memory_Z dt.sex_num dt.dist_BLAtoHPC];
 %my_output_var_names = {'Intercept'; 'Hemisphere'; 'IQ'; 'Memory'; 'Sex'; 'Distance to HPC'};%in order we determined above

%Experiment, IQ, Memory, Sex, Dist, age
 %my_predictors = [dt.Exp_nums dt.IQ_Z dt.Memory_Z dt.sex_num dt.dist_BLAtoHPC dt.age];
 %my_output_var_names = {'Intercept'; 'Experiment'; 'IQ'; 'Memory'; 'Sex'; 'Distance to HPC'; 'Age'};%in order we determined above

%Hemisphere, Memory, Sex, Dist, age
 %my_predictors = [dt.LeftStim dt.Memory_Z dt.sex_num dt.dist_BLAtoHPC];
 %my_output_var_names = {'Intercept'; 'Hemisphere'; 'Memory'; 'Sex'; 'Distance to HPC'};%in order we determined above

%Experiment, Memory, Sex, Dist, age
 %my_predictors = [dt.Exp_nums dt.Memory_Z dt.sex_num dt.dist_BLAtoHPC];
 %my_output_var_names = {'Intercept'; 'Experiment'; 'Memory'; 'Sex'; 'Distance to HPC'};%in order we determined above


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF PREDICTORS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

my_response_categories = categorical(dt.avg_responder_status);
%reorder to make non-responders last since last is reference category
my_response_categories = reordercats(my_response_categories, ...
  {'strong responder';'moderate responder'; 'anti-responder'; 'non-responder'});
my_response_names = categories(my_response_categories);
    % {'strong responder'  }  %SR
    % {'moderate responder'}  %MR
    % {'anti-responder'    }  %AR
    % {'non-responder'     }  %NR

%all the stats from mnrfit will be about odds relative to reference condition, which is Non-responder here
my_comparison_names = {'AN', 'MR', 'SR'}; %always compared to Non-responders
%let's label the rows of the output stats too


%use Matlab function mnrfit to do multinomial regression
%stats is a structure that contains all the info
%interactions are ON for this model, default for nominal models is ON
[~,dev,stats] = mnrfit(my_predictors, my_response_categories, 'Model', 'nominal', 'Interactions', 'on');

% Calculate AIC manually
nobs = size(my_predictors,2); % number of observations
k = size(stats.beta,1); % number of estimated parameters
L = -dev; % log-likelihood
aic_value = 2*k - 2*L + (2*k*(k+1))/(nobs-k-1);

%Calculate Goodness-of-fit which compares the given model to its intercept mode
gof = dev/stats.dfe; %uses the degrees of freedom estimate from the model

%the beta values are the "log odds of being in one category versus the reference category" per documentation.
my_beta_table = array2table(stats.beta, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
my_beta_stderr_table = array2table(stats.se, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
my_beta_t_table = array2table(stats.t, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
my_pval_table = array2table(stats.p, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);


%LET'S PLOT
figure('Name', 'Log Odds of Being Strong- vs Non-Responder', 'Color', 'w', 'WindowState', 'Maximized');

%for conveneience, let's plot a suplot for each predictor
num_predictors = size(my_predictors, 2); %take the second dimension of size (columns)
num_comparisons = length(my_comparison_names);

clr = [255,0,0; 
      128,0,128; 
      255,0,255]/255;

my_axes = NaN(1,num_predictors);%axes handles for linking Y/X axes

for sp = 1:num_predictors
  my_axes(sp) = subplot(1,num_predictors, sp);
  hold on
  tmp_yvals = my_beta_table{sp+1,:};%add 1 because first row is intercept
  tmp_errvals = my_beta_stderr_table{sp+1,:};%add 1 because first row is intercept
  tmp_pvals =  my_pval_table{sp+1,:};%add 1 because first row is intercept
  b = bar(1:num_comparisons, tmp_yvals, 1, 'FaceColor', 'flat', 'BarWidth',1);
  b.CData = clr;
  errorbar(1:num_comparisons, tmp_yvals, tmp_errvals,'k', 'LineStyle', 'none', 'LineWidth', 1)
  xticks(1:num_comparisons)
  xticklabels(my_comparison_names)
  set(gca, 'FontName', 'arial', 'FontSize', 15, 'LineWidth', 1)
  xlabel('Responder Comparison', 'FontWeight', 'bold', 'FontSize',15)
  ylabel('Log Odds Relative to Non-Responder (vNR)', 'FontWeight', 'bold', 'FontSize',20)
  title(my_output_var_names{sp+1},'FontSize', 15)%add 1 because first row is intercept
 
end%subplot loop

%make y axes be same across subplots
linkaxes(my_axes, 'xy')

my_ylim = get(gca, 'YLim');
tmp_y = my_ylim(1)+.4;%text just above bottom of figure

%let's go back through to put p values on figure above bars
for sp = 1:num_predictors
  subplot(my_axes(sp)), hold on
  for b = 1:num_comparisons
    tmpp = my_pval_table{sp+1,b};% p value for this bar
    text(b,tmp_y, sprintf('p =\n%.3f', tmpp), 'HorizontalAlignment', 'center', 'FontSize', 12)
  end%bar
end%subplot
%End of first figure




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF ANALYSIS 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%make sure to comment out previous analysis sections bc we reuse variable names


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN ANALYSIS 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%OUTCOME VARIABLE: QUARTILE-BASED NOMINAL WITH RAVLT ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                                     Predictors with RAVLT only                                                     %
% %                      You can change these two lines and the rest of the code will adapt                               %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % FYI the first row of output variables from mnrfit will be the intercept and subsequent rows will be our predictor names
% % Thus, we always have to include Intercept as the first output variable
% % We'll combine the predictors into one matrix
% 
% %Hemisphere, Experiment, IQ, RAVLT, Sex
% my_predictors = [dt.LeftStim dt.Exp_nums dt.IQ_Z dt.RAVLT_Z dt.sex_num];
% my_output_var_names = {'Intercept'; 'Hemisphere'; 'Experiment'; 'IQ'; 'RAVLT'; 'Sex'};%in order we determined above
% 
% %Hemisphere, IQ, RAVLT, Sex
% % my_predictors = [dt.LeftStim dt.IQ_Z dt.RAVLT_Z dt.sex_num];
% % my_output_var_names = {'Intercept'; 'Hemisphere'; 'IQ'; 'RAVLT'; 'Sex'};%in order we determined above
% 
% %Experiment, IQ, RAVLT, Sex
% %my_predictors = [dt.Exp_nums dt.IQ_Z dt.RAVLT_Z dt.sex_num];
% %my_output_var_names = {'Intercept'; 'Experiment'; 'IQ'; 'RAVLT'; 'Sex'};%in order we determined above
%  
% % Hemisphere,  RAVLT, Sex
% % my_predictors = [dt.LeftStim dt.RAVLT_Z dt.sex_num];
% % my_output_var_names = {'Intercept'; 'Hemisphere'; 'RAVLT'; 'Sex'};%in order we determined above
% 
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF PREDICTORS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% my_response_categories = categorical(dt.responder_status);
% %reorder to make non-responders last since last is reference category
% my_response_categories = reordercats(my_response_categories, ...
%   {'strong responder'; 'moderate responder'; 'anti-responder'; 'non-responder'});
% my_response_names = categories(my_response_categories);
% %     {'strong responder'  }  %SR
% %     {'moderate responder'}  %MR
% %     {'anti-responder'    }  %AR
% %     {'non-responder'     }  %NR
% 
% %all the stats from mnrfit will be about odds relative to reference condition, which is Non-responder here
% my_comparison_names = {'SRvNR', 'MRvNR', 'ARvNR'};
% %let's label the rows of the output stats too
% 
% 
% %use Matlab function mnrfit to do multinomial regression
% %stats is a structure that contains all the info
% %interactions are ON for this model, default for nominal models is ON
% [~,~,dev,stats] = mnrfit(my_predictors, my_response_categories, 'Model', 'nominal', 'Interactions', 'on');

% % Calculate AIC manually
% nobs = size(my_predictors,2); % number of observations
% k = size(stats.beta,1); % number of estimated parameters
% L = -dev; % log-likelihood
% aic_value = 2*k - 2*L + (2*k*(k+1))/(nobs-k-1);
% 
% %Calculate Goodness-of-fit which compares the given model to its intercept mode
% gof = dev/stats.dfe; %uses the degrees of freedom estimate from the model

% 
% %the beta values are the "log odds of being in one category versus the reference category" per documentation.
% my_beta_table = array2table(stats.beta, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% my_beta_stderr_table = array2table(stats.se, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% my_beta_t_table = array2table(stats.t, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% my_pval_table = array2table(stats.p, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% 
% 
% %LET'S PLOT
% figure('Name', 'Log Odds of Being Strong-, Moderate-, or Anti- vs Non-Responder', 'Color', 'w', 'WindowState', 'Maximized')
% 
% %for conveneience, let's plot a suplot for each predictor
% num_predictors = size(my_predictors, 2); %take the second dimension of size (columns)
% num_comparisons = length(my_comparison_names);
% 
% my_axes = NaN(1,num_predictors);%axes handles for linking Y/X axes
% 
% for sp = 1:num_predictors
%   my_axes(sp) = subplot(1,num_predictors, sp);
%   hold on
%   tmp_yvals = my_beta_table{sp+1,:};%add 1 because first row is intercept
%   tmp_errvals = my_beta_stderr_table{sp+1,:};%add 1 because first row is intercept
%   tmp_pvals =  my_pval_table{sp+1,:};%add 1 because first row is intercept
%   bar(1:num_comparisons, tmp_yvals, 1, 'FaceColor', [.6 .6 .6])
%   errorbar(1:num_comparisons, tmp_yvals, tmp_errvals,'k', 'LineStyle', 'none', 'LineWidth', 1)
%   xticks(1:num_comparisons)
%   xticklabels(my_comparison_names)
%   set(gca, 'FontName', 'arial', 'FontSize', 14, 'LineWidth', 1)
%   xlabel('Responder Comparison', 'FontWeight', 'bold', 'FontSize',17)
%   ylabel('Log Odds Relative to Non-Responder (vNR)', 'FontWeight', 'bold', 'FontSize',17)
%   title(my_output_var_names{sp+1},'FontSize', 20)%add 1 because first row is intercept
%  
% end%subplot loop
% 
% %make y axes be same across subplots
% linkaxes(my_axes, 'xy')
% 
% my_ylim = get(gca, 'YLim');
% tmp_y = my_ylim(1)+.4;%text just above bottom of figure
% 
% %let's go back through to put p values on figure above bars
% for sp = 1:num_predictors
%   subplot(my_axes(sp)), hold on
%   for b = 1:num_comparisons
%     tmpp = my_pval_table{sp+1,b};% p value for this bar
%     text(b,tmp_y, sprintf('p =\n%.3f', tmpp), 'HorizontalAlignment', 'center','FontSize', 15)
%   end%bar
%  end%subplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF ANALYSIS 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN ANALYSIS 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%OUTCOME VARIABLE: QUARTILE-BASED ORDINAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                      Predictors with composite memory score                                                      %
%                               You can change these two lines and the rest of the code will adapt                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FYI the first row of output variables from mnrfit will be the intercept and subsequent rows will be our predictor names
% Thus, we always have to include Intercept as the first output variable
%We'll combine the predictors into one matrix

%Hemisphere, Experiment, IQ, Memory, Sex
% my_predictors = [dt.LeftStim dt.Exp_nums dt.IQ_Z dt.Memory_Z dt.sex_num];
% my_output_var_names = {'Intercept'; 'Hemisphere'; 'Experiment'; 'IQ'; 'Memory'; 'Sex'};%in order we determined above

% %Hemisphere, IQ, Memory, Sex
% my_predictors = [dt.LeftStim dt.IQ_Z dt.Memory_Z dt.sex_num];
% my_output_var_names = {'Intercept'; 'Hemisphere'; 'IQ'; 'Memory'; 'Sex'};%in order we determined above

% %Experiment, IQ, Memory, Sex
% my_predictors = [dt.Exp_nums dt.IQ_Z dt.Memory_Z dt.sex_num];
% my_output_var_names = {'Intercept'; 'Experiment'; 'IQ'; 'Memory'; 'Sex'};%in order we determined above
 
% % Hemisphere,  Memory, Sex
% my_predictors = [dt.LeftStim dt.Memory_Z dt.sex_num];
% my_output_var_names = {'Intercept'; 'Hemisphere'; 'Memory'; 'Sex'};%in order we determined above


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF PREDICTORS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% my_response_categories = categorical(dt.responder_status, 'Ordinal', true);
% %reorder to make non-responders last since last is reference category
% my_response_categories = reordercats(my_response_categories, ...
%   {'strong responder'; 'moderate responder'; 'anti-responder'; 'non-responder'});
% my_response_names = categories(my_response_categories);
% %     {'strong responder'  }  %SR
% %     {'moderate responder'}  %MR
% %     {'anti-responder'    }  %AR
% %     {'non-responder'     }  %NR
% 
% %all the stats from mnrfit will be about odds relative to reference condition, which is Non-responder here
% my_comparison_names = {'SRvNR', 'MRvNR', 'ARvNR'};
% %let's label the rows of the output stats too
% 
% 
% %use Matlab function mnrfit to do multinomial regression
% %stats is a structure that contains all the info
% %interactions are OFF for this model, default for ordinal
% [~,~,dev,stats] = mnrfit(my_predictors, my_response_categories, 'Model', 'ordinal', 'Interactions', 'off');
% 
% 
% %THE CODE BREAKS HERE BC THE OUTPUT IS FORMATTED DIFFERENTLY THAT I CAN'T
% %INTERPRET THE OUTPUT AND WHAT MEANS WHAT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% 
% %the beta values are the "log odds of being in one category versus the reference category" per documentation.
% my_beta_table = array2table(stats.beta, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% my_beta_stderr_table = array2table(stats.se, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% my_beta_t_table = array2table(stats.t, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% my_pval_table = array2table(stats.p, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% 
% 
% 
% %LET'S PLOT
% figure('Name', 'Log Odds of Being Strong-, Moderate-, or Anti- vs Non-Responder', 'Color', 'w', 'WindowState', 'Maximized')
% 
% %for conveneience, let's plot a suplot for each predictor
% num_predictors = size(my_predictors, 2); %take the second dimension of size (columns)
% num_comparisons = length(my_comparison_names);
% 
% my_axes = NaN(1,num_predictors);%axes handles for linking Y/X axes
% 
% for sp = 1:num_predictors
%   my_axes(sp) = subplot(1,num_predictors, sp);
%   hold on
%   tmp_yvals = my_beta_table{sp+1,:};%add 1 because first row is intercept
%   tmp_errvals = my_beta_stderr_table{sp+1,:};%add 1 because first row is intercept
%   tmp_pvals =  my_pval_table{sp+1,:};%add 1 because first row is intercept
%   bar(1:num_comparisons, tmp_yvals, 1, 'FaceColor', [.6 .6 .6])
%   errorbar(1:num_comparisons, tmp_yvals, tmp_errvals,'k', 'LineStyle', 'none', 'LineWidth', 1)
%   xticks(1:num_comparisons)
%   xticklabels(my_comparison_names)
%   set(gca, 'FontName', 'arial', 'FontSize', 14, 'LineWidth', 1)
%   xlabel('Responder Comparison', 'FontWeight', 'bold', 'FontSize',17)
%   ylabel('Log Odds Relative to Non-Responder (vNR)', 'FontWeight', 'bold', 'FontSize',17)
%   title(my_output_var_names{sp+1},'FontSize', 20)%add 1 because first row is intercept
%  
% end%subplot loop
% 
% %make y axes be same across subplots
% linkaxes(my_axes, 'xy')
% 
% my_ylim = get(gca, 'YLim');
% tmp_y = my_ylim(1)+.4;%text just above bottom of figure
% 
% %let's go back through to put p values on figure above bars
% for sp = 1:num_predictors
%   subplot(my_axes(sp)), hold on
%   for b = 1:num_comparisons
%     tmpp = my_pval_table{sp+1,b};% p value for this bar
%     text(b,tmp_y, sprintf('p =\n%.3f', tmpp), 'HorizontalAlignment', 'center','FontSize', 15)
%   end%bar
% end%subplot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF ANALYSIS 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN ANALYSIS 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%OUTCOME VARIABLE: QUARTILE-BASED ORDINAL WITH RAVLT ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     Predictors with RAVLT only                                                     %
%                      You can change these two lines and the rest of the code will adapt                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FYI the first row of output variables from mnrfit will be the intercept and subsequent rows will be our predictor names
% Thus, we always have to include Intercept as the first output variable
% We'll combine the predictors into one matrix

%Hemisphere, Experiment, IQ, RAVLT, Sex
% my_predictors = [dt.LeftStim dt.Exp_nums dt.IQ_Z dt.RAVLT_Z dt.sex_num];
% my_output_var_names = {'Intercept'; 'Hemisphere'; 'Experiment'; 'IQ'; 'RAVLT'; 'Sex'};%in order we determined above

%Hemisphere, IQ, RAVLT, Sex
% my_predictors = [dt.LeftStim dt.IQ_Z dt.RAVLT_Z dt.sex_num];
% my_output_var_names = {'Intercept'; 'Hemisphere'; 'IQ'; 'RAVLT'; 'Sex'};%in order we determined above

%Experiment, IQ, RAVLT, Sex
%my_predictors = [dt.Exp_nums dt.IQ_Z dt.RAVLT_Z dt.sex_num];
%my_output_var_names = {'Intercept'; 'Experiment'; 'IQ'; 'RAVLT'; 'Sex'};%in order we determined above
 
% Hemisphere,  RAVLT, Sex
% my_predictors = [dt.LeftStim dt.RAVLT_Z dt.sex_num];
% my_output_var_names = {'Intercept'; 'Hemisphere'; 'RAVLT'; 'Sex'};%in order we determined above

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF PREDICTORS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% my_response_categories = categorical(dt.responder_status, 'Ordinal', true);
% %reorder to make non-responders last since last is reference category
% my_response_categories = reordercats(my_response_categories, ...
%   {'strong responder'; 'moderate responder'; 'anti-responder'; 'non-responder'});
% my_response_names = categories(my_response_categories);
% %     {'strong responder'  }  %SR
% %     {'moderate responder'}  %MR
% %     {'anti-responder'    }  %AR
% %     {'non-responder'     }  %NR
% 
% %all the stats from mnrfit will be about odds relative to reference condition, which is Non-responder here
% my_comparison_names = {'SRvNR', 'MRvNR', 'ARvNR'};
% %let's label the rows of the output stats too
% 
% 
% %use Matlab function mnrfit to do multinomial regression
% %stats is a structure that contains all the info
% %interactions are OFF for this model, default for ordinal
% [~,~,dev,stats] = mnrfit(my_predictors, my_response_categories, 'Model', 'ordinal', 'Interactions', 'off');
% 
% 
% %THE CODE BREAKS HERE BC THE OUTPUT IS FORMATTED DIFFERENTLY THAT I CAN'T
% %INTERPRET THE OUTPUT AND WHAT MEANS WHAT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% 
% %the beta values are the "log odds of being in one category versus the reference category" per documentation.
% my_beta_table = array2table(stats.beta, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% my_beta_stderr_table = array2table(stats.se, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% my_beta_t_table = array2table(stats.t, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% my_pval_table = array2table(stats.p, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% 
% %LET'S PLOT
% figure('Name', 'Log Odds of Being Strong-, Moderate-, or Anti- vs Non-Responder', 'Color', 'w', 'WindowState', 'Maximized')
% 
% %for conveneience, let's plot a suplot for each predictor
% num_predictors = size(my_predictors, 2); %take the second dimension of size (columns)
% num_comparisons = length(my_comparison_names);
% 
% my_axes = NaN(1,num_predictors);%axes handles for linking Y/X axes
% 
% for sp = 1:num_predictors
%   my_axes(sp) = subplot(1,num_predictors, sp);
%   hold on
%   tmp_yvals = my_beta_table{sp+1,:};%add 1 because first row is intercept
%   tmp_errvals = my_beta_stderr_table{sp+1,:};%add 1 because first row is intercept
%   tmp_pvals =  my_pval_table{sp+1,:};%add 1 because first row is intercept
%   bar(1:num_comparisons, tmp_yvals, 1, 'FaceColor', [.6 .6 .6])
%   errorbar(1:num_comparisons, tmp_yvals, tmp_errvals,'k', 'LineStyle', 'none', 'LineWidth', 1)
%   xticks(1:num_comparisons)
%   xticklabels(my_comparison_names)
%   set(gca, 'FontName', 'arial', 'FontSize', 14, 'LineWidth', 1)
%   xlabel('Responder Comparison', 'FontWeight', 'bold', 'FontSize',17)
%   ylabel('Log Odds Relative to Non-Responder (vNR)', 'FontWeight', 'bold', 'FontSize',17)
%   title(my_output_var_names{sp+1},'FontSize', 20)%add 1 because first row is intercept
%  
% end%subplot loop
% 
% %make y axes be same across subplots
% linkaxes(my_axes, 'xy')
% 
% my_ylim = get(gca, 'YLim');
% tmp_y = my_ylim(1)+.4;%text just above bottom of figure
% 
% %let's go back through to put p values on figure above bars
% for sp = 1:num_predictors
%   subplot(my_axes(sp)), hold on
%   for b = 1:num_comparisons
%     tmpp = my_pval_table{sp+1,b};% p value for this bar
%     text(b,tmp_y, sprintf('p =\n%.3f', tmpp), 'HorizontalAlignment', 'center','FontSize', 15)
%   end%bar
% end%subplot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF ANALYSIS 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% %make sex be a scalar var: female = 1; male = 0;
% dt.sex_num = strcmp(dt.sex, 'female');
% 
% %make Experiment names be numbers too
% %A better approach, if it's possible, would be to make
% % the experiment type be a categorical variable.
% dt.Exp_nums = zeros(numsubs,1);
% Exp_names = unique(dt.Experiment);
% Exp_names = Exp_names([2,1,3]);%reorder so it's now orig, duration, timing
% for ex = 1:length(Exp_names)%lenght is 3
%   whichrows = strcmp(dt.Experiment, Exp_names{ex});%compares the exp name in the loop with the exp name in the raw data
%   dt.Exp_nums(whichrows) = ex; %boolean var whichrows makes the ones (Trues) into the loop number (1, 2, or 3).
% end%ex
% 
% %hemisphere of stimulation
% dt.LeftStim = strcmp(dt.stim_hemisphere, 'L');%1 for left; 0 for right or bilateral
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BEGIN ANALYSIS 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%OUTCOME VARIABLE: QUARTILE-BASED NOMINAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                                     Predictors with composite memory score                                                      %
% %                              You can change these two lines and the rest of the code will adapt                               %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % FYI the first row of output variables from mnrfit will be the intercept and subsequent rows will be our predictor names
% % Thus, we always have to include Intercept as the first output variable
% % We'll combine the predictors into one matrix
% 
% %Hemisphere, Experiment, IQ, Memory, Sex
%  my_predictors = [dt.LeftStim dt.Exp_nums dt.IQ_Z dt.Memory_Z dt.sex_num];
%  my_output_var_names = {'Intercept'; 'Hemisphere'; 'Experiment'; 'IQ'; 'Memory'; 'Sex'};%in order we determined above
% 
% % %Hemisphere, IQ, Memory, Sex
% %  my_predictors = [dt.LeftStim dt.IQ_Z dt.Memory_Z dt.sex_num];
% %  my_output_var_names = {'Intercept'; 'Hemisphere'; 'IQ'; 'Memory'; 'Sex'};%in order we determined above
%  
% % %Experiment, IQ, Memory, Sex
% % my_predictors = [dt.Exp_nums dt.IQ_Z dt.Memory_Z dt.sex_num];
% % my_output_var_names = {'Intercept'; 'Experiment'; 'IQ'; 'Memory'; 'Sex'};%in order we determined above
%  
% % % Hemisphere,  Memory, Sex
% % my_predictors = [dt.LeftStim dt.Memory_Z dt.sex_num];
% % my_output_var_names = {'Intercept'; 'Hemisphere'; 'Memory'; 'Sex'};%in order we determined above
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF PREDICTORS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% my_response_categories = categorical(dt.responder_status);
% %reorder to make non-responders last since last is reference category
% my_response_categories = reordercats(my_response_categories, ...
%   {'strong responder'; 'moderate responder'; 'anti-responder'; 'non-responder'});
% my_response_names = categories(my_response_categories);
% %     {'strong responder'  }  %SR
% %     {'moderate responder'}  %MR
% %     {'anti-responder'    }  %AR
% %     {'non-responder'     }  %NR
% 
% %all the stats from mnrfit will be about odds relative to reference condition, which is Non-responder here
% my_comparison_names = {'SRvNR', 'MRvNR', 'ARvNR'};
% %let's label the rows of the output stats too
% 
% 
% %use Matlab function mnrfit to do multinomial regression
% %stats is a structure that contains all the info
% %interactions are ON for this model, default for nominal models is ON
% [~,dev,stats] = mnrfit(my_predictors, my_response_categories, 'Model', 'nominal', 'Interactions', 'on');
% 
% % Calculate AIC manually
% nobs = size(my_predictors,2); % number of observations
% k = size(stats.beta,1); % number of estimated parameters
% L = -dev; % log-likelihood
% aic_value = 2*k - 2*L + (2*k*(k+1))/(nobs-k-1);
% 
% %Calculate Goodness-of-fit which compares the given model to its intercept mode
% gof = dev/stats.dfe; %uses the degrees of freedom estimate from the model
% 
% %the beta values are the "log odds of being in one category versus the reference category" per documentation.
% my_beta_table = array2table(stats.beta, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% my_beta_stderr_table = array2table(stats.se, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% my_beta_t_table = array2table(stats.t, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% my_pval_table = array2table(stats.p, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% 
% 
% %LET'S PLOT
% figure('Name', 'Log Odds of Being Strong-, Moderate-, or Anti- vs Non-Responder', 'Color', 'w', 'WindowState', 'Maximized');
% 
% %for conveneience, let's plot a suplot for each predictor
% num_predictors = size(my_predictors, 2); %take the second dimension of size (columns)
% num_comparisons = length(my_comparison_names);
% 
% clr = [223 145 167; 
%     178 198 65; 
%     134 156 211] / 255;
% 
% my_axes = NaN(1,num_predictors);%axes handles for linking Y/X axes
% 
% for sp = 1:num_predictors
%   my_axes(sp) = subplot(1,num_predictors, sp);
%   hold on
%   tmp_yvals = my_beta_table{sp+1,:};%add 1 because first row is intercept
%   tmp_errvals = my_beta_stderr_table{sp+1,:};%add 1 because first row is intercept
%   tmp_pvals =  my_pval_table{sp+1,:};%add 1 because first row is intercept
%   b = bar(1:num_comparisons, tmp_yvals, 1, 'FaceColor', 'flat');
%   b.CData = clr;
%   errorbar(1:num_comparisons, tmp_yvals, tmp_errvals,'k', 'LineStyle', 'none', 'LineWidth', 1)
%   xticks(1:num_comparisons)
%   xticklabels(my_comparison_names)
%   set(gca, 'FontName', 'arial', 'FontSize', 14, 'LineWidth', 1)
%   xlabel('Responder Comparison', 'FontWeight', 'bold', 'FontSize',17)
%   ylabel('Log Odds Relative to Non-Responder (vNR)', 'FontWeight', 'bold', 'FontSize',17)
%   title(my_output_var_names{sp+1},'FontSize', 20)%add 1 because first row is intercept
%  
% end%subplot loop
% 
% %make y axes be same across subplots
% linkaxes(my_axes, 'xy')
% 
% my_ylim = get(gca, 'YLim');
% tmp_y = my_ylim(1)+.4;%text just above bottom of figure
% 
% %let's go back through to put p values on figure above bars
% for sp = 1:num_predictors
%   subplot(my_axes(sp)), hold on
%   for b = 1:num_comparisons
%     tmpp = my_pval_table{sp+1,b};% p value for this bar
%     text(b,tmp_y, sprintf('p =\n%.3f', tmpp), 'HorizontalAlignment', 'center','FontSize', 15)
%   end%bar
% end%subplot
% %End of first figure
% 
% 
% figure
% y = dt.avg_stim_dprime_diff;
% c = dt.responder_dummy;
% x = ones(1,length(y));
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF ANALYSIS 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% 
% %make sure to comment out previous analysis sections bc we reuse variable names
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN ANALYSIS 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%OUTCOME VARIABLE: QUARTILE-BASED NOMINAL WITH RAVLT ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %                                     Predictors with RAVLT only                                                     %
% % %                      You can change these two lines and the rest of the code will adapt                               %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % % FYI the first row of output variables from mnrfit will be the intercept and subsequent rows will be our predictor names
% % % Thus, we always have to include Intercept as the first output variable
% % % We'll combine the predictors into one matrix
% % 
% % %Hemisphere, Experiment, IQ, RAVLT, Sex
% % my_predictors = [dt.LeftStim dt.Exp_nums dt.IQ_Z dt.RAVLT_Z dt.sex_num];
% % my_output_var_names = {'Intercept'; 'Hemisphere'; 'Experiment'; 'IQ'; 'RAVLT'; 'Sex'};%in order we determined above
% % 
% % %Hemisphere, IQ, RAVLT, Sex
% % % my_predictors = [dt.LeftStim dt.IQ_Z dt.RAVLT_Z dt.sex_num];
% % % my_output_var_names = {'Intercept'; 'Hemisphere'; 'IQ'; 'RAVLT'; 'Sex'};%in order we determined above
% % 
% % %Experiment, IQ, RAVLT, Sex
% % %my_predictors = [dt.Exp_nums dt.IQ_Z dt.RAVLT_Z dt.sex_num];
% % %my_output_var_names = {'Intercept'; 'Experiment'; 'IQ'; 'RAVLT'; 'Sex'};%in order we determined above
% %  
% % % Hemisphere,  RAVLT, Sex
% % % my_predictors = [dt.LeftStim dt.RAVLT_Z dt.sex_num];
% % % my_output_var_names = {'Intercept'; 'Hemisphere'; 'RAVLT'; 'Sex'};%in order we determined above
% % 
% %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF PREDICTORS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % my_response_categories = categorical(dt.responder_status);
% % %reorder to make non-responders last since last is reference category
% % my_response_categories = reordercats(my_response_categories, ...
% %   {'strong responder'; 'moderate responder'; 'anti-responder'; 'non-responder'});
% % my_response_names = categories(my_response_categories);
% % %     {'strong responder'  }  %SR
% % %     {'moderate responder'}  %MR
% % %     {'anti-responder'    }  %AR
% % %     {'non-responder'     }  %NR
% % 
% % %all the stats from mnrfit will be about odds relative to reference condition, which is Non-responder here
% % my_comparison_names = {'SRvNR', 'MRvNR', 'ARvNR'};
% % %let's label the rows of the output stats too
% % 
% % 
% % %use Matlab function mnrfit to do multinomial regression
% % %stats is a structure that contains all the info
% % %interactions are ON for this model, default for nominal models is ON
% % [~,~,dev,stats] = mnrfit(my_predictors, my_response_categories, 'Model', 'nominal', 'Interactions', 'on');
% 
% % % Calculate AIC manually
% % nobs = size(my_predictors,2); % number of observations
% % k = size(stats.beta,1); % number of estimated parameters
% % L = -dev; % log-likelihood
% % aic_value = 2*k - 2*L + (2*k*(k+1))/(nobs-k-1);
% % 
% % %Calculate Goodness-of-fit which compares the given model to its intercept mode
% % gof = dev/stats.dfe; %uses the degrees of freedom estimate from the model
% 
% % 
% % %the beta values are the "log odds of being in one category versus the reference category" per documentation.
% % my_beta_table = array2table(stats.beta, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% % my_beta_stderr_table = array2table(stats.se, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% % my_beta_t_table = array2table(stats.t, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% % my_pval_table = array2table(stats.p, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% % 
% % 
% % %LET'S PLOT
% % figure('Name', 'Log Odds of Being Strong-, Moderate-, or Anti- vs Non-Responder', 'Color', 'w', 'WindowState', 'Maximized')
% % 
% % %for conveneience, let's plot a suplot for each predictor
% % num_predictors = size(my_predictors, 2); %take the second dimension of size (columns)
% % num_comparisons = length(my_comparison_names);
% % 
% % my_axes = NaN(1,num_predictors);%axes handles for linking Y/X axes
% % 
% % for sp = 1:num_predictors
% %   my_axes(sp) = subplot(1,num_predictors, sp);
% %   hold on
% %   tmp_yvals = my_beta_table{sp+1,:};%add 1 because first row is intercept
% %   tmp_errvals = my_beta_stderr_table{sp+1,:};%add 1 because first row is intercept
% %   tmp_pvals =  my_pval_table{sp+1,:};%add 1 because first row is intercept
% %   bar(1:num_comparisons, tmp_yvals, 1, 'FaceColor', [.6 .6 .6])
% %   errorbar(1:num_comparisons, tmp_yvals, tmp_errvals,'k', 'LineStyle', 'none', 'LineWidth', 1)
% %   xticks(1:num_comparisons)
% %   xticklabels(my_comparison_names)
% %   set(gca, 'FontName', 'arial', 'FontSize', 14, 'LineWidth', 1)
% %   xlabel('Responder Comparison', 'FontWeight', 'bold', 'FontSize',17)
% %   ylabel('Log Odds Relative to Non-Responder (vNR)', 'FontWeight', 'bold', 'FontSize',17)
% %   title(my_output_var_names{sp+1},'FontSize', 20)%add 1 because first row is intercept
% %  
% % end%subplot loop
% % 
% % %make y axes be same across subplots
% % linkaxes(my_axes, 'xy')
% % 
% % my_ylim = get(gca, 'YLim');
% % tmp_y = my_ylim(1)+.4;%text just above bottom of figure
% % 
% % %let's go back through to put p values on figure above bars
% % for sp = 1:num_predictors
% %   subplot(my_axes(sp)), hold on
% %   for b = 1:num_comparisons
% %     tmpp = my_pval_table{sp+1,b};% p value for this bar
% %     text(b,tmp_y, sprintf('p =\n%.3f', tmpp), 'HorizontalAlignment', 'center','FontSize', 15)
% %   end%bar
% %  end%subplot
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF ANALYSIS 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN ANALYSIS 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%OUTCOME VARIABLE: QUARTILE-BASED ORDINAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                                      Predictors with composite memory score                                                      %
% %                               You can change these two lines and the rest of the code will adapt                               %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %FYI the first row of output variables from mnrfit will be the intercept and subsequent rows will be our predictor names
% % Thus, we always have to include Intercept as the first output variable
% %We'll combine the predictors into one matrix
% 
% %Hemisphere, Experiment, IQ, Memory, Sex
% % my_predictors = [dt.LeftStim dt.Exp_nums dt.IQ_Z dt.Memory_Z dt.sex_num];
% % my_output_var_names = {'Intercept'; 'Hemisphere'; 'Experiment'; 'IQ'; 'Memory'; 'Sex'};%in order we determined above
% 
% % %Hemisphere, IQ, Memory, Sex
% % my_predictors = [dt.LeftStim dt.IQ_Z dt.Memory_Z dt.sex_num];
% % my_output_var_names = {'Intercept'; 'Hemisphere'; 'IQ'; 'Memory'; 'Sex'};%in order we determined above
% 
% % %Experiment, IQ, Memory, Sex
% % my_predictors = [dt.Exp_nums dt.IQ_Z dt.Memory_Z dt.sex_num];
% % my_output_var_names = {'Intercept'; 'Experiment'; 'IQ'; 'Memory'; 'Sex'};%in order we determined above
%  
% % % Hemisphere,  Memory, Sex
% % my_predictors = [dt.LeftStim dt.Memory_Z dt.sex_num];
% % my_output_var_names = {'Intercept'; 'Hemisphere'; 'Memory'; 'Sex'};%in order we determined above
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF PREDICTORS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % my_response_categories = categorical(dt.responder_status, 'Ordinal', true);
% % %reorder to make non-responders last since last is reference category
% % my_response_categories = reordercats(my_response_categories, ...
% %   {'strong responder'; 'moderate responder'; 'anti-responder'; 'non-responder'});
% % my_response_names = categories(my_response_categories);
% % %     {'strong responder'  }  %SR
% % %     {'moderate responder'}  %MR
% % %     {'anti-responder'    }  %AR
% % %     {'non-responder'     }  %NR
% % 
% % %all the stats from mnrfit will be about odds relative to reference condition, which is Non-responder here
% % my_comparison_names = {'SRvNR', 'MRvNR', 'ARvNR'};
% % %let's label the rows of the output stats too
% % 
% % 
% % %use Matlab function mnrfit to do multinomial regression
% % %stats is a structure that contains all the info
% % %interactions are OFF for this model, default for ordinal
% % [~,~,dev,stats] = mnrfit(my_predictors, my_response_categories, 'Model', 'ordinal', 'Interactions', 'off');
% % 
% % 
% % %THE CODE BREAKS HERE BC THE OUTPUT IS FORMATTED DIFFERENTLY THAT I CAN'T
% % %INTERPRET THE OUTPUT AND WHAT MEANS WHAT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% % 
% % %the beta values are the "log odds of being in one category versus the reference category" per documentation.
% % my_beta_table = array2table(stats.beta, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% % my_beta_stderr_table = array2table(stats.se, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% % my_beta_t_table = array2table(stats.t, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% % my_pval_table = array2table(stats.p, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% % 
% % 
% % 
% % %LET'S PLOT
% % figure('Name', 'Log Odds of Being Strong-, Moderate-, or Anti- vs Non-Responder', 'Color', 'w', 'WindowState', 'Maximized')
% % 
% % %for conveneience, let's plot a suplot for each predictor
% % num_predictors = size(my_predictors, 2); %take the second dimension of size (columns)
% % num_comparisons = length(my_comparison_names);
% % 
% % my_axes = NaN(1,num_predictors);%axes handles for linking Y/X axes
% % 
% % for sp = 1:num_predictors
% %   my_axes(sp) = subplot(1,num_predictors, sp);
% %   hold on
% %   tmp_yvals = my_beta_table{sp+1,:};%add 1 because first row is intercept
% %   tmp_errvals = my_beta_stderr_table{sp+1,:};%add 1 because first row is intercept
% %   tmp_pvals =  my_pval_table{sp+1,:};%add 1 because first row is intercept
% %   bar(1:num_comparisons, tmp_yvals, 1, 'FaceColor', [.6 .6 .6])
% %   errorbar(1:num_comparisons, tmp_yvals, tmp_errvals,'k', 'LineStyle', 'none', 'LineWidth', 1)
% %   xticks(1:num_comparisons)
% %   xticklabels(my_comparison_names)
% %   set(gca, 'FontName', 'arial', 'FontSize', 14, 'LineWidth', 1)
% %   xlabel('Responder Comparison', 'FontWeight', 'bold', 'FontSize',17)
% %   ylabel('Log Odds Relative to Non-Responder (vNR)', 'FontWeight', 'bold', 'FontSize',17)
% %   title(my_output_var_names{sp+1},'FontSize', 20)%add 1 because first row is intercept
% %  
% % end%subplot loop
% % 
% % %make y axes be same across subplots
% % linkaxes(my_axes, 'xy')
% % 
% % my_ylim = get(gca, 'YLim');
% % tmp_y = my_ylim(1)+.4;%text just above bottom of figure
% % 
% % %let's go back through to put p values on figure above bars
% % for sp = 1:num_predictors
% %   subplot(my_axes(sp)), hold on
% %   for b = 1:num_comparisons
% %     tmpp = my_pval_table{sp+1,b};% p value for this bar
% %     text(b,tmp_y, sprintf('p =\n%.3f', tmpp), 'HorizontalAlignment', 'center','FontSize', 15)
% %   end%bar
% % end%subplot
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF ANALYSIS 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN ANALYSIS 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%OUTCOME VARIABLE: QUARTILE-BASED ORDINAL WITH RAVLT ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                                     Predictors with RAVLT only                                                     %
% %                      You can change these two lines and the rest of the code will adapt                               %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % FYI the first row of output variables from mnrfit will be the intercept and subsequent rows will be our predictor names
% % Thus, we always have to include Intercept as the first output variable
% % We'll combine the predictors into one matrix
% 
% %Hemisphere, Experiment, IQ, RAVLT, Sex
% % my_predictors = [dt.LeftStim dt.Exp_nums dt.IQ_Z dt.RAVLT_Z dt.sex_num];
% % my_output_var_names = {'Intercept'; 'Hemisphere'; 'Experiment'; 'IQ'; 'RAVLT'; 'Sex'};%in order we determined above
% 
% %Hemisphere, IQ, RAVLT, Sex
% % my_predictors = [dt.LeftStim dt.IQ_Z dt.RAVLT_Z dt.sex_num];
% % my_output_var_names = {'Intercept'; 'Hemisphere'; 'IQ'; 'RAVLT'; 'Sex'};%in order we determined above
% 
% %Experiment, IQ, RAVLT, Sex
% %my_predictors = [dt.Exp_nums dt.IQ_Z dt.RAVLT_Z dt.sex_num];
% %my_output_var_names = {'Intercept'; 'Experiment'; 'IQ'; 'RAVLT'; 'Sex'};%in order we determined above
%  
% % Hemisphere,  RAVLT, Sex
% % my_predictors = [dt.LeftStim dt.RAVLT_Z dt.sex_num];
% % my_output_var_names = {'Intercept'; 'Hemisphere'; 'RAVLT'; 'Sex'};%in order we determined above
% 
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF PREDICTORS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % my_response_categories = categorical(dt.responder_status, 'Ordinal', true);
% % %reorder to make non-responders last since last is reference category
% % my_response_categories = reordercats(my_response_categories, ...
% %   {'strong responder'; 'moderate responder'; 'anti-responder'; 'non-responder'});
% % my_response_names = categories(my_response_categories);
% % %     {'strong responder'  }  %SR
% % %     {'moderate responder'}  %MR
% % %     {'anti-responder'    }  %AR
% % %     {'non-responder'     }  %NR
% % 
% % %all the stats from mnrfit will be about odds relative to reference condition, which is Non-responder here
% % my_comparison_names = {'SRvNR', 'MRvNR', 'ARvNR'};
% % %let's label the rows of the output stats too
% % 
% % 
% % %use Matlab function mnrfit to do multinomial regression
% % %stats is a structure that contains all the info
% % %interactions are OFF for this model, default for ordinal
% % [~,~,dev,stats] = mnrfit(my_predictors, my_response_categories, 'Model', 'ordinal', 'Interactions', 'off');
% % 
% % 
% % %THE CODE BREAKS HERE BC THE OUTPUT IS FORMATTED DIFFERENTLY THAT I CAN'T
% % %INTERPRET THE OUTPUT AND WHAT MEANS WHAT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% % 
% % %the beta values are the "log odds of being in one category versus the reference category" per documentation.
% % my_beta_table = array2table(stats.beta, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% % my_beta_stderr_table = array2table(stats.se, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% % my_beta_t_table = array2table(stats.t, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% % my_pval_table = array2table(stats.p, 'VariableNames', my_comparison_names, 'RowNames', my_output_var_names);
% % 
% % %LET'S PLOT
% % figure('Name', 'Log Odds of Being Strong-, Moderate-, or Anti- vs Non-Responder', 'Color', 'w', 'WindowState', 'Maximized')
% % 
% % %for conveneience, let's plot a suplot for each predictor
% % num_predictors = size(my_predictors, 2); %take the second dimension of size (columns)
% % num_comparisons = length(my_comparison_names);
% % 
% % my_axes = NaN(1,num_predictors);%axes handles for linking Y/X axes
% % 
% % for sp = 1:num_predictors
% %   my_axes(sp) = subplot(1,num_predictors, sp);
% %   hold on
% %   tmp_yvals = my_beta_table{sp+1,:};%add 1 because first row is intercept
% %   tmp_errvals = my_beta_stderr_table{sp+1,:};%add 1 because first row is intercept
% %   tmp_pvals =  my_pval_table{sp+1,:};%add 1 because first row is intercept
% %   bar(1:num_comparisons, tmp_yvals, 1, 'FaceColor', [.6 .6 .6])
% %   errorbar(1:num_comparisons, tmp_yvals, tmp_errvals,'k', 'LineStyle', 'none', 'LineWidth', 1)
% %   xticks(1:num_comparisons)
% %   xticklabels(my_comparison_names)
% %   set(gca, 'FontName', 'arial', 'FontSize', 14, 'LineWidth', 1)
% %   xlabel('Responder Comparison', 'FontWeight', 'bold', 'FontSize',17)
% %   ylabel('Log Odds Relative to Non-Responder (vNR)', 'FontWeight', 'bold', 'FontSize',17)
% %   title(my_output_var_names{sp+1},'FontSize', 20)%add 1 because first row is intercept
% %  
% % end%subplot loop
% % 
% % %make y axes be same across subplots
% % linkaxes(my_axes, 'xy')
% % 
% % my_ylim = get(gca, 'YLim');
% % tmp_y = my_ylim(1)+.4;%text just above bottom of figure
% % 
% % %let's go back through to put p values on figure above bars
% % for sp = 1:num_predictors
% %   subplot(my_axes(sp)), hold on
% %   for b = 1:num_comparisons
% %     tmpp = my_pval_table{sp+1,b};% p value for this bar
% %     text(b,tmp_y, sprintf('p =\n%.3f', tmpp), 'HorizontalAlignment', 'center','FontSize', 15)
% %   end%bar
% % end%subplot
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF ANALYSIS 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
