%Script to do multiple regression analyses on extended AMME dataset
%JRM 09/13/23




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    READ IN DATA AS A TABLE AND DO A BIT OF PREP WORK ON VARIABLES PRIOR TO REGRESSION    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;

%Read in csv file as a data table
dt = readtable('foranalysis_responder_firstsession_quartilebased.csv');
numsubs = height(dt);%number of subjects

%Let's Z-score the two memory scores and average them for an overall memory score
dt.RAVLT_Z = normalize(dt.RAVLT_dprime);
dt.ReyO_Z = normalize(dt.ReyO_CF_delayedrecall);
dt.Memory_Z = nanmean([dt.RAVLT_Z dt.ReyO_Z], 2);%nanmean because two patients are missing one memory score

%Let's also Z score the IQ scores so that the Beta scaling better aligns with other variable
dt.IQ_Z = normalize(dt.FSIQ);

%make sex be a categorical value
dt.sex = categorical(dt.sex);

%make Experiment names be numbers too
%A better approach, if it's possible, would be to make
% the experiment type be a categorical variable.
dt.Exp_nums = zeros(numsubs,1);
Exp_names = unique(dt.Experiment);
Exp_names = Exp_names([2,1,3]);%reorder so it's now orig, duration, timing
for ex = 1:length(Exp_names)
  whichrows = strcmp(dt.Experiment, Exp_names{ex});
  dt.Exp_nums(whichrows) = ex;
end%ex

%hemisphere of stimulation
dt.LeftStim = strcmp(dt.stim_hemisphere, 'L');%1 for left; 0 for right or bilateral

%IED activity 
dt.IEDday1 = strcmp(dt.IED_duringimageday1, 'yes'); %1 for yes; 0 fr no
dt.IEDday2 = strcmp(dt.IED_duringimageday2, 'yes'); %1 for yes; 0 fr no
dt.IEDbilateral = strcmp(dt.IED_bilateral, 'yes'); %1 for yes; 0 fr no

%create a new variable that is just the absolute value of the stim-nostim dprime difference
dt.abs_avg_stim_dprime_diff = abs(dt.avg_stim_dprime_diff);


% % Run correlation between distance and dprime at 1-day test
 [r, p] = corrcoef(dt.IED_freq, dt.abs_avg_stim_dprime_diff,'rows','complete');
 scatter(dt.IED_freq, dt.abs_avg_stim_dprime_diff, 'filled')
 xlabel('IED freq')
 ylabel('Abs value dprime at 1-day test')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 DO MULTIPLE REGRESSION                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We will use Matlab's fitlm function, which uses Wilkinson notatation, similar to many R packages
% The basic version is 'y~x1+x2+x3' for an intercept (default), three independent predictors (x's), 
%  and one dependent outcome (y), with no interaction terms
% To include all two and three way interaction terms, use 'y~x1*x2*x3'



fprintf('Fitting multiple regression model with signed d-prime differences\n')
modelnum = 1;%first one; signed y values
mymodels(modelnum).description = 'Signed d'' stim-nostim differences';
mymodels(modelnum).lm = fitlm(dt, 'avg_stim_dprime_diff ~ IQ_Z+ Memory_Z + sex + age + shortest_dist_MTL + IED_freq');
%print out some basic stats
mymodels(modelnum).lm.disp
%print out "ANOVA" style table at command line
mymodels(modelnum).lm.anova;
%-->Baseline memory is not significant in this model

%calculate AIC
loglikelihood = mymodels(modelnum).lm.LogLikelihood;
k= mymodels(modelnum).lm.NumPredictors;
AIC1 = 2*k - 2*loglikelihood;


fprintf('\n\nFitting multiple regression model with absolute value d-prime differences\n')
modelnum = 2;%second one; absolute y values
mymodels(modelnum).description = 'Absolute d'' stim-nostim differences';
mymodels(modelnum).lm = fitlm(dt, 'abs_avg_stim_dprime_diff ~ IQ_Z + Memory_Z + sex + age + shortest_dist_MTL+ IED_freq');
%print out some basic stats
mymodels(modelnum).lm.disp

%print out "ANOVA" style table at command line
mymodels(modelnum).lm.anova

%-->Baseline memory now is significant predictor of absolute dprime differences. 

% %calculate AIC
% loglikelihood = mymodels(modelnum).lm.LogLikelihood;
% k= mymodels(modelnum).lm.NumPredictors;
% AIC2 = 2*k - 2*loglikelihood;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                          PLOT                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure('Name', 'Mutiple Regression With Extended AMME dataset; signed and absolute dprime differences', 'Color', 'w', 'WindowState', 'Maximized')

%for conveneience, let's plot a suplot for each predictor
%let's plot a separate row for each model; top = signed, bottom = absolute
num_predictors = mymodels(1).lm.NumPredictors;%assume same for all models
num_models = length(mymodels);%one for signed data; one for absolute values of stim-nostim dprime differences

my_axes = NaN(num_models,num_predictors);%axes handles for linking Y/X axes

for modelnum = 1:num_models
  for pnum = 1:num_predictors
    my_axes(modelnum, pnum) = subplot(num_models,num_predictors,6*(modelnum-1)+pnum);
    %the LinearModel object has a built-in plot adjusted response function
    tmpfighandle = mymodels(modelnum).lm.plotAdjustedResponse(mymodels(modelnum).lm.PredictorNames{pnum});
    ylabel('Adjusted d'' Differences', 'FontName', 'arial', 'FontWeight', 'bold', 'FontSize', 10)
    xlabel(mymodels(modelnum).lm.PredictorNames{pnum}, 'FontName', 'arial', 'FontWeight', 'bold', 'FontSize', 10, 'interpreter', 'none')
    %set(gca, 'FontName', 'arial', 'FontSize', 10, 'LineWidth', 1)
    %title(mymodels(modelnum).description, 'FontName', 'arial', 'FontWeight', 'bold', 'FontSize', 12)
    legend off
    %make the fitted line black
    tmpfighandle(2).Color = 'k';
    %make underlying lmplot marker none to overwrite it with a gscatter
    tmpfighandle(1).Color = 'none';
    %set(tmpfighandle(1), 'Marker', 'x'); % if we wanted to change the 

    hold on
    %use group scatter and overlay it to recolor circles red for women, blue for men
    p = gscatter(tmpfighandle(1).XData,tmpfighandle(1).YData,dt.sex,'rb', 'o', 10, 'off'); %rb is the order of red,blue bc it accounted females first
    %set the legend outside the plot
    %p = findobj(gcf,'tag','legend'); 
    %set(p,'location','northeastoutside');
    %p(1).LineWidth = 2; %change the linewidt of the points - looks kinda weird
    %p(2).LineWidth = 2;

  end%predictors
  linkaxes(my_axes(modelnum,:), 'y')%make y axes be same across subplots of same model
  
  tmpanova = mymodels(modelnum).lm.anova;
  for pnum = 1:num_predictors
    %print p value on plot
    tmppval = tmpanova{pnum,end};
    subplot(num_models,num_predictors,6*(modelnum-1)+pnum);
    hold on
    myylim = ylim;
    myxlim = xlim;
    text(myxlim(2)*.95, myylim(2)*.95, sprintf('p=%.4f', tmppval), 'Color', 'black',  ...
      'FontName', 'arial', 'FontWeight', 'bold', 'FontSize', 12, ...
      'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
  end
  
  
end%models






