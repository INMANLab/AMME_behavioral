
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats
from scipy.stats import f_oneway, norm, uniform, kruskal, mannwhitneyu
import researchpy as rp
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.graphics.gofplots import qqplot
import scikit_posthocs as sp



#Define the dataframe and make a boxplot.
df = pd.read_csv('foranalysis_responder_firstsession_quartilebased.csv')
df = df.dropna()
print (df)


#create a mask 
strong_responder = df['responder_dummy'] == 1
moderate_responder = df['responder_dummy'] == 0.5
nonresponder = df['responder_dummy'] == 0
antiresponder = df['responder_dummy'] == -1


################### FSIQ RESULTS ###################################

#Define each group for ANOVA
#create a mask 
strong_responder = df['responder_dummy'] == 1
moderate_responder = df['responder_dummy'] ==0.5
nonresponder = df['responder_dummy'] == 0
antiresponder = df['responder_dummy'] == -1

#filter df by the mask
sresp_FSIQ = df[strong_responder]
mresp_FSIQ = df[moderate_responder]
nresp_FSIQ = df[nonresponder]
anresp_FSIQ  = df[antiresponder]
#print(resp_FSIQ)

#create ANOVA variables
sresp_FSIQ = sresp_FSIQ['FSIQ']
mresp_FSIQ = mresp_FSIQ['FSIQ']
nresp_FSIQ = nresp_FSIQ['FSIQ']
anresp_FSIQ = anresp_FSIQ['FSIQ']

print(df.columns)
print(df['responder_status'])

#print(resp_FSIQ)
#print(nresp_FSIQ)
#print(anresp_FSIQ)

#checking normality: if p < 0.05 then not normally distributed data
#shapiro = stats.shapiro(df['FSIQ'])
#print(shapiro)

#checking homogeneity of variance: if p < 0.05 then the sample is heterogeneous
#levene = stats.levene(resp_FSIQ, anresp_FSIQ)
#print(levene)

#checking distribution normality 
#fig, axs = plt.subplots(1, 1, figsize =(10, 7), tight_layout = True)
#axs.hist(df['FSIQ'])
#plt.xlabel("FSIQ")
#plt.show()

#checking normality assumption with QQ plot: if not close to the line it's not normally distributed
#plt.rcParams['figure.figsize'] = [10, 7]
#plt.rc('font', size=14)
#qqplot(df['FSIQ'],norm,fit=True,line="45")
#plt.show()



#perform ANOVA
fvalue, pvalue = stats.f_oneway(sresp_FSIQ, mresp_FSIQ, nresp_FSIQ, anresp_FSIQ)
print('FSIQ ANOVA')
print(fvalue, pvalue)

#run post-hoc test
tukey = pairwise_tukeyhsd(endog=df['FSIQ'],
                          groups=df['responder_status'],
                          alpha=0.05)

print(tukey)

#independent t-test
summary, results= rp.ttest(sresp_FSIQ, anresp_FSIQ)
print('FSIQ t-test strong resp vs antiresp') 
print(summary, results)

summary, results= rp.ttest(sresp_FSIQ, nresp_FSIQ)
print('FSIQ t-test strong resp vs nonresp') 
print(summary, results)

summary, results= rp.ttest(nresp_FSIQ, anresp_FSIQ)
print('FSIQ t-test nonresp vs antiresp') 
print(summary, results)

summary, results= rp.ttest(sresp_FSIQ, mresp_FSIQ)
print('FSIQ t-test strong resp vs moderate resp') 
print(summary, results)

cor= stats.pearsonr(df['FSIQ'], df['avg_stim_dprime_diff'])
print("corr=", cor)


sns.barplot(data = df, x = 'responder_status', y = 'FSIQ', palette='Blues_r', errorbar=('ci', 95), order=['strong responder', 'moderate responder', 'non-responder', 'anti-responder'])
sns.stripplot(data = df, x = 'responder_status', y = 'FSIQ', size=10, color='black', alpha=.7, order=['strong responder', 'moderate responder', 'non-responder', 'anti-responder'])
#plt.title('FSIQ by Responder Status')
plt.ylim([df['FSIQ'].min()-5, df['FSIQ'].max()+5])
plt.show()

#scatterplot or linear model plot 
#sns.lmplot(x="FSIQ", y="avg_stim_dprime_diff", data=df)
#plt.xlim([df['FSIQ'].min()-1, df['FSIQ'].max()+1])
#plt.show()

sns.scatterplot(x="FSIQ", y="avg_stim_dprime_diff", data=df)
plt.ylabel('Stimulation dprime difference')
plt.xlabel('FSIQ')
plt.show()


################### RAVLT RESULTS ###################################


#filter df by the mask
sresp_RAVLT = df[strong_responder]
mresp_RAVLT = df[moderate_responder]
nresp_RAVLT = df[nonresponder]
anresp_RAVLT  = df[antiresponder]

#create variables
sresp_RAVLT = sresp_RAVLT['RAVLT_dprime']
mresp_RAVLT = mresp_RAVLT['RAVLT_dprime']
nresp_RAVLT= nresp_RAVLT['RAVLT_dprime']
anresp_RAVLT = anresp_RAVLT['RAVLT_dprime']

print(sresp_RAVLT)
print(mresp_RAVLT)
print(nresp_RAVLT)
print(anresp_RAVLT)

#perform ANOVA
fvalue, pvalue = stats.f_oneway(sresp_RAVLT, mresp_RAVLT,nresp_RAVLT, anresp_RAVLT)
print('RAVLT ANOVA')
print(fvalue, pvalue)

#run post-hoc test
tukey = pairwise_tukeyhsd(endog=df['RAVLT_dprime'],
                          groups=df['responder_status'],
                          alpha=0.05)

print(tukey)

#independent t-test
summary, results= rp.ttest(sresp_RAVLT, anresp_RAVLT)
print('RAVLT t-test strong resp vs antiresp') 
print(summary, results)

summary, results= rp.ttest(sresp_RAVLT, mresp_RAVLT)
print('RAVLT t-test strong resp vs moderate resp') 
print(summary, results)

summary, results= rp.ttest(sresp_RAVLT, nresp_RAVLT)
print('RAVLT t-test sresp vs nonresp') 
print(summary, results)

summary, results= rp.ttest(nresp_RAVLT, anresp_RAVLT)
print('RAVLT t-test nonresp vs antiresp') 
print(summary, results)

cor= stats.pearsonr(df['RAVLT_dprime'], df['avg_stim_dprime_diff'])
print("corr=", cor)

sns.barplot(data = df, x = 'responder_status', y = 'RAVLT_dprime', palette='Blues_r', errorbar=('ci', 95),order=['strong responder', 'moderate responder', 'non-responder', 'anti-responder'])
sns.swarmplot(data = df, x = 'responder_status', y = 'RAVLT_dprime', size=10, color='black', alpha=.7,order=['strong responder', 'moderate responder', 'non-responder', 'anti-responder'])
#plt.axhline(15, ls='-', color='k')
#plt.title('RAVLT Dprime by Responder Status')
plt.ylim([df['RAVLT_dprime'].min(), df['RAVLT_dprime'].max()+5])
plt.show()

#scatterplot or linear model plot 
#sns.lmplot(x="RAVLT_dprime", y="avg_stim_dprime_diff", data=df)
#plt.xlim([df['RAVLT_dprime'].min()-.2, df['RAVLT_dprime'].max()+.2])
#plt.show()

#simple scatter plot
sns.scatterplot(x="RAVLT_dprime", y="avg_stim_dprime_diff", data=df)
plt.ylabel('Stimulation dprime difference')
plt.xlabel('RAVLT dprime')
plt.show()



################### REY-O COMPLEX FIGURE TASK DELAYED RECALL RESULTS ###################################


#filter df by the mask
sresp_reyo = df[strong_responder]
mresp_reyo = df[moderate_responder]
nresp_reyo = df[nonresponder]
anresp_reyo = df[antiresponder]

#create ANOVA variables
sresp_reyo = sresp_reyo['ReyO_CF_delayedrecall']
mresp_reyo = mresp_reyo['ReyO_CF_delayedrecall']
nresp_reyo = nresp_reyo['ReyO_CF_delayedrecall']
anresp_reyo = anresp_reyo['ReyO_CF_delayedrecall']

print(sresp_reyo)
print(mresp_reyo)
print(nresp_reyo)
print(anresp_reyo)


#perform ANOVA
fvalue, pvalue = stats.f_oneway(sresp_reyo, mresp_reyo, nresp_reyo, anresp_reyo)
print('Rey-O ANOVA')
print(fvalue, pvalue)

#run post-hoc test
tukey = pairwise_tukeyhsd(endog=df['ReyO_CF_delayedrecall'],
                          groups=df['responder_status'],
                          alpha=0.05)

print(tukey)

#independent t-test
summary, results= rp.ttest(sresp_reyo, anresp_reyo)
print('Rey-o t-test strong resp vs antiresp') 
print(summary, results)

summary, results= rp.ttest(sresp_reyo, mresp_reyo)
print('Rey-o t-test strong resp vs moderate resp') 
print(summary, results)

summary, results= rp.ttest(sresp_reyo, nresp_reyo)
print('Rey-o t-test strong resp vs nonresp') 
print(summary, results)

summary, results= rp.ttest(nresp_reyo, anresp_reyo)
print('Rey-o t-test nonresp vs antiresp') 
print(summary, results)

summary, results= rp.ttest(mresp_reyo, anresp_reyo)
print('Rey-o t-test moderate resp vs antiresp') 
print(summary, results)


sns.barplot(data = df, x = 'responder_status', y = 'ReyO_CF_delayedrecall', palette='Blues_r', errorbar=('ci', 95),order=['strong responder', 'moderate responder', 'non-responder', 'anti-responder'])
sns.stripplot(data = df, x = 'responder_status', y = 'ReyO_CF_delayedrecall', size=10, color='black', alpha=.7,order=['strong responder', 'moderate responder', 'non-responder', 'anti-responder'])
#plt.title('Rey-O CF delayed recall by Responder Status')
plt.ylabel('Percent Correct Rey-O CF delayed recall')
plt.ylim([df['ReyO_CF_delayedrecall'].min(), df['ReyO_CF_delayedrecall'].max()+.1])
plt.show()

#scatterplot or linear model plot 
#sns.lmplot(x="ReyO_CF_delayedrecall", y="avg_stim_dprime_diff", data=df)
#plt.xlim([df['ReyO_CF_delayedrecall'].min()-.01, df['ReyO_CF_delayedrecall'].max()+.01])
#plt.show()

#simple scatter plot
sns.scatterplot(x="ReyO_CF_delayedrecall", y="avg_stim_dprime_diff", data=df)
plt.ylabel('Stimulation dprime difference')
plt.xlabel('Rey-O CF delayed recall')
plt.show()


cor= stats.pearsonr(df['ReyO_CF_delayedrecall'], df['avg_stim_dprime_diff'])
print("pearson corr=", cor)



################### COMBINED BASELINE MEMORY SCORE (Z-SCORED) ###################################


#filter df by the mask
sresp_mem = df[strong_responder]
mresp_mem = df[moderate_responder]
nresp_mem = df[nonresponder]
anresp_mem = df[antiresponder]

#create ANOVA variables
sresp_mem = sresp_mem['Memory']
mresp_mem = mresp_mem['Memory']
nresp_mem = nresp_mem['Memory']
anresp_mem = anresp_mem['Memory']

print(sresp_mem)
print(mresp_mem)
print(nresp_mem)
print(anresp_mem)


#perform ANOVA
fvalue, pvalue = stats.f_oneway(sresp_mem, mresp_mem, nresp_mem, anresp_mem)
print('Memory ANOVA')
print(fvalue, pvalue)

#run post-hoc test
tukey = pairwise_tukeyhsd(endog=df['Memory'],
                          groups=df['responder_status'],
                          alpha=0.05)

print(tukey)

#independent t-test
summary, results= rp.ttest(sresp_mem, anresp_mem)
print('Memory t-test strong resp vs antiresp') 
print(summary, results)

summary, results= rp.ttest(sresp_mem, mresp_mem)
print('Memory t-test strong resp vs moderate resp') 
print(summary, results)

summary, results= rp.ttest(sresp_mem, nresp_mem)
print('Memory t-test strong resp vs nonresp') 
print(summary, results)

summary, results= rp.ttest(nresp_mem, anresp_mem)
print('Memory t-test nonresp vs antiresp') 
print(summary, results)

summary, results= rp.ttest(mresp_mem, anresp_mem)
print('Memory t-test moderate resp vs antiresp') 
print(summary, results)


sns.barplot(data = df, x = 'responder_status', y = 'Memory', palette='Blues_r', errorbar=('ci', 95),order=['strong responder', 'moderate responder', 'non-responder', 'anti-responder'])
sns.stripplot(data = df, x = 'responder_status', y = 'Memory', size=10, color='black', alpha=.7,order=['strong responder', 'moderate responder', 'non-responder', 'anti-responder'])
#plt.title('Memory by Responder Status')
plt.ylabel('Memory Composite score (RAVLT + Rey-O)')
plt.ylim([df['Memory'].min()+.5, df['Memory'].max()+.5])
plt.show()

#scatterplot or linear model plot 
#sns.lmplot(x="Memory", y="avg_stim_dprime_diff", data=df)
#plt.xlim([df['Memory'].min()-.1, df['Memory'].max()+.1])
#plt.show()

cor= stats.pearsonr(df['Memory'], df['avg_stim_dprime_diff'])
print("pearson corr=", cor)

#simple scatter plot
sns.scatterplot(x="Memory", y="avg_stim_dprime_diff", data=df)
plt.ylabel('Stimulation dprime difference')
plt.xlabel('Baseline Memory')
plt.show()


################### AGE RESULTS ###################################


#filter df by the mask
sresp_age = df[strong_responder]
mresp_age = df[moderate_responder]
nresp_age = df[nonresponder]
anresp_age = df[antiresponder]

#create ANOVA variables
sresp_age = sresp_age['age']
mresp_age = mresp_age['age']
nresp_age= nresp_age['age']
anresp_age = anresp_age['age']

print(sresp_age)
print(mresp_age)
print(nresp_age)
print(anresp_age)


#perform ANOVA
fvalue, pvalue = stats.f_oneway(sresp_age, mresp_age, nresp_age, anresp_age)
print('age ANOVA')
print(fvalue, pvalue)

#run post-hoc test
tukey = pairwise_tukeyhsd(endog=df['age'],
                          groups=df['responder_status'],
                          alpha=0.05)

print(tukey)

#independent t-test
summary, results= rp.ttest(sresp_age, anresp_age)
print('age t-test strong resp vs antiresp') 
print(summary, results)

summary, results= rp.ttest(sresp_age, mresp_age)
print('age t-test strong resp vs moderate resp') 
print(summary, results)

summary, results= rp.ttest(sresp_age, nresp_age)
print('age t-test strong resp vs nonresp') 
print(summary, results)

summary, results= rp.ttest(nresp_age, anresp_age)
print('age t-test nonresp vs antiresp') 
print(summary, results)

summary, results= rp.ttest(nresp_age, mresp_age)
print('age t-test nonresp vs moderate resp') 
print(summary, results)

cor= stats.pearsonr(df['age'], df['avg_stim_dprime_diff'])
print("pearson corr=", cor)

sns.barplot(data = df, x = 'responder_status', y = 'age', palette='Blues_r', errorbar=('ci', 95),order=['strong responder', 'moderate responder', 'non-responder', 'anti-responder'])
sns.stripplot(data = df, x = 'responder_status', y = 'age', size=10, color='black', alpha=.7,order=['strong responder', 'moderate responder', 'non-responder', 'anti-responder'])
#plt.title('Age by Responder Status')
plt.ylim([df['age'].min()-3, df['age'].max()+20])
plt.show()

#scatterplot or linear model plot 
#sns.lmplot(x="age", y="avg_stim_dprime_diff", data=df)
#plt.xlim([df['age'].min()-1, df['age'].max()+1])
#plt.show()

#simple scatter plot
sns.scatterplot(x="age", y="avg_stim_dprime_diff", data=df)
plt.ylabel('Stimulation dprime difference')
plt.xlabel('Age')
plt.show()


####### DISTANCE FROM BLA STIM ELECTRODES TO CLOSEST HPC VOXEL ############

#filter df by the mask
sresp_dist = df[strong_responder]
mresp_dist = df[moderate_responder]
nresp_dist = df[nonresponder]
anresp_dist = df[antiresponder]

#create ANOVA variables
sresp_dist = sresp_dist['dist_BLAtoHPC']
mresp_dist = mresp_dist['dist_BLAtoHPC']
nresp_dist = nresp_dist['dist_BLAtoHPC']
anresp_dist = anresp_dist['dist_BLAtoHPC']

print(sresp_dist)
print(mresp_dist)
print(nresp_dist)
print(anresp_dist)


#perform ANOVA
fvalue, pvalue = stats.f_oneway(sresp_dist, mresp_dist, nresp_dist, anresp_dist)
print('Dist BLA to HPC ANOVA')
print(fvalue, pvalue)

#run post-hoc test
tukey = pairwise_tukeyhsd(endog=df['dist_BLAtoHPC'],
                          groups=df['responder_status'],
                          alpha=0.05)

print(tukey)

#independent t-test
summary, results= rp.ttest(sresp_dist, anresp_dist)
print('dist t-test strong resp vs antiresp') 
print(summary, results)

summary, results= rp.ttest(sresp_dist, mresp_dist)
print('dist t-test strong resp vs moderate resp') 
print(summary, results)

summary, results= rp.ttest(sresp_dist, nresp_dist)
print('dist t-test strong resp vs nonresp') 
print(summary, results)

summary, results= rp.ttest(nresp_dist, anresp_dist)
print('dist t-test nonresp vs antiresp') 
print(summary, results)

summary, results= rp.ttest(nresp_dist, mresp_dist)
print('dist t-test nonresp vs moderate resp') 
print(summary, results)

cor= stats.pearsonr(df['dist_BLAtoHPC'], df['avg_stim_dprime_diff'])
print("pearson corr=", cor)

sns.barplot(data = df, x ='responder_status', y = 'dist_BLAtoHPC', palette='Blues_r', errorbar=('ci', 95),order=['strong responder', 'moderate responder', 'non-responder', 'anti-responder'])
sns.stripplot(data = df, x = 'responder_status', y = 'dist_BLAtoHPC', size=10, color='black', alpha=.7,order=['strong responder', 'moderate responder', 'non-responder', 'anti-responder'])
#plt.title('Distance from BLA to HPC by Responder')
#plt.ylim([df['age'].min(), df['age'].max()])
plt.show()

#scatterplot or linear model plot 
#sns.lmplot(x="dist_BLAtoHPC", y="avg_stim_dprime_diff", data=df)
#plt.xlim([df['dist_BLAtoHPC'].min()-.1, df['dist_BLAtoHPC'].max()+.1])
#plt.show()

#simple scatter plot
sns.scatterplot(x="dist_BLAtoHPC", y="avg_stim_dprime_diff", data=df)
plt.ylabel('Stimulation dprime difference')
plt.xlabel('Distance to HPC')
plt.show()


##################### END #################################################



'''
####### non-parametric tests ###########


#non-parametric anova - Kruskal
#fvalue, pvalue = stats.f_oneway(resp_RAVLT, nresp_RAVLT, anresp_RAVLT)
statistic, pvalue = stats.kruskal(resp_RAVLT, nresp_RAVLT, anresp_RAVLT)
print('ANOVA')
print(statistic, pvalue)

#perform non-parametric post hoc test
summary= sp.posthoc_dunn([resp_RAVLT, nresp_RAVLT, anresp_RAVLT], p_adjust = 'sidak')
print(summary)


#independent non-parametric test
statistic, pvalue = stats.mannwhitneyu(resp_RAVLT, anresp_RAVLT)
print('Mann-Whitney t-test resp vs antiresp') 
print(statistic, pvalue)

statistic, pvalue = stats.mannwhitneyu(resp_RAVLT, nresp_RAVLT)
print('Mann-Whitney t-test resp vs nonresp') 
print(statistic, pvalue)

statistic, pvalue = stats.mannwhitneyu(nresp_RAVLT, anresp_RAVLT)
print('Mann-Whitney t-test nonresp vs antiresp') 
print(statistic, pvalue)
'''
