#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np

#read csv file into pandas dataframe
df1 = pd.read_csv('sequencing_data_biochem1.csv') 

#replace 0's with NaN so that cycles with no reads in all 4 columns will return NaN by idxmax() function
df1.replace(to_replace={0.000: np.nan}, inplace=True)

#select columns for second cycle of biochem1 and rename columns to match format of call_1
second_cycle = df1[['A_2','C_2','G_2','T_2']].rename(columns={'A_2':'A', 'C_2':'C', 'G_2':'G', 'T_2':'T'}) 
#return column name of max value of each cycle, NaN if 0
df1['call_2']=second_cycle.idxmax(axis=1) 


#replace NaN in base call column with N to indicate no reads
df1['call_2'].replace(to_replace={np.nan: 'N'}, inplace=True) 
#convert NaN back to 0.000 in rest of dataframe
df1.fillna(0, inplace=True) 
#with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    #print(df1)



#compare ref_1 to call_1 and record 1 for every mismatch or 0 if the nt readings match
df1.insert(6, 'mismatch_1', np.where(df1['ref_1'] != df1['call_1'], 1, 0)) 
#compare ref_2 to call_2 and record 1 for every mismatch or 0 if the nt readings match
df1['mismatch_2'] = np.where(df1['ref_2'] != df1['call_2'], 1, 0) 



#calculate error for both cycles and total error for biochem1
#done by summing 1's in mismatch column and dividing by total reads and multiplying by 100

biochem1_error1 = (sum(df1['mismatch_1']))/len(df1['ref_1'])*100 
biochem1_error2 = (sum(df1['mismatch_2']))/len(df1['ref_1'])*100
biochem1_totalerror = (sum(df1['mismatch_1']) + sum(df1['mismatch_2']))/(2*len(df1['ref_1']))*100

print(biochem1_error1, biochem1_error2, biochem1_totalerror)



#Biochem 2

df2 = pd.read_csv('sequencing_data_biochem2.csv')

df2.replace(to_replace = {0.000: np.nan}, inplace=True)
second_cycle1 = df2[['A_2','C_2','G_2','T_2']].rename(columns={'A_2':'A', 'C_2':'C', 'G_2':'G', 'T_2':'T'})
df2['call_2'] = second_cycle1.idxmax(axis=1)

#replace NaN with N
df2['call_2'].replace(to_replace = {np.nan: 'N'}, inplace=True)
#convert NaN back to 0.000
df2.fillna(0, inplace=True) 

#create dictionary to compensate for software error
mydict = {'C':'G', 'T':'A', 'A': 'C', 'G':'T'}

#replace values in both call columns with correct reading
df2['call_1'] = df2['call_1'].replace(mydict)
df2['call_2'] = df2['call_2'].replace(mydict)
#df2



df2.insert(6, 'mismatch_1', np.where(df2['ref_1'] != df2['call_1'], 1, 0))
df2['mismatch_2'] = np.where(df2['ref_2'] != df2['call_2'], 1, 0)



#calculate error for both cycles and total error for biochem1
#done by summing 1's in mismatch column and dividing by total reads and multiplying by 100

biochem_2_error_1 = (sum(df2['mismatch_1']))/len(df2['ref_1'])*100
biochem_2_error_2 = (sum(df2['mismatch_2']))/len(df2['ref_1'])*100
biochem_2_totalerror = (sum(df2['mismatch_1']) + sum(df2['mismatch_2']))/(2*len(df2['ref_1']))*100

print(biochem_2_error_1, biochem_2_error_2, biochem_2_totalerror)



##Further Analysis

####Below code tests to see what percentage of mismatches were result of no reads (N)

percent_nans1 = (sum(np.where(df1['call_1'] == 'N', 1, 0)))/len(df1['ref_1'])*100
percent_nans2 = (sum(np.where(df1['call_2'] == 'N', 1, 0)))/len(df1['ref_2'])*100
percent_nans3 = (sum(np.where(df2['call_1'] == 'N', 1, 0)))/len(df2['ref_1'])*100
percent_nans4 = (sum(np.where(df2['call_2'] == 'N', 1, 0)))/len(df2['ref_2'])*100

print(percent_nans1, percent_nans2 ,percent_nans3 ,percent_nans4) #both cycles of both biochemistries return 2% N's, so this is not a factor that contributes to varying accuracies between two methods of sequencing



###generate summary chart for biochemistry 1
##find average signal intensity values for correct base calls

#subset the correct reads from cycle 1 of biochem 1
cycle1 = df1.iloc[:,1:6].loc[df1['mismatch_1'] == 0]
#rename columns to match those of cycle 2
cycle1.rename(columns={'A_1':'A', 'C_1':'C', 'G_1':'G', 'T_1':'T', 'call_1': 'call'}, inplace=True)

#subset the correct reads from cycle 2 of biochem 1
cycle2 = df1.iloc[:,8:13].loc[df1['mismatch_2'] == 0]
#renaame columns to match those of cycle 1
cycle2.rename(columns={'A_2':'A', 'C_2':'C', 'G_2':'G', 'T_2':'T', 'call_2': 'call'}, inplace=True)

#add the second cycle to the bottom of the first cycle
biochem1 = cycle1.append(cycle2, ignore_index=True)


#find the max signal intensity for all rows, representing the signal intensity for correct reads
biochem1['max'] = biochem1.max(axis=1)
#remove signal intensity reading columns
biochem1.drop(columns=['A', 'C', 'G', 'T'], inplace=True)
#find the mean of signal intensities for all four bases
biochem1['mean'] = biochem1.groupby(['call']).transform('mean')
#remove max signal intensity column, and remove duplicates
#only the average max intensity for each base is left
biochem1 = biochem1.drop(columns=['max']).drop_duplicates()
#sort alphabetically
biochem1.sort_values(by=['call'], inplace=True)
#rename call column to nucleotide
biochem1.rename(columns={'call':'nucleotide'}, inplace=True)


##find average signal intensity for incorrect base calls and incorrect base call counts
#subset the incorrect reads from cycle 1 of biochem 1
error1 = df1.iloc[:,0:5].loc[((df1['mismatch_1'] == 1) & (df1['call_1'] != 'N'))].reset_index(drop=True)
error1.rename(columns={'ref_1':'ref','A_1':'A', 'C_1':'C', 'G_1':'G', 'T_1':'T'}, inplace=True)

#subset the incorrect reads from cycle 2 of biochem 1
error2 = df1.iloc[:,7:12].loc[((df1['mismatch_2'] == 1) & (df1['call_2'] != 'N'))].reset_index(drop=True)
error2.rename(columns={'ref_2':'ref', 'A_2':'A', 'C_2':'C', 'G_2':'G', 'T_2':'T'}, inplace=True)

#add the second cycle to the bottom of the first cycle
biochem1_errors = error1.append(error2, ignore_index=True)

#function to return the signal intensity of the correct base call
def find_errors(df):
    ref = df.iloc[:,0]

    my_arr = []
    i = 0
    for j in ref:
        my_arr.append(df.loc[i,j])
        i += 1

    return(my_arr)

#create new column to contain signal intensity of intended base call
biochem1_errors['error'] = find_errors(biochem1_errors)
#remove signal intensity reading columns
biochem1_errors.drop(columns=['A', 'C', 'G', 'T'], inplace=True)
#find the mean of signal intensities for all four bases
biochem1_errors['mean'] = biochem1_errors.groupby(['ref']).transform('mean')
#count the number of errors for each nucleotide
biochem1_errors['error_counts'] = biochem1_errors.groupby(['ref'])['error'].transform('count')
#remove erroneous signal intensity column, and remove duplicates
#only the average max intensity for each base is left
biochem1_errors = biochem1_errors.drop(columns=['error']).drop_duplicates()
#sort alphabetically
biochem1_errors.sort_values(by=['ref'], inplace=True)
#rename call column to nucleotide
biochem1_errors.rename(columns={'ref':'nucleotide'}, inplace=True)


#merge two dataframes together to generate summary chart
summary1 = pd.merge(biochem1, biochem1_errors, on='nucleotide')
#count the total number of each nucleotides on the reference genome and append the total to summary chart
total_counts = df1['ref_1'].value_counts() + df1['ref_2'].value_counts()
summary1['total_counts'] = total_counts.values
summary1.to_csv('summary1.csv', index=False)
print(summary1)




###generate summary chart for biochemistry 2
##find average signal intensity values for correct base calls

#subset the correct reads from cycle 1 of biochem 2
cycle3 = df2.iloc[:,1:6].loc[df2['mismatch_1'] == 0]
cycle3.rename(columns={'A_1':'A', 'C_1':'C', 'G_1':'G', 'T_1':'T', 'call_1': 'call'}, inplace=True)

#subset the correct reads from cycle 2 of biochem 2
cycle4 = df2.iloc[:,8:13].loc[df2['mismatch_2'] == 0]
cycle4.rename(columns={'A_2':'A', 'C_2':'C', 'G_2':'G', 'T_2':'T', 'call_2': 'call'}, inplace=True)

#append second cycle to bottom of first cycle
biochem2 = cycle3.append(cycle4, ignore_index=True)

#find the max signal intensity for all rows, representing the signal intensity for correct reads
biochem2['max'] = biochem2.max(axis=1)
#remove signal intensity reading columns
biochem2.drop(columns=['A', 'C', 'G', 'T'], inplace=True)
#find the mean of signal intensities for all four bases
biochem2['mean'] = biochem2.groupby(['call']).transform('mean')
#remove max signal intensity column, and remove duplicates
#only the average max intensity for each base is left
biochem2 = biochem2.drop(columns=['max']).drop_duplicates()
#sort alphabetically
biochem2.sort_values(by=['call'], inplace=True)
#rename call column to nucleotide
biochem2.rename(columns={'call':'nucleotide'}, inplace=True)


##find average signal intensity for incorrect base calls and incorrect base call counts
#subset the incorrect reads from cycle 1 of biochem 2
error3 = df2.iloc[:,0:5].loc[((df2['mismatch_1'] == 1) & (df2['call_1'] != 'N'))].reset_index(drop=True)
error3.rename(columns={'ref_1':'ref','A_1':'A', 'C_1':'C', 'G_1':'G', 'T_1':'T'}, inplace=True)

#subset the incorrect reads from cycle 2 of biochem 2
error4 = df2.iloc[:,7:12].loc[((df2['mismatch_2'] == 1) & (df2['call_2'] != 'N'))].reset_index(drop=True)
error4.rename(columns={'ref_2':'ref', 'A_2':'A', 'C_2':'C', 'G_2':'G', 'T_2':'T'}, inplace=True)

#add the second cycle to the bottom of the first cycle
biochem2_errors = error3.append(error4, ignore_index=True)

#create new column to contain signal intensity of intended base call
biochem2_errors['error'] = find_errors(biochem2_errors)
#remove signal intensity reading columns
biochem2_errors.drop(columns=['A', 'C', 'G', 'T'], inplace=True)
#find the mean of signal intensities for all four bases
biochem2_errors['mean'] = biochem2_errors.groupby(['ref']).transform('mean')
#count the number of errors for each nucleotide
biochem2_errors['error_counts'] = biochem2_errors.groupby(['ref'])['error'].transform('count')
#remove erroneous signal intensity column, and remove duplicates
#only the average max intensity for each base is left
biochem2_errors = biochem2_errors.drop(columns=['error']).drop_duplicates()
#sort alphabetically
biochem2_errors.sort_values(by=['ref'], inplace=True)
#rename call column to nucleotide
biochem2_errors.rename(columns={'ref':'nucleotide'}, inplace=True)


#merge two dataframes together to generate summary chart
summary2 = pd.merge(biochem2, biochem2_errors, on='nucleotide')
#count the total number of each nucleotides on the reference genome and append the total to summary chart
total_counts = df2['ref_1'].value_counts() + df2['ref_2'].value_counts()
summary2['total_counts'] = total_counts.values
summary2.to_csv('summary2.csv', index=False)
print(summary2)

