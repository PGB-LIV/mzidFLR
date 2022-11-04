#Calculate false detection rate statistics and plots
#-Calls Extract_DF
#	Loads CSV file to dataframe for FDR calculations
#Writes folder for downstream analysis eg. "FDR_0.01_PTM_score_0" FDR filter 0.01, no PTM score filter
#output="FDR_output.csv","FDR.jpg","FDR_filter.jpg","FDR_score.jpg"

	
import TPP_reusable.Extract_DF as Extract_DF
import numpy as np
import matplotlib.pylab as plt
import pandas as pd

#calculate FDR
def calculateFDR(results_file,output,PXD,mod):
    #extract results to df
    df=Extract_DF.extract_PTMprophet_IDent_df(results_file,PXD,mod)
    df['Score']=df['Score'].astype(float)
    df=df.sort_values(by='Score',ascending=False)
    df=df.reset_index(drop=True)
    #set decoy and target
    df['Protein_count'] = df['Protein'].str.count(":")+1
    df['Decoy_count'] = df['Protein'].str.count("DECOY")
    df['Decoy'] = np.where(df['Protein_count']==df['Decoy_count'],1,0)
    df=df.drop(columns=['Protein_count', 'Decoy_count'])
    #target_count = row_count-decoy_count
    df['decoy_count']=df['Decoy'].cumsum()
    df['row_count']=(df.index)+1
    df['target_count']=df['row_count']-df['decoy_count']
    #FDR = decoy_count/target_count
    df['FDR']= df['decoy_count']/df['target_count']
    df['FDR'] = df['FDR'].astype(float)

    #warning for no decoys
    if df['FDR'].max==0:
        sys.exit("No decoys found - check decoy prefix")

    #q_val = min FDR at this position or lower
    df['q_value']=df['FDR']
    df['q_value']=df.iloc[::-1]['FDR'].cummin()

    #return as csv
    df.to_csv(output,index=False)

    #FDR plots
    ax1=df.plot.scatter(x='FDR',y='row_count')
    ax2=df.plot.line(x='q_value',y='row_count')
    plt.ylabel("PSM count")
    plt.savefig('FDR.jpg',dpi=300)

    df.plot.line(x='Score', y='FDR',ylim=(0,1),xlim=(1,0))
    plt.xlabel("Peptide Probabilty")
    plt.ylabel("Global FDR")
    plt.savefig('FDR_score.jpg',dpi=300)
    plt.close('all')


#precursor mass tolerance plot
def ppm_error(results_file):
    df = pd.read_csv(results_file)
    df=df.loc[df['mass_diff'] <=0.9]
    df.plot.scatter(x='Score', y='ppm_error', s=0.5)
    plt.ylabel("PPM error")
    plt.xlabel("PSM probability")
    plt.savefig('PPM_error.jpg', dpi=300)
    df.to_csv("error.csv", index=False)

    df = df.loc[df['q_value'] <= 0.01]
    df=df.loc[df['mass_diff'] <=0.9]
    df.plot.scatter(x='Score', y='ppm_error', s=0.5)
    plt.ylabel("PPM error")
    plt.xlabel("PSM probability")
    plt.savefig('PPM_error_FDR.jpg', dpi=300)
    df.to_csv("error_FDR.csv", index=False)