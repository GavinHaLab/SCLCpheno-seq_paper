import sys, re, time
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# compared to v2:
# allowed specifying number of trim steps to perform (will usually run with 0 since i have just been using the 720 bp window size anwyays)
# compared to v1:
# adding rand index computation
# automatically making cluster 1 the "on" cluster and cluster 1 the "off" cluster

_,tfbs_profile_fn,sample_label_fn,num_trim_steps,out_file, = sys.argv

tmp = pd.read_csv(tfbs_profile_fn,index_col=0,dtype={'patient_id':'string','subtype':'object'})
tmp = tmp.set_index([ x for x in tmp.columns if not re.match('^-{0,1}[0-9]*$', x) ])
tmp.columns = tmp.columns.astype(int)
tmp = tmp.sort_index(axis=1)
max_window_size = 720
tmp = tmp.reindex( columns=np.arange(-int(max_window_size/2),int(max_window_size/2),15) ) # will index from -360 to +345 which is desired behavior because windows are left-indexed
print('loaded TFBS profiles')

labels = pd.read_csv(sample_label_fn).set_index('unique_library_id') * 1

l = []
km = KMeans(n_clusters=2, random_state=0)
histo_fltr_list = [ (lambda df:df,'all_pdx'), (lambda df:df.reindex(['SCLC','sclc'],level='histology',axis=0),'sclc_only') ]
for histo_fltr,histo_fltr_name in histo_fltr_list:
    g1 = histo_fltr(tmp)
    for trim in range(int(num_trim_steps)+1):
        start_time = time.time()
        ctr = 0
        g2 = g1
        if trim > 0:
            g2 = g2.reindex( columns=g2.columns[2*trim:-2*trim] )
        num_features = g2.shape[1]
        for (site_name,site_id),g3 in g2.groupby(['site_name','unique_ID']):
            tf = site_name.split('.')[0]
            g3 = g3[ ~g3.index.get_level_values(tf).isna() ]
            g3 = g3[ g3.var(axis=1)>0 ] # remove samples with invariant coverage profiles
            g3 = g3.dropna(axis=1,how='all') # remove positions that are nan for all samples (aligned REST sites?)
            if min(g3.shape)<3:
                continue
            exp = g3.index.get_level_values(tf).values
            this_labels = labels.reindex(g3.index.get_level_values('unique_library_id').tolist())
            num_features_actual = g3.shape[1]
            num_samples = g3.shape[0]
            pred = km.fit_predict(g3)
            if exp[pred==0].mean() > exp[pred==1].mean():
                pred[ pred == 0 ] = -1
                pred[ pred == 1 ] = 0
                pred[ pred == -1 ] = 1
            _,pvalue = mannwhitneyu( exp[pred==0], exp[pred==1] )
            result =  pd.Series( { 'site_name':site_name,
                                    'unique_ID':site_id,
                                    'histo_fltr':histo_fltr_name, 
                                    'num_samples':num_samples, 
                                    'num_features':num_features, 
                                    'num_features_actual':num_features_actual,
                                    'exp_min':exp.min(), 
                                    'exp_max':exp.max(), 
                                    'exp_var':exp.var(), 
                                    'pred_0_num':(pred==0).sum(), 
                                    'pred_0_mean':exp[pred==0].mean(), 
                                    'pred_1_mean':exp[pred==1].mean(),
                                    'adj_rand_index': adjusted_rand_score( this_labels[tf].dropna(), pred[this_labels[tf].notna()] ), # drop na ground truth labels when tabulating rand index...
                                    'pvalue':pvalue } )
            result = pd.concat( [ result, pd.Series( pred, index=g3.index.get_level_values('unique_library_id') ) ] )
            l.append(result)
            ctr+=1
            if (ctr)%100==0:
                print('processed {:,d} genes in {:.0f} seconds for histo_fltr_name={} and trim={}...'.format(ctr,time.time()-start_time,histo_fltr_name,trim),end='\r') # 
        print('processed {:,d} genes in {:.0f} seconds for histo_fltr_name={} and trim={}...DONE'.format(ctr,time.time()-start_time,histo_fltr_name,trim))
    print('done with histo_fltr_name={}...'.format(histo_fltr_name))
    
tmp = pd.concat( l, axis=1 ).T
tmp = pd.concat( [ pd.concat([g1,pd.Series(multipletests(g1['pvalue'])[1],index=g1.index,name='pvalue_adj')],axis=1) for n1,g1 in tmp.groupby(['histo_fltr','num_features']) ] )
tmp.to_csv(out_file)
print('DONE')