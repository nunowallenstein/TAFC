import numpy as np
import pandas as pd
import csv
#file reading
datasettrain=pd.read_csv('training.csv',header=None,low_memory=False)
datasettest=pd.read_csv('test.csv',header=None,low_memory=False)

header=datasettrain.iloc[0,:]
#header=[datasettrain.iloc[0,0],'PC1','PC2','PC3','PC4',datasettrain.iloc[0,-2],datasettrain.iloc[0,-1]]


weights_train=datasettrain.iloc[1:,-2]
Y_train=datasettrain.iloc[1:,-1]
Event_train=datasettrain.iloc[1:,0]
X_train=datasettrain.iloc[1:,1:-2].values
header_test=datasettest.iloc[0,:]

#header_test=[datasettest.iloc[0,0],'PC1','PC2','PC3','PC4']
X_test=datasettest.iloc[1:,1:].values
Event_test=datasettest.iloc[1:,0]

#to normalize the data so that we have a faster convergence

from sklearn.preprocessing import StandardScaler
sc =StandardScaler()
X_train=sc.fit_transform(X_train)
X_test=sc.transform(X_test)
"""
from sklearn.decomposition import PCA
pca=PCA(n_components=4)
X_train=pca.fit_transform(X_train)
X_test=pca.transform(X_test)
explained_variance=pca.explained_variance_ratio_
"""


df1=pd.DataFrame(Event_train).reset_index(drop=True)
df2=pd.DataFrame(X_train)
df3=pd.DataFrame(weights_train).reset_index(drop=True)
df4=pd.DataFrame(Y_train).reset_index(drop=True)

df5=pd.concat([df1,df2,df3,df4],axis=1,ignore_index=True)
df6=pd.DataFrame(header)
Lista_train=pd.concat([df6.T,df5],join="outer", ignore_index=True)


Lista_train.to_csv('proctrain.csv',index=False,index_label=False,header=False, float_format="%.7f",decimal='.')

df7=pd.DataFrame(X_test)
df8=pd.DataFrame(Event_test).reset_index(drop=True)
df9=pd.concat([df8,df7],axis=1,ignore_index=True)
df10=pd.DataFrame(header_test)
Lista_test=pd.concat([df10.T,df9],join="outer",ignore_index=True)
Lista_test.to_csv('proctest.csv',index=False,index_label=False,header=False, float_format="%.7f",decimal='.')

