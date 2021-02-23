import streamlit as st
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
import altair as alt
from lmfit import Minimizer,Parameters, report_fit



#def get_data_file(dfpath):
#    with open(dfpath,'r') as f:
#        return f.read()
#
def fluor_counts_predict(Kd,bound_fluor,bline_fluor,pconc,sconc):
    bTerm=1 - pconc/sconc - Kd/sconc #treat protein as "ligand"
    squareTerm=np.power(bTerm,2) - 4*pconc*sconc
    theta=0.5*(-bTerm - np.sqrt(squareTerm))
    fluor=bline_fluor + bound_fluor*(1-theta)
    return fluor

def pulldown_l2lossfunc(kdatup,counts,pconc,sconcs):#
    Kd,bound_fluor,bline_fluor=kdatup
    pred_counts=[]
    for sconc in sconcs:
        pred_counts.append(fluor_counts_predict(Kd,bound_fluor,bline_fluor,pconc,sconc))
    delta = [pred_counts - y for pred_counts,y in zip(pred_counts,counts)]
    return np.dot(delta,delta)

def lm_pdminimizer(params,concs,fluor_counts,pconc):
    bTerm= - (pconc + concs + params['Kd'])
    squareTerm=np.power(bTerm,2) - 4*pconc*concs
    complex= 0.5*(-bTerm - np.sqrt(squareTerm))
    theta=complex/pconc
    pred_counts=params['bline_fluor'] + params['bound_fluor']*(1-theta)
    return fluor_counts-pred_counts

st.sidebar.title("Pull-down Binding Fit")
data_file = st.sidebar.file_uploader('Upload a data file')
if data_file is not None:
    datadf=pd.read_csv(data_file,names=['conc','counts'])
    assert (datadf.shape[0]>1 and datadf.shape[1]==2)
else: #for debugging
    datadf=pd.read_csv('data/example_pulldown.txt',names=['conc','counts'])

scale_type=st.sidebar.radio('scale',['linear','log'],index=1)

#melty=df.reset_index().melt(id_vars=['conc'],value_vars=['counts'])
#chart=alt.Chart(melty)
#
#lchart=chart.mark_line().encode(x=alt.X('conc',title='Ligand Concentration (mM)',scale=alt.Scale(type=scale_type)),
#                       y=alt.Y('counts',title='Concentration (uM)',axis=alt.Axis()))
#                       color='molecule')#,tooltip=alt.Tooltip('hi','complex'))
params=Parameters()
params.add('Kd',value=1)
params.add('bound_fluor',3000)
params.add('bline_fluor',500)
lmpdminner=Minimizer(lm_pdminimizer,params,fcn_args=(datadf.conc.values,datadf.counts.values,0.5e-6))
result=lmpdminner.minimize()

chart=alt.Chart(datadf)
lchart=chart.mark_circle(size=100).encode(alt.X('conc'),alt.Y('counts'))
exdf=pd.DataFrame({'predconc':[0,3,10],'predcounts':[3800,3100,1500]})
exchart=alt.Chart(exdf)
exchart.mark_line().encode(alt.X('predconc'),alt.Y('predcounts'))

st.altair_chart(lchart+exchart)#,use_container_width=True)
#st.line_chart(df)
#st.beta_expander

#if data_file and st.sidebar.button('Load'):
#    image = get_data_file(data_file)
#    with st.beta_expander('Selected Image', expanded = True):
#        st.image(image, use_column_width = True)
