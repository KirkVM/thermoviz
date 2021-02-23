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
#def fluor_counts_predict(Kd,bound_fluor,bline_fluor,pconc,sconc):
#    bTerm=1 - pconc/sconc - Kd/sconc #treat protein as "ligand"
#    squareTerm=np.power(bTerm,2) - 4*pconc*sconc
#    theta=0.5*(-bTerm - np.sqrt(squareTerm))
#    fluor=bline_fluor + bound_fluor*(1-theta)
#    return fluor
#
#def pulldown_l2lossfunc(kdatup,counts,pconc,sconcs):#
#    Kd,bound_fluor,bline_fluor=kdatup
#    pred_counts=[]
#    for sconc in sconcs:
#        pred_counts.append(fluor_counts_predict(Kd,bound_fluor,bline_fluor,pconc,sconc))
#    delta = [pred_counts - y for pred_counts,y in zip(pred_counts,counts)]
#    return np.dot(delta,delta)

def predict_fluorcounts(Kd,bound_fluor,bline_fluor,pconc,sconc):
    bTerm= - (pconc + sconc + Kd)
    squareTerm=np.power(bTerm,2) - 4*pconc*sconc
    complex= 0.5*(-bTerm - np.sqrt(squareTerm))
    theta=complex/pconc
    pred_counts=bline_fluor + bound_fluor * (1-theta)
    return pred_counts

def lm_pdminimizer(params,concs,fluor_counts,pconc):
    pred_counts=predict_fluorcounts(params['Kd'],params['bound_fluor'],params['bline_fluor'],pconc,concs)
    return fluor_counts-pred_counts

st.sidebar.title("Pull-down Binding Fit")
data_file = st.sidebar.file_uploader('Upload a data file')
if data_file is not None:
    datadf=pd.read_csv(data_file,names=['conc','counts'])
    assert (datadf.shape[0]>1 and datadf.shape[1]==2)
else: #for debugging
    datadf=pd.read_csv('data/example_pulldown.txt',names=['conc','counts'])

#Kd=st.sidebar.select_slider('Kd (mL/mg)',options=[x for x in np.power(10,np.linspace(-6,0,100))],format_func=lambda x:f'{x*1e3:.2g}')
col1,col2=st.sidebar.beta_columns([4,1])
with col1:
    Kd=st.number_input('Kd (mg/mL)',value=2.3,min_value=0.001,max_value=100.)
with col2:
    Kd_hold=st.checkbox('Kd Hold?')
col1,col2=st.sidebar.beta_columns([4,1])
with col1:
    bound_fluor=st.number_input('protein fluorescence',value=3500.,min_value=2000.,max_value=5000.) 
with col2:
    protfluor_hold=st.checkbox('Fluor Hold?')
col1,col2=st.sidebar.beta_columns([4,1])
with col1:
    bline_fluor=st.number_input('baseline fluorescence',value=500.,min_value=0.,max_value=1000.) 
with col2:
    blinefluor_hold=st.checkbox('Bline Hold?')
#Kd=st.sidebar.number_input('Kd (mg/mL)',value=2.3,min_value=0.001,max_value=100.) 
#Kd_hold=st.sidebar.checkbox('Hold?')

params=Parameters()
params.add('Kd',value=1.,min=.001,max=100)
params.add('bound_fluor',value=3000,min=2000,max=5000)
params.add('bline_fluor',value=500,min=0,max=1000)
lmpdminner=Minimizer(lm_pdminimizer,params,fcn_args=(datadf.conc.values,datadf.counts.values,0.5e-6))
result=lmpdminner.minimize()
fitconcs=np.linspace(0,10,50)
fitKd=result.params['Kd'].value

#    bound_fluor=st.number_input('protein fluorescence',value=3500.,min_value=2000.,max_value=5000.) 
#    bline_fluor=st.number_input('baseline fluorescence',value=500.,min_value=0.,max_value=1000.) 
#    protfluor_hold=st.checkbox('Fluor Hold?')
#    blinefluor_hold=st.checkbox('Bline Hold?')
#Kd=st.sidebar.number_input('Kd (mg/mL)',value=2.3,min_value=0.001,max_value=100.) 
#Kd_hold=st.sidebar.checkbox('Hold?')

params=Parameters()
params.add('Kd',value=1.,min=.001,max=100)
params.add('bound_fluor',value=3000,min=2000,max=5000)
params.add('bline_fluor',value=500,min=0,max=1000)
lmpdminner=Minimizer(lm_pdminimizer,params,fcn_args=(datadf.conc.values,datadf.counts.values,0.5e-6))
result=lmpdminner.minimize()
fitconcs=np.linspace(0,10,50)
fitKd=result.params['Kd'].value

fitfluors=predict_fluorcounts(result.params['Kd'].value,
                              result.params['bound_fluor'].value,
                              result.params['bline_fluor'].value,
                              0.5e-6,fitconcs)
fitdf=pd.DataFrame({'conc':fitconcs,'fitcounts':fitfluors})

dchart=alt.Chart(datadf).mark_circle(size=100).encode(alt.X('conc',scale=alt.Scale(type='linear')),
                                                      alt.Y('counts')
                                                      )
fchart=alt.Chart(fitdf).mark_line().encode(alt.X('conc'),alt.Y('fitcounts'))
xchart=alt.layer(dchart,fchart)
st.altair_chart(xchart,use_container_width=True)
#st.text_area(report_fit(result))
dumpstr=f"Kd={result.params['Kd'].value}\n"
dumpstr+=f'protein fluor. counts={result.params["bound_fluor"].value}\n'
dumpstr+=f'bline fluor.={result.params["bline_fluor"].value}\n'
st.text(dumpstr)


#combodf=pd.concat([datadf,fitdf])
#melty=pd.melt(combodf,id_vars=['conc'])
#chart=alt.Chart(melty)
#cchart=chart.mark_circle(size=100).encode(alt.X('conc'),alt.Y('value',field='counts'))
#st.altair_chart(cchart)



#exdf=pd.DataFrame({'predconc':[0,3,10],'predcounts':[3800,3100,1500]})
#exchart=alt.Chart(exdf)
#exchart.mark_line().encode(alt.X('predconc'),alt.Y('predcounts'))
#
#st.altair_chart(lchart+exchart)#,use_container_width=True)
#st.line_chart(df)
#st.beta_expander

#if data_file and st.sidebar.button('Load'):
#    image = get_data_file(data_file)
#    with st.beta_expander('Selected Image', expanded = True):
#        st.image(image, use_column_width = True)

