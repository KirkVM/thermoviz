import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import altair as alt

#chart_data = pd.DataFrame(
#     np.random.randn(20, 3),
#     columns=['a', 'b', 'c'])
#st.slider('Kd',min_value=0.001,max_value=0.01)
Kd=st.sidebar.select_slider('Kd (mM)',options=[x for x in np.power(10,np.linspace(-6,0,100))],format_func=lambda x:f'{x*1e3:.2g}')
mtot=st.sidebar.select_slider('Protein Concentration (uM)',options=[x*1e-6 for x in np.power(10,np.linspace(-4,4,100))],
                               format_func=lambda x:f'{x*1e6:.2g}')
scale_type=st.sidebar.radio('scale',['linear','log'],index=1)
ltot_min=0
ltot_max=0.1
ltot=np.power(10,np.linspace(-6,-1,100))
bterm=-(mtot+ltot+Kd)
complex=(-bterm-np.sqrt(np.power(bterm,2)-4*ltot*mtot))/2
ligand=ltot-complex
macro=mtot-complex
frac_lbound=complex/mtot
concdf=pd.DataFrame({'complex':complex,
                     'ligand':ligand,
                     'protein':macro,
                     'ltot':ltot,
                     'frac_bound':frac_lbound})
melty=concdf.reset_index().melt(id_vars=['ltot'],value_vars=['complex','protein'],var_name='molecule',value_name='y')
melty.y=melty.y.apply(lambda x:x*1e6)
chart=alt.Chart(melty)
lchart=chart.mark_line().encode(x=alt.X('ltot',title='Ligand Concentration (mM)',scale=alt.Scale(type=scale_type)),
                       y=alt.Y('y',title='Concentration (uM)',axis=alt.Axis()),
                       color='molecule')#,tooltip=alt.Tooltip('hi','complex'))
thetachart=alt.Chart(concdf).mark_line(color='black').\
           encode(alt.X('ltot',title='Ligand Concentration (mM)',scale=alt.Scale(type=scale_type)),
                  alt.Y('frac_bound',title='Theta (fraction protein bound)'))

st.altair_chart(thetachart,use_container_width=True)
st.altair_chart(lchart,use_container_width=True)


#chart=alt.layer(lchart,pchart)
#st.altair_chart(chart, use_container_width=True)

#chart=alt.Chart(concdf)
#lchart=chart.mark_line(color='blue').encode(alt.X('ltot',scale=alt.Scale(type=scale_type)),
#                       y=alt.Y('complex',title='Concentration (uM)'))#,tooltip=alt.Tooltip('hi','complex'))







#chart=lchart+pchart
#chart.encode(alt.Y())
#cchart=chart.mark_line().encode(x='ltot',y='complex')
#st.altair_chart(cchart, use_container_width=True)
#st.altair_chart(pchart, use_container_width=True)




#plt.plot('ltot','complex',data=concdf)
#plt.show()
#ltot_conc=[x for x in range(100)]
#ptot_conc=np.linspace(0,1,100)
#bound_conc=[15 for x in range(100)]
#concdict={'ltot_conc':ltot_conc,'ptot_conc':ptot_conc,'bound_conc':bound_conc}
#concdf=pd.DataFrame(concdict)
#
#
#st.line_chart(concdf)
#st.slider()