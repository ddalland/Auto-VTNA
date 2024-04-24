import pandas as pd
import matplotlib
matplotlib.use('Qt5Agg')
#import matplotlib.font_manager
#import matplotlib.pyplot as plt
data_VTNA2=pd.read_excel('data_VTNA2_aa.xlsx',sheet_name=None)
bull_data=pd.read_excel('Bull_CH_data_2.xlsx',sheet_name=None)
#data_THF=pd.read_excel('THF_30.xlsx',sheet_name=None)
#import auto_vtna
from __init__ import plot_data
from __init__ import plot_data_together
from __init__ import plot_data_MB
#from Normal_VTNA import Normal_VTNA
#from auto_vtna.Normal_VTNA import Normal_VTNA
#from Automatic_VTNA import Automatic_VTNA
#selection_VTNA2_a={"RM":{"RM5":{'omissions':[6,5]},"hallo_hallo_hallo":{'omissions':[5,6,4]},"RM3":{'omissions':[4,5,6]}}, "normalised_species":{'RBr':1,'4-Me-Aniline':1}, "output_species":'Prod'}
#auto_bull1=Automatic_VTNA(bull_data,bull_d1,constraint=None,multi_threading=True)
#auto_bull1.plot_orders_vs_overlay()
# del data_THF['THF_30_1']
#del data_VTNA2['RM6']
# del data_VTNA2['RM1']
#vtna=Normal_VTNA(data_VTNA2,selection_VTNA2_a)
#vtna.plot_VTNA_with_overlay_score(DP_scaler=1,xtick_rotation=20,constraint='j',deg=1,fit_metric='RMSE',show_legend=True,legend_outside=True,size_scaler=1.2)
#vtna.plot_VTNA_with_overlay_score(constraint='origin',show_legend=True,legend_outside=True)
#vtna.plot_VTNA(legend_outside=True,show_legend=False)
#vtna.plot_VTNA_with_overlay_score(size_scaler=0.9,constraint=None,deg=2,fit_metric='through_origin',show_fit_function=True,show_overlay_score=False)
#plt.show()
plot_data_MB(data_VTNA2,['RBr','Styrene', 'Prod'],[1,1,1],legend_outside=True,plot_mode='scrollable',linewidth=0,DP_scaler=2,ylim=2)
plot_data(data_VTNA2,plot_mode='together',legend_outside=False,ylim=2)
plot_data(bull_data,plot_mode='together',legend_outside=True)
#plot_data_MB(data_VTNA2,['RBr','Prod'],[1,1],plot_mode='scrollable',fig_size_scaler=1,legend_outside=True)
#plot_data_MB(data_VTNA2,['RBr','Prod'],[1,1],plot_mode='sc',fig_size_scaler=1,legend_outside=True)
#plot_data_MB(data_THF,['A','P'],[1,1],plot_mode='together',fig_size_scaler=1,legend_outside=False)
#plot_data(data_VTNA2,plot_mode='separate',legend_outside=True)
plot_data(data_VTNA2,plot_mode='together',legend_outside=True,linewidth=0,DP_scaler=2)
#plot_data(data_VTNA2,plot_mode='scrf',legend_outside=True)
#plot_data(data_THF,plot_mode='together')
#plot_data_together(data_VTNA2,['RBr'],legend_outside=True)
#plot_data_2(data_VTNA2,plot_mode='scrollable',fig_size=1)
#vtna2_a=Automatic_VTNA(data_VTNA2,selection_VTNA2_a,constraint=None,fixed_order_species='RBr')
#vtna2_a=Automatic_VTNA(data_VTNA2,selection_VTNA2_a,constraint='monotonic')
#vtna2_b=Automatic_VTNA(data_VTNA2,selection_VTNA2_a,constraint=None,fixed_order_species='RBr',deg=7)
#vtna2_c=Automatic_VTNA(data_VTNA2,selection_VTNA2_a,fixed_order_species='RBr',deg=1,constraint='origin')
#vtna2_a.plot_orders_vs_overlay(interval=True,points=True)
#vtna2_b.plot_orders_vs_overlay(interval=True)
#vtna2_c.plot_orders_vs_overlay(interval=True)
#vtna2_analysis_2=Automatic_VTNA(data_VTNA2,selection_VTNA2_a,score_interval=0.3,deg=5,constraint=None,fixed_order_species='RBr')