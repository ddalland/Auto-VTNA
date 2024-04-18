from auto_vtna.Normal_VTNA import Normal_VTNA
data_VTNA2=pd.read_excel('data_VTNA2_aa.xlsx',sheet_name=None)
selection_VTNA2_a={'RM': {'RM4': {'omissions': [5, 6]},
  'RM5': {'omissions': [6, 5]},
  'RM2': {'omissions': [5, 6, 4]},
  'RM3': {'omissions': [4, 5, 6]},
  'RM6': {'omissions': [6]}},
 'normalised_species': {'RBr': 1.033,
  '4-Me-Aniline': 1.043},
 'output_species': 'Styrene'}
vtna=Normal_VTNA(data_VTNA2,selection_VTNA2_a)
vtna.plot_VTNA()