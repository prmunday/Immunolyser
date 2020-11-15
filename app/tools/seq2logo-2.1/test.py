
# import plotly.offline as pyo

# import plotly.graph_objs as go
# pyo.init_notebook_mode()

import pandas as pd
import numpy as np
# import re

import seqlogo as sl

psm = pd.read_table('../tools/seq2logo/seq2logo-2.1.all/seq2logo-2.1/output.txt',skiprows=2, sep=' ')

np.random.seed(42)

random_ppm = np.random.dirichlet(np.ones(4), size=6)
ppm = sl.Ppm(random_ppm)
sl.seqlogo(ppm, ic_scale = False, format = 'svg', size = 'medium')