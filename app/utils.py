import plotly
import plotly.graph_objs as go
import json
import re

# The following method reutrns a histogram of list passed in JSON form.
def create_plot(data, normalized=False):

    data = [
        go.Histogram(
            x=data['Length'],
            histnorm = 'probability' if normalized else None
        )       
    ]

    graphJSON = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)

    return graphJSON

    
# The following method filters the data to remove contamination.
# This method is specific to a PEAKS output file for peptides.
def filterPeaksFile(data, dropPTM=True, minLen=1, maxLen=133):
    
#   Removing rows with no accession identity
    data = data.dropna(subset=['Accession'])
    
#   Dropping the peptides with Post translational modifications
    if dropPTM:
        data = data[data.apply(lambda x : re.search(r'[(].+[)]',x['Peptide']) == None,axis=1)]

#   Removing contamincation founf from accession number 
    data = data[data.apply(lambda x : str(x['Accession']).find('CONTAM') == -1,axis=1)]
    
#   Filtering on the basis of the peptide length
    data = data[data.apply(lambda x : x['Length'] in range(minLen,maxLen),axis=1)]
    
    return data