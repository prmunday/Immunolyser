import plotly
import plotly.graph_objs as go
import json
import re
import pandas as pd
from numpy import std
import os

# The following method reutrns a histogram of list passed in JSON form.
def plot_lenght_distribution(samples, hist="percent"):

    # Initializing plotly figure
    fig = go.Figure()

    # Appeding all samples
    for sample_name,sample in samples.items():

        # len_data = sample.getCombniedData()['Length']

        data = []
        
        for replicateData in sample.values():
            data.append(replicateData)
        
        len_data = pd.concat(data)['Length']

        peptideProportion = {}

        if hist=='percent':
            # Storing the proportions of each n-mer
            for replicate, data in sample.items():
                peptideProportion[replicate] = data.groupby('Length').count()['Peptide']/data.shape[0]*100

            title = 'The frequency distribution of the peptide lengths'
            yaxis_label = '% Peptides'
        else:
            for replicate, data in sample.items():
                peptideProportion[replicate] = data.groupby('Length').count()['Peptide']

            title = 'The density distribution of the peptide lengths'
            yaxis_label = 'Numebr of Peptides'

        # Combining arrays to further calculate the standard deviation
        errors = pd.concat(peptideProportion, axis=1).apply(lambda x : std(x), axis=1)



        fig.add_trace(go.Histogram(
                            x=len_data,
                            histnorm = hist,
                            name = sample_name,
                            error_y=dict( 
                                type='data', 
                                # array=sample.getPeptideLengthError()/2,
                                array=errors/2,
                                color='green', 
                                thickness=1, 
                                width=3, 
                            )
                        )
        )

    # Adding labels in the chart
    fig.update_layout(
            title= title,
            xaxis = dict(title='<i>Length</i>'),
            yaxis = dict(title='<i>'+ yaxis_label +'</i>')
        )
    
    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

    return graphJSON

    
# The following method filters the data to remove contamination.
# This method is specific to a PEAKS output file for peptides.
def filterPeaksFile(samples, dropPTM=True, minLen=1, maxLen=133):
    
    for file_name,sample in samples.items():
#       Removing rows with no accession identity
        temp = sample.dropna(subset=['Accession'])
        
#       Dropping the peptides with Post translational modifications
        if dropPTM:
            temp = temp[temp.apply(lambda x : re.search(r'[(].+[)]',x['Peptide']) == None,axis=1)]

#       Removing contamincation founf from accession number 
        temp = temp[temp.apply(lambda x : str(x['Accession']).find('CONTAM') == -1,axis=1)]
        
#       Filtering on the basis of the peptide length
        temp = temp[temp.apply(lambda x : x['Length'] in range(minLen,maxLen),axis=1)]
    
        samples[file_name] = temp

    return samples

def saveNmerData(location, samples, peptideLength = 9):

    for file_name, data in samples.items():
        for replicate_name, replicate_data in data.items():
            replicate_data[replicate_data.Length == 9]['Peptide'].to_csv(os.path.join(location, file_name, replicate_name[:-4]+'.txt'), header=False, index=False)


def getSeqLogosImages(samples_data):

    seqlogos = {}

    for sample,replicates in samples_data.items():
        seqlogos[sample] = [replicate[:-4]+'-001.jpg' for replicate in replicates.keys()]

    return seqlogos