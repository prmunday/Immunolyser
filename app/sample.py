import pandas as pd

# The folloing class will be 'Sample' class. This will contain all information about one sample.
# This will also keep the replicated, if any, of the sample in one place.
# This will also have the methods associated with a single sample
class Sample:
    def __init__(self, name,data={}):
        assert isinstance(data, dict), 'Sample should be of dictionary type. Where key is the name of the replicate and value is the data itself'
        self.data = data
        self.nReplicates = len(data)
        self.name = name
    
#   Method to add replicates
    def addReplicate(self, name, data):
        if name not in self.data.keys():
            self.data[name] = data
            self.nReplicates += 1
        else:
            print('Replicate {} already exists in the sample'.format(name))
    
#   Method to return combined data
    def getCombniedData(self):
        data = []
        
        for replicateData in self.data.values():
            data.append(replicateData)
        
        return(pd.concat(data))

#   This method is to get min/max proportion of n-mers for the purpose of length distribution
#   This method is exclusively fot PEAKS data file
    def getPeptideLengthError(self):

        # This dictionary is to store grouped data for every replicate
        peptideProportion = {}

        # Storing the proportions of each n-mer
        for replicate, data in sample.items():
            peptideProportion[replicate] = data.groupby('Length').count()['Peptide']/data.shape[0]*100

        # Combining arrays to further calculate the standard deviation
        errors = pd.concat(peptideProportion, axis=1).apply(lambda x : std(x), axis=1)

        return(errors)