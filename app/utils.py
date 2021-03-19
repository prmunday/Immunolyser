import plotly
import plotly.graph_objs as go
import json
import re
import pandas as pd
from numpy import std
import os
import glob
from subprocess import call, Popen
from distutils.dir_util import copy_tree
import shutil
from zipfile import ZipFile
from os.path import basename
from app import app

project_root = os.path.dirname(os.path.realpath(os.path.join(__file__, "..")))
data_mount = app.config['IMMUNOLYSER_DATA']

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
    
#     control_peptides = list()
#     # Filtering the control data first
#     for file_name,sample in control_data.items():
# #       Removing rows with no accession identity
#         temp = sample.dropna(subset=['Accession'])
        
# #       Dropping the peptides with Post translational modifications
#         if dropPTM:
#             temp = temp[temp.apply(lambda x : re.search(r'[(].+[)]',x['Peptide']) == None,axis=1)]

# #       Removing contamincation founf from accession number 
#         temp = temp[temp.apply(lambda x : str(x['Accession']).find('CONTAM') == -1,axis=1)]
        
# #       Filtering on the basis of the peptide length
#         temp = temp[temp.apply(lambda x : x['Length'] in range(minLen,maxLen),axis=1)]
    
#         control_data[file_name] = temp

#     # Finding the control peptides from control data
#     for control_replicate, control_replicate_data in control_data.items():
#         control_peptides.extend(control_data[control_replicate]['Peptide'].to_list())

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

#       Filtering the control peptides out
        # temp = temp[temp.apply(lambda x : x['Peptide'] not in control_peptides,axis=1)]        
    
        samples[file_name] = temp

    return samples

def saveNmerData(location, samples, peptideLength = 9):

    for file_name, data in samples.items():
        for replicate_name, replicate_data in data.items():

            if type(peptideLength) == int:
                replicate_data[replicate_data.Length == peptideLength]['Peptide'].to_csv(os.path.join(location, file_name, replicate_name[:-4]+'_'+str(peptideLength)+'mer.txt'), header=False, index=False)
            else:
                replicate_data[replicate_data['Length'].between(peptideLength[0], peptideLength[1], inclusive=True)]['Peptide'].to_csv(os.path.join(location, file_name, replicate_name[:-4]+'_'+str(peptideLength[0])+'to'+str(peptideLength[1])+'mer.txt'), header=False, index=False)


def getSeqLogosImages(samples_data):

    seqlogos = {}

    # This approach has to be modified as the names of logos are derived from the data input,
    # not from the seq2logo results.
    for sample,replicates in samples_data.items():
        seqlogos[sample] = [replicate[:-4]+'-001.jpg' for replicate in replicates.keys()]

    return seqlogos

def getGibbsImages(taskId, samples_data):

    gibbsImages = {}

    # This approach has to be modified as the cluster are picked from the files(JPG) present in results.
    # It should be linked with gibbscluster directly to get the results.
    for sample,replicates in samples_data.items():
        # seqlogos[sample] = [replicate[:-4]+'-001.jpg' for replicate in replicates.keys()]

        gibbsImages[sample] = dict()

        for replicate in replicates.keys():
            bar_plot = [os.path.basename(x) for x in glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate[:-4]}/images/*.JPG')]
            clusters = [os.path.basename(x) for x in glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate[:-4]}/logos/*.jpg')]
            
            gibbsImages[sample][replicate[:-4]] = dict()
            gibbsImages[sample][replicate[:-4]][bar_plot[0]] = clusters
            print(gibbsImages)

            # gibbsImages[sample][replicate][bar_plot] = clusters

    return gibbsImages

#def setUpWsl():
    #print('Setting up the WSL(Windows ) for windows system.')
    
    # Following script is to set up a WSL in windows system to linux platform dependent tools.
    # This is not automated yet.

# This method will generate the binding predictions.
# User can select the binding prediction to be used.
# User can enter the names of alleles of interest.
def generateBindingPredictions(taskId, alleles, method):
    
    # Reading the alleles which can be accepted by Anthem
    if method=='ANTHEM':
        with open('app/tools/Anthem-master/source/lenghHLA.json') as f:
            data = json.load(f)

    print('Generating Binding Predictions for task {} for {} alleles.'.format(taskId,alleles))

    # Converting HLAs syntax from HLA-A02:01 to HLA-A*02:01
    temp = list()
    allelesForAnthem = ""
    for allele in alleles.split(','):
        temp.append('HLA-{}*{}:{}'.format(allele[4],allele[5:7],allele[8:]))
        
    allelesForAnthem = ",".join(temp)
    del temp

    for sample in os.listdir('{}/{}'.format(data_mount,taskId)):
        for replicate in os.listdir('{}/{}/{}'.format(data_mount,taskId,sample)):
            if sample != 'Control':
                if replicate[-12:] == '8to14mer.txt':
                    if(method=='MixMHCpred'):
                        call(['./app/tools/MixMHCpred/MixMHCpred', '-i', '{}/{}/{}/{}'.format(data_mount,taskId,sample,replicate), '-o', 'app/static/images/{}/{}/MixMHCpred/{}/{}'.format(taskId,sample,replicate[:-13],replicate), '-a', alleles ])

                    elif(method=='NetMHCpan'):
                        f = open('app/static/images/{}/{}/NetMHCpan/{}/{}'.format(taskId,sample,replicate[:-13],replicate), 'w')
                        p = Popen(['./app/tools/netMHCpan-4.1/netMHCpan', '-p', '{}/{}/{}/{}'.format(data_mount,taskId,sample,replicate), '-a', alleles], stdout=f)
                        output, err = p.communicate(b"input data that is passed to subprocess' stdin")
                        f.close()     

                elif replicate[-7:] == 'mer.txt':
                    if(method=='ANTHEM'):

                        # Here, the peptides which can be accpeted by Anthem will be carry forwarded
                        temp = list()
                        for allele in allelesForAnthem.split(','):
                            
                            if replicate[-8:-7] in ['8','9']:
                                nmer = replicate[-8:-7]
                            else:
                                nmer = replicate[-9:-7]

                            if allele in data[nmer]:
                                temp.append(allele)

                        filteredAllelesForAnthem = ",".join(temp)
                        del temp

                        os.chdir(os.path.join(project_root,'app/tools/Anthem-master'))

                        # Have to find a better way to read output of anthem.
                        # For now reading the files in newly generated folder and deleting them.

                        # Current folder
                        current_files = os.listdir()
                        
                        if filteredAllelesForAnthem != "":
                            call(['../../../lenv/bin/python3','sware_b_main.py', '--HLA', filteredAllelesForAnthem, '--mode', 'prediction', '--peptide_file', '{}/{}/{}/{}'.format(data_mount,taskId,sample,replicate)])

                            present_files = os.listdir()
                            data_folder = list(set(present_files)-set(current_files))

                            # Copying the output file to the destination
                            if replicate[-8:-7] in ['8','9']:
                                copy_tree(data_folder[0], '../../static/images/{}/{}/ANTHEM/{}'.format(taskId,sample,replicate[:-9]))

                            else:
                                copy_tree(data_folder[0], '../../static/images/{}/{}/ANTHEM/{}'.format(taskId,sample,replicate[:-10]))
        
                            # Deleting the output file
                            shutil.rmtree(data_folder[0])

                            os.chdir(project_root)


        if method=='ANTHEM' and sample!='Control':

            os.chdir('{}/app/static/images/{}/{}/ANTHEM/'.format(project_root,taskId,sample))
        
            for replicate in os.listdir('./'):
            #     if replicate[-4:] == '.txt':

                zip_file = ZipFile(replicate+'.zip', 'w')
                for folderName, subfolders, filenames in os.walk(replicate):
                    for filename in filenames:
                        #create complete filepath of file in directory
                        filePath = os.path.join(folderName, filename)
                        # Add file to zip
                        if filename[-3:]=='txt':
                            zip_file.write(filePath, basename(filePath))
                zip_file.close()

            os.chdir(project_root)

def saveBindersData(taskId, alleles, method):

    # Taking out the peptides from the uploaded control data to tag binders present in control group.
    control_peptides = set()
    control_replicates = glob.glob(f'{data_mount}/{taskId}/Control/*8to14mer.txt')

    if len(control_replicates) != 0:
        for control_replicate in control_replicates:
            f = open(control_replicate,'r')
            
            for peptide in f.readlines():
                control_peptides.add(peptide.replace("\n",""))
            
            f.close()

    for sample in os.listdir('{}/{}'.format(data_mount,taskId)):
        for replicate in os.listdir('{}/{}/{}'.format(data_mount,taskId,sample)):
            if sample != 'Control' and replicate[-12:] == '8to14mer.txt':

                # Initialsing the allele and binders collection
                alleles_dict = {}
                
                # ANTHEM case
                if method == 'ANTHEM':
                    # Looping through all the result files of one replicate in ANTHEM results
                    for nmer_file in glob.glob('app/static/images/{}/{}/ANTHEM/{}/*.txt'.format(taskId,sample,replicate[:-13])):
                        # reading nmer file
                        f = open(nmer_file,'r')
                        lines = f.read()
                        f.close()

                        # Using regex taking out alleles and peptides.
                        alleles = list()
                        alleles = re.findall(r'(?P<allele>[A-Z]{3}-[A-Z]\S\d{2}:\d{2}).*[\n] .*[\n]\S*\s*\S*\s*(?P<peptides>(?:.*\.\d*\s)*)',lines)

                        res_cols = ['Peptide', 'Binder', 'Score']


                        # For each allele in one nmer file...
                        for i in range(len(alleles)):

                            # if coming across new allele, adding it into the result dictionary.
                            if alleles[i][0] not in alleles_dict.keys():
                                alleles_dict[alleles[i][0]] = list()


                            # for each peptide under one allele...
                            for peptide in alleles[i][1].split('\n'):

                                # checking if it is a binder
                                if 'yes' in peptide:

                                    # adding into the results list if a binder
                                    alleles_dict[alleles[i][0]].append(peptide.split())

                        
                        
                    # Saving the allele-binders dictionary as text files for ANTHEM case.

                    for allele, binders in alleles_dict.items():
                        
                        alleles_dict[allele] = pd.DataFrame(binders, columns=res_cols)
                        
                        alleles_dict[allele]['Binding Level'] = ""
                        alleles_dict[allele]['Control'] = ""
                        
                        # Keeping strong and weak binders only
                        alleles_dict[allele] = alleles_dict[allele][alleles_dict[allele].apply(lambda x : float(x['Score'])<2,axis=1)]

                        # Tagging each binder as SB(Strong binder) or WB(Weak binder)
                        # alleles_dict[allele]['Binding Level'] = alleles_dict[allele]['Score'].apply(lambda x : 'SB' if float(x)>0.95 else 'WB')
                        alleles_dict[allele]['Binding Level'] = "B"

                        # Tagging binders present in control group
                        temp = alleles_dict[allele]
                        alleles_dict[allele]['Control'] = temp['Peptide'].apply(lambda x : 'Y' if x in control_peptides else '')
                        
                        allele_fromatted = allele[4:].replace("*","").replace(":","")

                        # Sorting the binders from stong to weak binding level and saving it.
                        alleles_dict[allele].sort_values(by=['Score'],ascending=False)[["Peptide","Score","Binding Level","Control"]].to_csv('app/static/images/{}/{}/{}/{}/binders/{}/{}_{}_{}_binders.csv'.format(taskId,sample,method,replicate[:-13],allele_fromatted,replicate[:-13],allele_fromatted,method), index=False)

                # MixMHCpred case
                if method == 'MixMHCpred':
                    f = pd.read_csv('app/static/images/{}/{}/MixMHCpred/{}/{}'.format(taskId,sample,replicate[:-13],replicate),skiprows=11,sep='\t')

                    f['Binding Level'] = ""
                    f['Control'] = ""
                        
                    # Keeping strong and weak binders only
                    f = f[f.apply(lambda x : float(x['%Rank_bestAllele'])<2,axis=1)]

                    # Tagging each binder as SB(Strong binder) or WB(Weak binder)
                    f['Binding Level'] = f['%Rank_bestAllele'].apply(lambda x : 'SB' if float(x)<0.5 else 'WB')
                        
                    # Tagging binders present in control group
                    f['Control'] = f['Peptide'].apply(lambda x : 'Y' if x in control_peptides else '')

                    for allele in alleles.split(','):
                        f[f['BestAllele'] == allele].sort_values(by=['%Rank_bestAllele'])[['Peptide','%Rank_bestAllele','Binding Level','Control']].to_csv('app/static/images/{}/{}/{}/{}/binders/{}/{}_{}_{}_binders.csv'.format(taskId,sample,method,replicate[:-13],allele,replicate[:-13],allele,method), index=False)

                # netMHCpan case
                if method == 'NetMHCpan':
                    f = open('app/static/images/{}/{}/NetMHCpan/{}/{}'.format(taskId,sample,replicate[:-13],replicate),'r')

                    # Storing all entries as string
                    lines = list(f.readlines())
                    f.close() 
                    
                    # Tags in output files of NetMHCpan
                    border = '---------------------------------------------------------------------------------------------------------------------------\n'
                    header = ' Pos         MHC        Peptide      Core Of Gp Gl Ip Il        Icore        Identity  Score_EL %Rank_EL BindLevel\n'

                    allele_dict = {}

                    # Reading data of each allele from output files.
                    for i in range(len(lines)):
                        if lines[i] == border and lines[i+1] == header:
                            allele = lines[i+3].split()[1][4:].replace("*","").replace(":","")
                            allele_dict[allele] = list()
                            
                            i = i+3

                            # Adding all data of one allele till the next allele data starts.
                            while(lines[i]!=border):
                                allele_dict[allele].append(lines[i])
                                i = i+1
                        else:
                            continue
                    
                    # Keeping only required fields
                    for allele, binders in allele_dict.items():
                        temp = list()
                        
                        for binder in binders:
                            temp.append(binder.split()[:13])
                            
                        allele_dict[allele] = temp.copy()
                        
                        del temp

                    # Furhter filtration of the result files.
                    for allele,binders in allele_dict.items():
                        
                        # Keeping only required fields
                        allele_dict[allele] = pd.DataFrame(binders,columns=header.split()[:13])
                        
                        # Making the allele name to same as user input's format
                        allele_dict[allele]['MHC'] = allele
                        
                        temp = allele_dict[allele]
                        # Keeping strong and weak binders only
                        allele_dict[allele] = temp[temp.apply(lambda x : float(x['%Rank_EL'])<2,axis=1)]
                        
                        temp = allele_dict[allele]
                        # Tagging each binder as SB(Strong binder) or WB(Weak binder)
                        allele_dict[allele]['Binding Level'] = temp['%Rank_EL'].apply(lambda x : 'SB' if float(x)<0.5 else 'WB')
                        
                        # Initialsing the column to show if peptide present in control/blank sample.
                        allele_dict[allele]['Control'] = ""

                        
                        # Tagging binders present in control group
                        temp = allele_dict[allele]
                        allele_dict[allele]['Control'] = temp['Peptide'].apply(lambda x : 'Y' if x in control_peptides else '')
                        
                        # Sorting the binders from stong to weak binding level and saving it.
                        allele_dict[allele].sort_values(by=['%Rank_EL'])[["Peptide","%Rank_EL","Binding Level","Control"]].to_csv('app/static/images/{}/{}/{}/{}/binders/{}/{}_{}_{}_binders.csv'.format(taskId,sample,method,replicate[:-13],allele,replicate[:-13],allele,method), index=False)


def getPredictionResuslts(taskId,alleles,methods,samples):

    alleles = alleles.split(',')

    predicted_binders = {}

    for sample in samples:
        if sample !='Control':

            os.chdir(project_root)
            predicted_binders[sample] = {}
            
            for allele in alleles:
                predicted_binders[sample][allele] = {}

                for method in methods:
                    predicted_binders[sample][allele][method] = {}

    for sample in os.listdir(f'{data_mount}/{taskId}/'):
        
        for method in methods:
            
            for replicate in os.listdir(f'{data_mount}/{taskId}/{sample}/'):
                if replicate[-12:] == '8to14mer.txt':

                    for allele in alleles:
                        os.chdir('app/')
                        temp = glob.glob(f'static/images/{taskId}/{sample}/{method}/{replicate[:-13]}/binders/{allele}/*.csv')
                        if len(temp) != 0:
                            predicted_binders[sample][allele][method][replicate[:-13]]= temp[0]
                        os.chdir(project_root)
                        
    return predicted_binders
