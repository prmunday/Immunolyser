from pandas.tseries.offsets import Second
from app import app
from flask import render_template, request, jsonify
from app import sample
from app.forms import InitialiserForm, ParentForm
import pandas as pd
import os
import subprocess
from werkzeug.utils import secure_filename
from app.utils import plot_lenght_distribution, filterPeaksFile, saveNmerData, getSeqLogosImages, getGibbsImages, generateBindingPredictions, saveBindersData, getPredictionResuslts, getPredictionResusltsForUpset, findNumberOfPeptidesInCore
from app.sample import Sample
import time
from pathlib import Path
import glob
import shutil
from app.Pepscan import PepScan
import re
from collections import Counter,OrderedDict

project_root = os.path.dirname(os.path.realpath(os.path.join(__file__, "..")))

# Experiment ID
TASK_COUNTER = 0

# DEMO Task ID
DEMO_TASK_ID = "202208221407254"

data_mount = app.config['IMMUNOLYSER_DATA']

@app.route("/initialiser", methods=["POST", "GET"])
def initialiser():    
    samples = []

    # Have to take this input from user
    maxLen = 30
    minLen = 5

    taskId = getTaskId()
    dirName = os.path.join(data_mount, taskId)
    try:
        # Create target Directory
        os.makedirs(dirName)
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")


    data = {}
    control = list()

    # Creating folders to store images
    for key, value  in request.files.items():
        sample_name = request.form[key]

        # Skipping Control Data
        # if sample_name == "Control":
        #     continue

        # Creating sub directories to store sample data
        try:
            # for seqlogos
            path_for_logos = os.path.join('app', 'static', 'images', taskId, sample_name, 'seqlogos')
            if not os.path.exists(path_for_logos):
                # os.makedirs(directory)
                Path(path_for_logos).mkdir(parents=True, exist_ok=True)
                print("Directory Created") 
        except FileExistsError:
            print("Directory already exists")    

    # Saving the data and loading into the dictionary
    for key, value  in request.files.items():
        sample_name = request.form[key]

        # Creating sub directories to store sample data
        try:
            os.mkdir(os.path.join(dirName, sample_name))
            print("Directory Created") 
        except FileExistsError:
            print("Directory already exists")    

        # Not including the control group in data dict 
        # if sample_name != "Control":
        data[sample_name] = list()

        replicates = request.files.getlist(key)
        print("replics name : {}".format(replicates))
        for replicate in replicates:
            file_filename = secure_filename(replicate.filename)

            # If there is no control file uploaded then there is no point to save it.
            if file_filename != "":
                replicate.save(os.path.join(dirName, sample_name, file_filename))

            # Not including the control group in sample data dict
            # if sample_name != "Control":
                data[sample_name].append(file_filename)
            # elif file_filename !="":
                # control.append(file_filename)

        # If control data is not uploaded, then deleting the sample from the data dictionary
        temp = data.copy()
        for sample_name, replicates in temp.items():
            if len(replicates) == 0:
                data.pop(sample_name)

    if len(data) == 0:
        return render_template("initialiser.html", initialiser=True)

#     # Storing in classes
#     # sample1 = Sample(sample_one_name)
#     # sample2 = Sample(sample_two_name)

#     # Or in variables
#     sample1 = {}
#     sample2 = {}

    alleles_unformatted = request.form.get('alleles')

    # saving alleles selected in a file
    allele_file = open(os.path.join('app', 'static', 'images', taskId, "selectedalleles.txt"), "w")
    allele_file.write(alleles_unformatted)
    allele_file.close()

    # Prediction tools selected by the user
    predictionTools = request.form.getlist('predictionTools')
    print("Prediction tools selected: {}".format(predictionTools))

    # Creating directories to store binding prediction results
    for sample, replicates in data.items():
        for predictionTool in predictionTools:
            for replicate in replicates:
                if alleles_unformatted != "":
                    for allele in alleles_unformatted.split(','):
                        try:
                            if sample != 'Control':
                                path = os.path.join('app', 'static', 'images', taskId, sample, predictionTool, replicate[:-4], 'binders',allele)
                            else:
                                path = os.path.join('app', 'static', 'images', taskId, sample)
                            if not os.path.exists(path):
                                # os.makedirs(directory)
                                Path(path).mkdir(parents=True, exist_ok=True)
                                print("Directory Created : {}".format(path))
                        except FileExistsError:
                            print("Directory already exists {}".format(path))

    # Converting alleles from A0203 format to HLA-A02:03
    alleles = ''
    if alleles_unformatted != "":
        temp = list()
        for allele in alleles_unformatted.split(','):
            temp.append('HLA-{}:{}'.format(allele[:3],allele[3:]))
            
        alleles = ",".join(temp)
        del temp
                
    sample_data = {}
    # control_data = {}

    # Loading sample data in pandas frames
    for sample_name, file_names in data.items():
        sample_data[sample_name] = dict()
        for replicate in file_names:
            sample_data[sample_name][replicate] = pd.read_csv(os.path.join(dirName, sample_name, replicate))

    # Loading control data in pandas frames
    # for control_replicate in control:
        # control_data[control_replicate] = pd.read_csv(os.path.join(dirName, "Control", control_replicate))

    # Have to later add the user input for length
    for sample_name, sample in sample_data.items():
        sample_data[sample_name] = filterPeaksFile(sample, minLen=minLen, maxLen=maxLen)


    bar_percent = plot_lenght_distribution(sample_data, hist='percent')
    bar_density = plot_lenght_distribution(sample_data, hist='density')

    # Saving 8 to 14 nmers for mhc1 predictions
    saveNmerData(dirName, sample_data, peptideLength=[8,14], unique = True)

    for i in range(8,14):
        saveNmerData(dirName, sample_data, peptideLength=i, unique = True)

    # Saving 9mer of each file including duplicates. This will be used for the peptide distribution graph
    saveNmerData(dirName, sample_data, peptideLength=9, unique=False)

    # # Calling script to generate sequence logos
    subprocess.call('sudo python3 {} {} {}'.format(os.path.join('app','seqlogo.py'), taskId, data_mount), shell=True)

    # # Method to return names of png files of seqlogos
    # # This value is supposed to be returned from saveNmerDate method but for now writting
    # # temporary script to return names of seqlogos pngs files in a dictionary.
    
    seqlogos = getSeqLogosImages(sample_data)
    # seqlogos = {}
    
    # gibbsImages = {}
    
    # # Calling script to generate gibbsclusters
    subprocess.call('sudo python3 {} {} {}'.format(os.path.join('app', 'gibbscluster.py'), taskId, data_mount), shell=True)

    # # Getting names of the gibbscluster
    gibbsImages = getGibbsImages(taskId, sample_data)

    # Generating binding predictions
    if alleles!="":    
        for predictionTool in predictionTools:
            generateBindingPredictions(taskId, alleles, predictionTool)

    # Fetching the binders from the results
    if alleles!="":    
        for predictionTool in predictionTools:
            saveBindersData(taskId, alleles_unformatted, predictionTool)


    predicted_binders = None
    # Method to get the prediction results
    if alleles!="":    
        predicted_binders = getPredictionResuslts(taskId,alleles_unformatted,predictionTools,sample_data.keys())
    
    upsetLayout = getPredictionResusltsForUpset(taskId,alleles_unformatted,predictionTools,sample_data.keys())

    return render_template('analytics.html', taskId=taskId, peptide_percent=bar_percent, peptide_density=bar_density, seqlogos = seqlogos, gibbsImages = gibbsImages, analytics=True,predicted_binders=predicted_binders, predictionTools = predictionTools,upsetLayout=upsetLayout)

@app.route("/analytics")
def analytics():

    return render_template("error.html",analytics=True, msg = 'initialiser')

@app.route('/<taskId>')
def getExistingReport(taskId):

    global DEMO_TASK_ID
    demo = False
    # Static ID for the demo
    if str(taskId) == DEMO_TASK_ID:
        demo = True
        pass
    elif str(taskId).isnumeric() == False:
        return 'No task id recieved to generate old report'

    predictionTools = ['MixMHCpred','ANTHEM','NetMHCpan']
    
    data = {}
    maxLen = 30
    minLen = 5
    sample_data = {}
    dirName = os.path.join(data_mount, taskId)
    predicted_binders = None
    with open(os.path.join('app', 'static', 'images', taskId, "selectedalleles.txt")) as f:
        alleles_unformatted = f.readline() 

    samples =[ f.name for f in os.scandir(dirName) if f. is_dir()]

    # Saving the data and loading into the dictionary
    for sample_name  in samples:

        # Not including the control group in data dict 
        # if sample_name != "Control":
        data[sample_name] = list()

        filenames = os.listdir(os.path.join(dirName,sample_name))
        replicates = [ filename for filename in filenames if filename.endswith( ".csv" ) ]

        for file_filename in replicates:
            data[sample_name].append(file_filename)
            
        # If control data is not uploaded, then deleting the sample from the data dictionary
        temp = data.copy()
        for sample_name, replicates in temp.items():
            if len(replicates) == 0:
                data.pop(sample_name)


    # Loading sample data in pandas frames
    for sample_name, file_names in data.items():
        sample_data[sample_name] = dict()
        for replicate in file_names:
            sample_data[sample_name][replicate] = pd.read_csv(os.path.join(dirName, sample_name, replicate))

    # Loading control data in pandas frames
    # for control_replicate in control:
        # control_data[control_replicate] = pd.read_csv(os.path.join(dirName, "Control", control_replicate))

    # Have to later add the user input for length
    for sample_name, sample in sample_data.items():
        sample_data[sample_name] = filterPeaksFile(sample, minLen=minLen, maxLen=maxLen)


    bar_percent = plot_lenght_distribution(sample_data, hist='percent')
    bar_density = plot_lenght_distribution(sample_data, hist='density')

    seqlogos = getSeqLogosImages(sample_data)
    gibbsImages = getGibbsImages(taskId, sample_data)

    if alleles_unformatted != '':
        predicted_binders = getPredictionResuslts(taskId,alleles_unformatted,predictionTools,sample_data.keys())

    upsetLayout = getPredictionResusltsForUpset(taskId,alleles_unformatted,predictionTools,sample_data.keys())

    return render_template('analytics.html', taskId=taskId, analytics=True, demo=demo, peptide_percent=bar_percent, peptide_density=bar_density, seqlogos =seqlogos, gibbsImages=gibbsImages, upsetLayout=upsetLayout, predicted_binders=predicted_binders,predictionTools=predictionTools)

@app.route("/feedback", methods=["POST", "GET"])
def feedback():

    if 'feedback' not in request.form.keys():
        return render_template("feedback.html", feedback=None)

    if len(request.form.get('feedback'))==0:
        return render_template("feedback.html", feedback=None)
        

    data_mount = app.config['IMMUNOLYSER_DATA']

    dirName = os.path.join(data_mount,'feedback')
    try:
        # Create target Directory
        os.makedirs(dirName)
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")

    feedback_file = open(os.path.join(dirName,'feedback.txt'), 'a')

    feedback_file.write(request.form.get('feedback')+'\n\n---Feedback---\n\n')
    # Close the file
    feedback_file.close()

   
    
    return render_template("feedback.html", feedback=True)

# Method to manage experiment ID
def getTaskId():
    global TASK_COUNTER
    TASK_COUNTER = TASK_COUNTER+1
    task_Id = time.strftime("%Y%m%d%H%M%S")+str(TASK_COUNTER)

    return task_Id

# This method is to create the bar graphs for an input file not created already
@app.route("/api/generateGibbs", methods=["POST"])
def createGibbsBar():
    
    cluster = request.form['cluster']
    taskId = request.form['taskId']
    replicate = request.form['replicate']
    sample = request.form['sample']

    barLocation = glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/images/gibbs.KLDvsCluster.barplot.JPG')

    if len(barLocation) ==1:
        barLocation = barLocation[0][4:]

    else:

        # Deleting previous meta files present
        dirpath = f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}'
        if os.path.exists(dirpath) and os.path.isdir(dirpath):
            shutil.rmtree(dirpath)

        subprocess.call('sudo python3 {} {} {} {} {}'.format(os.path.join('app', 'gibbsclusterBarGraph.py'), taskId, data_mount, sample, replicate), shell=True)

        barLocation = glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/images/gibbs.KLDvsCluster.barplot.JPG')

        if len(barLocation) ==1:
            barLocation = barLocation[0][4:]
        
        else: 
            barLocation = f'/static/others/gibbsBarNotFound.JPG'
        
    if len(cluster) == 0:
        bestCluster = pd.read_table(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/images/gibbs.KLDvsClusters.tab')
        bestCluster = bestCluster[bestCluster.columns].sum(axis=1).idxmax()

        seqClusters = [[x[4:], "Could not be calculated"] for x in glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/logos/gibbs_logos_*of{bestCluster}*.jpg')]

    else:
        seqClusters = [x[4:] for x in glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/logos/gibbs_logos_*of{cluster}*.jpg')]

        if len(seqClusters) != int(cluster):
            subprocess.call('sudo python3 {} {} {} {} {} {}'.format(os.path.join('app', 'gibbsclusterSeqLogo.py'), taskId, data_mount, sample, replicate, int(cluster)), shell=True)
            seqClusters = [[x[4:], "Could not be calculated"] for x in glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/logos/gibbs_logos_*of{cluster}*.jpg')]

    # Finding the number of records used for the cluster
    findNumberOfPeptidesInCore(seqClusters, taskId, sample, replicate)

    return {barLocation:seqClusters}

# This method is used to get the found binders in different combinations
@app.route("/api/getBinders", methods=["POST"])
def getBinders():

    tool = request.form['tool']
    taskId = request.form['taskId']
    allele = request.form['allele']
    listonly = request.form['list']
    replicates = request.form['replicates']

    predictionTools = ['ANTHEM','NetMHCpan','MixMHCpred']

    samples = {}


    replicates = replicates.split(',')
    
    for i in replicates:
        
        entries = i.split(';')

        if entries[0] == allele:

            if entries[1] not in samples.keys():
                samples[entries[1]] = []
                samples[entries[1]].append(entries[2])
            
            else:
                samples[entries[1]].append(entries[2])


    if listonly == "":

        return 'bindersFile'
    
    else:
        res = []
        for sample,replicates in samples.items():
            if tool !="":
                binder_files = []
            else:
                binder_files = {}
                for i in predictionTools:
                    binder_files[i] = []

            for replicate in replicates:
                
                if tool == "":

                    binding = getPredictionResuslts(alleles=allele,taskId=taskId,methods=predictionTools,samples=[sample])

                    for method in predictionTools:

                        try:
                            binder_files[method].append(binding[sample][allele][method][replicate])
                        except KeyError:
                            continue

                else:
                    binding = getPredictionResuslts(alleles=allele,taskId=taskId,methods=[tool],samples=[sample])
                    try:
                        binder_files.append(binding[sample][allele][tool][replicate])
                    except KeyError:
                        continue

            binders = []

            if tool == "":
                binders = {}

                for i,j in binder_files.items():
                    binders[i] = []

                    for k in j:
                        binders[i].extend(pd.read_csv(os.path.join('app',k))['Peptide'].to_list())    

                first = set(binders[predictionTools[0]]).intersection(set(binders[predictionTools[1]]))
                second = set(binders[predictionTools[1]]).intersection(set(binders[predictionTools[2]]))
                third = set(binders[predictionTools[0]]).intersection(set(binders[predictionTools[2]]))

                binders = first.union(second).union(third)

            else:
                for i in binder_files:
                    binders.extend(pd.read_csv(os.path.join('app',i))['Peptide'].to_list())
                binders = set(binders)

            res.append({'name':sample,'elems':list(binders)})
    
    return jsonify(res)

@app.route("/api/getSeqLogo", methods=["POST"])

def getSeqLogo():

    name = request.form['name']

    # removing '(' and ')' from name if present (causing problem for subprocess call)
    print(name)
    name = str(name).replace('∩','and').replace('(','').replace(')','').replace(' ','_').strip()


    taskId = request.form['taskId']    
    elems = request.form['elems']

    # print(type(elems))
    # peptides = json.loads(elems)
    peptides = elems.split(',')
    peptides_location_forseqlogo = os.path.join(project_root,'app','static','images',taskId,'selected-9mer-binders-for-seqlogo.txt')
    binders_location = os.path.join(project_root,'app','static','images',taskId,'selected-binders.txt')

    peptides = pd.DataFrame(peptides)
    peptides.columns = ['peptide']

    total_peptides = peptides.shape[0]
    # Saving all binders which can be downloaded in text files
    peptides.to_csv(binders_location,index=False,header=False)

    peptides =  peptides[peptides.peptide.apply(lambda x: True if len(x)==9 else False)]

    nine_mers = peptides.shape[0]

    # Saving 9mer binders to be used for seq2logo generation
    peptides.to_csv(peptides_location_forseqlogo,index=False,header=False)

    seqLogoLocation = os.path.join(project_root,'app','static','images',taskId,'seqLogoApi')

    print('sudo python3 {} {} {} {} {} {}'.format(os.path.join('app','seqlogoAPI.py'),peptides_location_forseqlogo,seqLogoLocation,name,nine_mers,total_peptides))

    subprocess.call('sudo python3 {} {} {} {} {} {}'.format(os.path.join('app','seqlogoAPI.py'),peptides_location_forseqlogo,seqLogoLocation,name,nine_mers,total_peptides), shell=True)

    return os.path.join('static','images',taskId,'seqLogoApi-001.jpg')

@app.route("/help")
def help():
    
    # Checking to platform, if it is windows, wsl will be initialised
    #if request.user_agent.platform =='windows':
        #setUpWsl()
    return render_template("help.html", help=True)

@app.route("/pepscanner", methods=["GET"])
def pepscanner():

    return render_template("pepscanner.html", pepscanner=True,pep=False)

@app.route("/api/pepscanner", methods=["POST"])
def generatePepscanner():

    taskId = getTaskId()

    
    dirName = os.path.join('app', 'static', 'images', taskId,'protienandepeptides')
    try:
        # Create target Directory
        os.makedirs(dirName)
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")


    # Extracing passed file and the peptides list
    uploaded_file = request.files['file']
    peptides = request.form['peptides']    
    fileName = uploaded_file.filename.replace('C:\\fakepath\\',"")

    # Deleing any existing heatmap
    myfile=os.path.join(project_root,'app/static/images/pepscanner.png')

    ## If file exists, delete it ##
    if os.path.isfile(myfile):
        os.remove(myfile)
    else:    ## Show an error ##
        print("Error: %s file not found" % myfile)


    # Attaching reference proteome file
    ref_proteome = os.path.join(project_root,'app','references data','uniprot-proteome_UP000005640.fasta')

    # Ignore this file for now
    accessionsidfile = os.path.join(project_root,'app','references data','accessionids.csv')

    # Input peptide file
    peptides_file = os.path.join(project_root,'app','static','images',taskId,fileName)
    
    if uploaded_file.filename != '':
        uploaded_file.save(peptides_file)

    scanner = PepScan()

    inputFile = pd.read_csv(peptides_file)
    metadata = findMostOccuringAccessionIds(inputFile, taskId, fileName)

    scanner.search_proteome(peptide_file=peptides_file, proteome_file=ref_proteome,accessionsids=accessionsidfile)

    if peptides != 'Type proteins of interest (optional; human proteins only)':
        # Getting the list of peptides entered (and preprocessing, e.g., removing empty strings)
        peptides = peptides.replace(' ','').split(',')
        while("" in peptides) :
            peptides.remove("")
    else:
        peptides = list(metadata['top_protiens'].keys())

        if len(peptides) > 5:
            peptides = peptides[:5]

    print('Peptides passed for pepscanner: {}'.format(peptides))

    scanner.peptide_dist(peptides, taskId)

    metadata['taskId'] = taskId

    return jsonify(metadata)

def findMostOccuringAccessionIds(inputFile, taskId, inputFileName):
    
    accessions = inputFile['Accession'].to_list()
    accessionIds = []
    metadata = {}
    
    for i in accessions:
        for j in str(i).split(':'):
            found = re.search(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",j)

            if found:
                accessionIds.append(found[0])
    

    metadata['unique_peptides'] = len(set(accessionIds))
    metadata['top_protiens'] = {}

    if metadata['unique_peptides'] >= 10:
        accessionIds = Counter(accessionIds).most_common(10)
    else:
        accessionIds = Counter(accessionIds).most_common(metadata['unique_peptides'])

    # Reading mapping file
    mapping = pd.read_csv(os.path.join(project_root,'app','references data','proteinmapping.csv'))
        
    odict = OrderedDict(accessionIds)

    for key, value in odict.items():
        
        if (len (mapping[mapping['Protein'] == key])>0):
            gn = mapping.loc[mapping['Protein'] == key, "GN"].iloc[0]
            species = mapping.loc[mapping['Protein'] == key, "Species"].iloc[0]
        else:
            gn = ""
            species = ""

        odict[key]= [value, gn, species]

    metadata['top_protiens'] = odict

    for accessiondId in metadata['top_protiens'].keys():
        fileName = inputFileName+ ' ' + accessiondId + '.csv'
        subFile = inputFile[inputFile['Accession'].str.contains(accessiondId, na=False)]
        subFile.to_csv(os.path.join(project_root,'app','static','images',taskId,'protienandepeptides',fileName))

    return metadata

# Following method return the pre-run job using the specified task id.
@app.route("/demo")
def demo():
    global DEMO_TASK_ID
    return getExistingReport(DEMO_TASK_ID)