from app import app
from flask import render_template, request, send_file
from flask import current_app
from app import sample
from app.forms import InitialiserForm, ParentForm
import pandas as pd
import os
import subprocess
from werkzeug.utils import secure_filename
from app.utils import plot_lenght_distribution, filterPeaksFile, saveNmerData, getSeqLogosImages, getGibbsImages, generateBindingPredictions, saveBindersData, getPredictionResuslts
from app.sample import Sample
import time
from pathlib import Path
import glob
import shutil

project_root = os.path.dirname(os.path.realpath(os.path.join(__file__, "..")))

# Experiment ID
TASK_COUNTER = 0

data_mount = app.config['IMMUNOLYSER_DATA']

@app.route("/")
@app.route("/index")
@app.route("/home")
def index():
    
    # Checking to platform, if it is windows, wsl will be initialised
    #if request.user_agent.platform =='windows':
        #setUpWsl()
    return render_template("index.html", index=True)

@app.route("/initialiser", methods=["POST", "GET"])
def initialiser():    
    samples = []

    # Have to take this input from user
    maxLen = 30
    minLen = 1

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
    saveNmerData(dirName, sample_data, peptideLength=[8,14])

    for i in range(8,14):
        saveNmerData(dirName, sample_data, peptideLength=i)


    # # Calling script to generate sequence logos
    # subprocess.call('sudo python3 {} {} {}'.format(os.path.join('app','seqlogo.py'), taskId, data_mount), shell=True)

    # # Method to return names of png files of seqlogos
    # # This value is supposed to be returned from saveNmerDate method but for now writting
    # # temporary script to return names of seqlogos pngs files in a dictionary.
    
    # seqlogos = getSeqLogosImages(sample_data)
    seqlogos = {}
    
    # gibbsImages = {}
    
    # # Calling script to generate gibbsclusters
    subprocess.call('sudo python3 {} {} {}'.format(os.path.join('app', 'gibbscluster.py'), taskId, data_mount), shell=True)

    # # Getting names of the gibbscluster
    gibbsImages = getGibbsImages(taskId, sample_data)
# 
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
    

    return render_template('analytics.html', taskId=taskId, peptide_percent=bar_percent, peptide_density=bar_density, seqlogos = seqlogos, gibbsImages = gibbsImages, analytics=True,predicted_binders=predicted_binders, predictionTools = predictionTools)

@app.route("/analytics")
def analytics():

    return render_template("error.html",analytics=True, msg = 'initialiser')


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

        seqClusters = [x[4:] for x in glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/logos/gibbs_logos_*of{bestCluster}*.jpg')]

    else:
        seqClusters = [x[4:] for x in glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/logos/gibbs_logos_*of{cluster}*.jpg')]

        if len(seqClusters) != int(cluster):
            subprocess.call('sudo python3 {} {} {} {} {} {}'.format(os.path.join('app', 'gibbsclusterSeqLogo.py'), taskId, data_mount, sample, replicate, int(cluster)), shell=True)
            seqClusters = [x[4:] for x in glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/logos/gibbs_logos_*of{cluster}*.jpg')]

    return {barLocation:seqClusters}