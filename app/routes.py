from app import app
from flask import render_template, request, redirect, url_for, flash
from flask import current_app
from app.forms import InitialiserForm, ParentForm
import pandas as pd
import os
import subprocess
from werkzeug.utils import secure_filename
from app.utils import plot_lenght_distribution, filterPeaksFile, saveNmerData, getSeqLogosImages, getGibbsImages, generateBindingPredictions, saveBindersData, getPredictionResuslts
from app.sample import Sample
import time
from pathlib import Path

# Experiment ID
TASK_COUNTER = 0

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

    data_mount = app.config['IMMUNOLYSER_DATA']

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

    # characters_not_dir_name = ['\\','/',':','*','?','"','<','>','|',' ']

    # Creating folders to store images
    for key, value  in request.files.items():
        sample_name = request.form[key]

        # sample_name = sample_name.replace('')

        # for character_not_dir_name in characters_not_dir_name:
        #     sample_name = sample_name.replace(character_not_dir_name,'_') 
            
        # Sk ipping Control Data
        if sample_name == "Control":
            continue

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
        if sample_name != "Control":
            data[sample_name] = list()

        replicates = request.files.getlist(key)
        print("replics name : {}".format(replicates))
        for replicate in replicates:
            file_filename = secure_filename(replicate.filename)

            # If there is no control file uploaded then there is no point to save it.
            if file_filename != "":
                replicate.save(os.path.join(dirName, sample_name, file_filename))

            # Not including the control group in sample data dict
            if sample_name != "Control":
                data[sample_name].append(file_filename)
            elif file_filename !="":
                control.append(file_filename)

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
                            path = os.path.join('app', 'static', 'images', taskId, sample, predictionTool, replicate[:-4], 'binders',allele)
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
    control_data = {}

    # Loading sample data in pandas frames
    for sample_name, file_names in data.items():
        sample_data[sample_name] = dict()
        for replicate in file_names:
            sample_data[sample_name][replicate] = pd.read_csv(os.path.join(dirName, sample_name, replicate))

    # Loading control data in pandas frames
    for control_replicate in control:
        control_data[control_replicate] = pd.read_csv(os.path.join(dirName, "Control", control_replicate))

    # Have to later add the user input for length
    for sample_name, sample in sample_data.items():
        sample_data[sample_name] = filterPeaksFile(sample, minLen=minLen, maxLen=maxLen, control_data= control_data)


    bar_percent = plot_lenght_distribution(sample_data, hist='percent')
    bar_density = plot_lenght_distribution(sample_data, hist='density')

    # Saving 8 to 14 nmers for mhc1 predictions
    saveNmerData(dirName, sample_data, peptideLength=[8,14])

    for i in range(8,14):
        saveNmerData(dirName, sample_data, peptideLength=i)


    # Calling script to generate sequence logos
    subprocess.call('sudo python3 {} {} {}'.format(os.path.join('app','seqlogo.py'), taskId, data_mount), shell=True)

    # Method to return names of png files of seqlogos
    # This value is supposed to be returned from saveNmerDate method but for now writting
    # temporary script to return names of seqlogos pngs files in a dictionary.
    
    seqlogos = getSeqLogosImages(sample_data)

    # Calling script to generate gibbsclusters
    subprocess.call('sudo python3 {} {} {}'.format(os.path.join('app', 'gibbscluster.py'), taskId, data_mount), shell=True)

    # Getting names of the gibbscluster
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
    

    return render_template('analytics.html', taskId=taskId, peptide_percent=bar_percent, peptide_density=bar_density, seqlogos = seqlogos, gibbsImages = gibbsImages, analytics=True,predicted_binders=predicted_binders, predictionTools = predictionTools)

@app.route("/analytics")
def analytics():

    predicted_binders = getPredictionResuslts('202103092206303','B1301,C0303,A0202',['NetMHCpan', 'ANTHEM', 'MixMHCpred'],['VMM1','Inferon','Tmpi'])

    return render_template("temp.html", analytics=True,predicted_binders=predicted_binders, predictionTools = ['NetMHCpan', 'ANTHEM', 'MixMHCpred'])

# Method to manage experiment ID
def getTaskId():
    global TASK_COUNTER
    TASK_COUNTER = TASK_COUNTER+1
    task_Id = time.strftime("%Y%m%d%H%M%S")+str(TASK_COUNTER)

    return task_Id
