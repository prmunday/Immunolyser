from app import app
from flask import render_template, request, redirect, url_for, flash
from app.forms import InitialiserForm, ParentForm
import pandas as pd
import os
import subprocess
from werkzeug.utils import secure_filename
from app.utils import plot_lenght_distribution, filterPeaksFile, saveNmerData, getSeqLogosImages
from app.sample import Sample
import time
from pathlib import Path

# Experiment ID
EXPERIMENT_ID = 0

@app.route("/")
@app.route("/index")
@app.route("/home")
def index():
    temp = {'a':[1,2,3], 'b':[4,4,6]}
    return render_template("index.html", index=True, temp =temp )

@app.route("/initialiser", methods=["POST", "GET"])
def initialiser():
    experiment_id = getExperminetId()
    
    samples = []

    # Have to take this input from user
    maxLen = 30
    minLen = 1

    # Creating directory for the experiment
    # Create directory

    experiment_id_date = str(experiment_id)+'-ID-'+time.strftime("%Y-%m-%d")
    dirName = os.path.join('data', experiment_id_date)
    try:
        # Create target Directory
        os.makedirs(dirName)
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")


    data = {}

    # Creating folders to store images
    for key, value  in request.files.items():
        sample_name = request.form[key]
    
        # Creating sub directories to store sample data
        try:
            # for seqlogos
            path_for_logos = os.path.join('app', 'static', 'images', experiment_id_date, sample_name, 'seqlogos')
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

        data[sample_name] = list()
        replicates = request.files.getlist(key)
        print("replics name : {}".format(replicates))
        for replicate in replicates:
            file_filename = secure_filename(replicate.filename)
            replicate.save(os.path.join(dirName, sample_name, file_filename))
            data[sample_name].append(file_filename)

    if len(data) == 0:
        return render_template("initialiser.html", initialiser=True)

#     # Storing in classes
#     # sample1 = Sample(sample_one_name)
#     # sample2 = Sample(sample_two_name)

#     # Or in variables
#     sample1 = {}
#     sample2 = {}

    sample_data = {}

    for sample_name, file_names in data.items():
        sample_data[sample_name] = dict()
        for replicate in file_names:
            sample_data[sample_name][replicate] = pd.read_csv(os.path.join(dirName, sample_name, replicate))

    # Have to later add the user input for length
    for sample_name, sample in sample_data.items():
        sample_data[sample_name] = filterPeaksFile(sample, minLen=minLen, maxLen=maxLen)


    bar_percent = plot_lenght_distribution(sample_data, hist='percent')
    bar_density = plot_lenght_distribution(sample_data, hist='density')

    # Saving 9 mers data for sequence logos
    saveNmerData(dirName, sample_data, peptideLength=9)

    # Calling script to generate sequence logos
    subprocess.call('python app\\seqlogo.py {}'.format(experiment_id_date), shell=True)

    # Method to return names of png files of seqlogos
    # This value is supposed to be returned from saveNmerDate method but for now writting
    # temporary script to return names of seqlogos pngs files in a dictionary.
    
    seqlogos = getSeqLogosImages(sample_data)

    return render_template('analytics.html', experiment_id_date=experiment_id_date,peptide_percent=bar_percent, peptide_density=bar_density, seqlogos = seqlogos, analytics=True)
    # return render_template("initialiser.html", form=form, initialiser=True)

@app.route("/analytics")
def analytics():
    return render_template("analytics.html", analytics=True, iframe = 'data/6-ID-2020-11-13/Inferon/seqlogos/peptide_AB190613_1106_IFN_PEAKS10_DT9-001.png')

# Method to manage experiment ID
def getExperminetId():
    global EXPERIMENT_ID
    EXPERIMENT_ID = EXPERIMENT_ID+1

    return EXPERIMENT_ID