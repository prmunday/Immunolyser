from app import app
from flask import render_template, request, redirect, url_for
from app.forms import InitialiserForm
import pandas as pd
import os
from werkzeug.utils import secure_filename
from app.utils import plot_lenght_distribution, filterPeaksFile
from app.sample import Sample
import time

# Experiment ID
EXPERIMENT_ID = 0

@app.route("/")
@app.route("/index")
@app.route("/home")
def index():
    return render_template("index.html", index=True )

@app.route("/initialiser", methods=["POST", "GET"])
def initialiser():
    form = InitialiserForm()
    if form.validate_on_submit():
        
        experiment_id = getExperminetId()
        
        sample_one_name = form.sample_one_name.data
        sample_two_name = form.sample_two_name.data
        # Have to take this input from user
        maxLen = 30
        minLen = 1

        # Creating directory for the experiment
        # Create directory
        dirName = os.path.join('data', str(experiment_id)+'-ID-'+time.strftime("%Y-%m-%d"))
        try:
            # Create target Directory
            os.makedirs(dirName)
            print("Directory " , dirName ,  " Created ") 
        except FileExistsError:
            print("Directory " , dirName ,  " already exists")

        # Creating sub directories to store sample data
        try:
            os.mkdir(os.path.join(dirName, sample_one_name))
            print("Directory Created") 
        except FileExistsError:
            print("Directory already exists")

        try:
            os.mkdir(os.path.join(dirName, sample_two_name))
            print("Directory Created") 
        except FileExistsError:
            print("Directory already exists")        

        sample_one = []
        for sample in form.sample_one.data:
            file_filename = secure_filename(sample.filename)
            sample.save(os.path.join(dirName, sample_one_name, file_filename))
            sample_one.append(file_filename)

        sample_two = []
        for sample in form.sample_two.data:
            file_filename = secure_filename(sample.filename)
            sample.save(os.path.join(dirName, sample_two_name, file_filename))
            sample_two.append(file_filename)

        
        # Storing in classes
        # sample1 = Sample(sample_one_name)
        # sample2 = Sample(sample_two_name)

        # Or in variables
        sample1 = {}
        sample2 = {}

        for replicate in sample_one:
            sample1[replicate] = pd.read_csv(os.path.join(dirName, sample_one_name, replicate))

        for replicate in sample_two:
            sample2[replicate] = pd.read_csv(os.path.join(dirName, sample_two_name, replicate))   

        # Have to later add the user input for length
        sample1 = filterPeaksFile(sample1, minLen=minLen, maxLen=maxLen)
        sample2 = filterPeaksFile(sample2, minLen=minLen, maxLen=maxLen)


        bar = plot_lenght_distribution({sample_one_name:sample1, sample_two_name:sample2})

        return render_template('analytics.html', plot=bar, analytics=True)
    return render_template("initialiser.html", form=form, initialiser=True)

@app.route("/analytics")
def analytics():
    return render_template("analytics.html", analytics=True)

# Method to manage experiment ID
def getExperminetId():
    global EXPERIMENT_ID
    EXPERIMENT_ID = EXPERIMENT_ID+1

    return EXPERIMENT_ID