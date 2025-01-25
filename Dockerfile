FROM python:3.12-slim

# Set the working directory in the container
WORKDIR /app

# Install build tools and other required dependencies
RUN apt-get update && apt-get install -y \
    git \
    tar \
    build-essential \
    wget \
    libssl-dev \
    libbz2-dev \
    libreadline-dev \
    libsqlite3-dev \
    zlib1g-dev && \
    wget https://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs926/ghostscript-9.26.tar.gz && \
    tar -xzf ghostscript-9.26.tar.gz && \
    cd ghostscript-9.26 && \
    ./configure && make && make install && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Clone the repository
RUN git clone https://github.com/prmunday/Immunolyser /app/Immunolyser

# Change to the repository directory
WORKDIR /app/Immunolyser

# Checkout the develop branch
RUN git fetch --all && git checkout feature/docker_setup

# Make images under static folder
RUN mkdir -p /app/Immunolyser/app/static/images

# Copy the seq2logo tar.gz file from the local tools folder to the container
COPY /tools/seq2logo-2.1.all.tar.gz /app/Immunolyser/app/tools/

# Create a tools directory and extract the tar.gz file there
RUN mkdir -p /app/Immunolyser/app/tools && \
    tar -xzf /app/Immunolyser/app/tools/seq2logo-2.1.all.tar.gz -C /app/Immunolyser/app/tools && \
    rm /app/Immunolyser/app/tools/seq2logo-2.1.all.tar.gz

# Create a virtual environment for Python 3
RUN python3 -m venv lenv

# Install dependencies for Python 2 and Python 3
RUN /bin/bash -c "source lenv/bin/activate && \
    pip install -r requirements_python2.txt && \
    pip install -r requirements_python3.txt"

# Install Celery
RUN /bin/bash -c "pip install celery"

# Run the hotfix script
RUN /bin/bash -c "python hotfix_package_files.py"

# Download and extract Python 2.7.18
RUN wget https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz && \
    tar xvf Python-2.7.18.tgz

# Build and install Python 2.7.18
WORKDIR Python-2.7.18
RUN ./configure && \
    make && \
    make install

# Install pip for Python 2.7
RUN wget https://bootstrap.pypa.io/pip/2.7/get-pip.py && \
    python2 get-pip.py

# Install numpy for Python 2.7
RUN python2 -m pip install numpy

# Change to the repository directory
WORKDIR /app/Immunolyser

# Expose Flask and Celery ports
EXPOSE 5000
EXPOSE 5555

# Set a default environment variable for IMMUNOLYSER_DATA
ENV IMMUNOLYSER_DATA=/app/data

# Copy entrypoint script
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

# Set the entrypoint script
ENTRYPOINT ["/entrypoint.sh"]
