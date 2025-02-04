FROM python:3.12-slim

# Set the working directory in the container
WORKDIR /app

# Install build tools, R, tcsh, and other required dependencies
RUN apt-get update && apt-get install -y \
    git \
    tar \
    build-essential \
    wget \
    libssl-dev \
    libbz2-dev \
    libreadline-dev \
    libsqlite3-dev \
    zlib1g-dev \
    r-base \
    man-db \
    ncompress \
    tcsh \
    && \
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

# Copy the seq2logo tar.gz file from the local tools folder to the container
COPY /tools/seq2logo-2.1.all.tar.gz /app/Immunolyser/app/tools/

# Create a tools directory and extract the tar.gz file there
RUN mkdir -p /app/Immunolyser/app/tools && \
    tar -xzf /app/Immunolyser/app/tools/seq2logo-2.1.all.tar.gz -C /app/Immunolyser/app/tools && \
    rm /app/Immunolyser/app/tools/seq2logo-2.1.all.tar.gz

# Copy the gibbscluster tar.gz file to the container
COPY /tools/gibbscluster-2.0f.Linux.tar.gz /app/Immunolyser/app/tools/

# Extract gibbscluster tar.gz
RUN mkdir -p /app/Immunolyser/app/tools && \
    tar -xvf /app/Immunolyser/app/tools/gibbscluster-2.0f.Linux.tar.gz -C /app/Immunolyser/app/tools && \
    rm /app/Immunolyser/app/tools/gibbscluster-2.0f.Linux.tar.gz

# Update GIBBS path in the gibbscluster file
RUN sed -i 's|setenv\s*GIBBS .*|setenv GIBBS /app/Immunolyser/app/tools/gibbscluster-2.0|' /app/Immunolyser/app/tools/gibbscluster-2.0/gibbscluster

# Comment out the line containing '$resdir .= "/$prefix";' in GibbsCluster-2.0e_SA.pl; command update for seq2logo and gibbs
RUN sed -i \
    -e 's|^\(\s*\$resdir .= "/\$prefix";\)|# \1  # Comment or remove this line|' \
    -e 's|^\(my \$barplot = "\$resdir/images/\$prefix.gibbs.KLDvsCluster.barplot.png";\)|my \$barplot = "\$resdir/images/gibbsKLDvsCluster.barplot.JPG";|' \
    -e '530s@.*@$cmd .= "$seq2logo -f $corefile -o $logofile -I 2 --format [JPEG] -b $wlc -C 2 -t $title \&>/dev/null; ";@' \
    /app/Immunolyser/app/tools/gibbscluster-2.0/GibbsCluster-2.0e_SA.pl

# Copy the netMHCpan tar.gz file to the container
COPY /tools/netMHCpan-4.1b.Linux.tar.gz /app/Immunolyser/app/tools/

# Uncompress and untar the netMHCpan package
RUN mkdir -p /app/Immunolyser/app/tools && \
    cat /app/Immunolyser/app/tools/netMHCpan-4.1b.Linux.tar.gz | gunzip | tar xvf - -C /app/Immunolyser/app/tools && \
    rm /app/Immunolyser/app/tools/netMHCpan-4.1b.Linux.tar.gz && \
    man -d /app/Immunolyser/app/tools/netMHCpan-4.1/netMHCpan.1 | compress > /app/Immunolyser/app/tools/netMHCpan-4.1/netMHCpan.Z && \
    mkdir -p /app/Immunolyser/app/tools/netMHCpan-4.1/tmp

# Download and untar the data.tar.gz file from the CBS website
RUN wget https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/data.tar.gz -P /app/Immunolyser/app/tools/netMHCpan-4.1 && \
    tar -xvf /app/Immunolyser/app/tools/netMHCpan-4.1/data.tar.gz -C /app/Immunolyser/app/tools/netMHCpan-4.1 && \
    rm /app/Immunolyser/app/tools/netMHCpan-4.1/data.tar.gz

# Copy the netMHCIIpan tar.gz file to the container
COPY /tools/netMHCIIpan-4.3e.Linux.tar.gz /app/Immunolyser/app/tools/

# Uncompress and untar the netMHCIIpan package
RUN mkdir -p /app/Immunolyser/app/tools && \
    tar -xvf /app/Immunolyser/app/tools/netMHCIIpan-4.3e.Linux.tar.gz -C /app/Immunolyser/app/tools && \
    rm /app/Immunolyser/app/tools/netMHCIIpan-4.3e.Linux.tar.gz && \
    man -d /app/Immunolyser/app/tools/netMHCIIpan-4.3/netMHCIIpan.1 | compress > /app/Immunolyser/app/tools/netMHCIIpan-4.3/netMHCIIpan.Z

# Update netMHCIIpan configuration to use the correct NMHOME path
RUN sed -i 's|setenv\s*NMHOME\s*/tools/src/netMHCIIpan-4.3|setenv NMHOME ${PWD}/app/tools/netMHCIIpan-4.3|' \
    /app/Immunolyser/app/tools/netMHCIIpan-4.3/netMHCIIpan

# Update netMHCpan configuration to use the correct NMHOME and TMPDIR paths
RUN sed -i \
    -e 's|setenv\s*NMHOME\s*/net/sund-nas.win.dtu.dk/storage/services/www/packages/netMHCpan/4.1/netMHCpan-4.1|setenv NMHOME ${PWD}/app/tools/netMHCpan-4.1|' \
    -e 's|setenv\s*TMPDIR\s*/tmp|setenv TMPDIR $NMHOME/tmp|' \
    /app/Immunolyser/app/tools/netMHCpan-4.1/netMHCpan

# Clone MixMHCpred repository
RUN git clone https://github.com/GfellerLab/MixMHCpred.git /app/Immunolyser/app/tools/MixMHCpred && \
    chmod +x /app/Immunolyser/app/tools/MixMHCpred/MixMHCpred

RUN wget https://github.com/GfellerLab/MixMHC2pred/releases/download/v2.0.2.2/MixMHC2pred-2.0.zip -P /app/Immunolyser/app/tools && \
    unzip -o /app/Immunolyser/app/tools/MixMHC2pred-2.0.zip -d /app/Immunolyser/app/tools/MixMHC2pred-2.0 && \
    rm /app/Immunolyser/app/tools/MixMHC2pred-2.0.zip

# Download Alleles_list_Mouse.txt and put it in PWMdef directory
RUN wget http://ec2-18-188-210-66.us-east-2.compute.amazonaws.com:4000/data/Alleles_lists/Alleles_list_Mouse.txt -P /app/Immunolyser/app/tools/MixMHC2pred-2.0/PWMdef

# Install mhcflurry
RUN pip install mhcflurry

# Fetch mhcflurry downloads
RUN mhcflurry-downloads fetch

# Create a virtual environment for Python 3
RUN python3 -m venv lenv

# Install dependencies for Python 2 and Python 3
RUN /bin/bash -c "source lenv/bin/activate && \
    pip install -r requirements_python2.txt && \
    pip install -r requirements_python3.txt"

# Install Celery, sqlacheny
RUN /bin/bash -c "pip install celery"
RUN /bin/bash -c "pip install SQLAlchemy==2.0.31"

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

# Create the results.sqlite file
RUN touch /results.sqlite

# Set the entrypoint script
ENTRYPOINT ["/entrypoint.sh"]
