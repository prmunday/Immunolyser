# Use the official Python 3.12 slim image as a parent image
FROM python:3.12-slim

# Set the working directory in the container
WORKDIR /app

# Install git and other required dependencies
RUN apt-get update && apt-get install -y git

# Clone the repository
RUN git clone https://github.com/prmunday/Immunolyser /app/Immunolyser

# Change to the repository directory
WORKDIR /app/Immunolyser

# Checkout the develop branch
RUN git fetch --all && git checkout feature/docker_setup

# Create a virtual environment
RUN python3 -m venv lenv

# Activate the virtual environment and install dependencies
RUN /bin/bash -c "source lenv/bin/activate && \
    pip install -r requirements_python2.txt && \
    pip3 install -r requirements_python3.txt && \
    pip3 install celery"

# Run the hotfix script
RUN /bin/bash -c "source lenv/bin/activate && python hotfix_package_files.py"

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
