# Immunolyser 2.0 Local Version Setup Guide

## Prerequisites

1. **Docker Desktop**  
   You’ll need Docker Desktop installed on your machine. Docker is used to run applications in containers, ensuring compatibility across systems.  
   - **Download Docker Desktop**: [Docker Desktop Download](https://www.docker.com/products/docker-desktop)
   - **System Requirements**: Windows (10/11), macOS, or Linux. You will need administrator privileges to run Docker.

## Setup Instructions

### 1. **Download the Dockerfile**
   - [Download the Dockerfile here](https://github.com/prmunday/Immunolyser/blob/develop/Dockerfile) and save it to your local directory.

### 2. **Create Required Folders**
   In the same directory where you downloaded the Dockerfile, create the following folders:
   - `data_files`
   - `data_files_2`
   - `tools`

### 3. **Download and Place Tools**
   To make the Docker setup work, you'll need to download several tools and place them in the `tools` folder. Follow the steps below:

   - **Seq2Logo**:
     - [Download Seq2Logo 2.0](https://services.healthtech.dtu.dk/services/Seq2Logo-2.0/).
     - Click on the "Download" tab and select version 2.1.
     - Download the `.gz` file and save it in the `tools` folder, next to the Dockerfile.

   - **GibbsCluster 2.0f**:
     - [Download GibbsCluster 2.0f](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=gibbscluster&version=2.0&packageversion=2.0f&platform=Linux).
     - Save the `.gz` file in the `tools` folder.

   - **NetMHCpan 4.1b**:
     - [Download NetMHCpan 4.1b](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=netMHCpan&version=4.1&packageversion=4.1b&platform=Linux).
     - Save the `.gz` file in the `tools` folder.
     - Ensure the file is named `netMHCpan-4.1b.Linux.tar.gz` (remove any extra numbers added by the system).

   - **NetMHCIIpan 4.3e**:
     - [Download NetMHCIIpan 4.3e](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=netMHCIIpan&version=4.3&packageversion=4.3e&platform=Linux).
     - Save the `.gz` file in the `tools` folder as `netMHCIIpan-4.3e.Linux.tar.gz`.

### 4. **Build and Run Docker**
   After you've set up the tools, you can build and run the Docker container. Follow these steps:

   1. Open **Terminal** (macOS/Linux) or **Command Prompt** (Windows).
   2. Navigate to the directory where the Dockerfile is located.
   3. Run the following commands:
      ```sh
      docker-compose build
      docker-compose up
      ```
   - Docker will take **5 to 10 minutes** to build the container, depending on your system.
   - Make sure you have enough disk space for the build process.

### 5. **Check for Errors**
   If any errors occur during the Docker setup, check that all the files are in place and named correctly.

---

### Additional Notes:

- **Supported Systems**:  
  This guide supports the following operating systems:
  - **Windows** (10/11)
  - **macOS**
  - **Linux**

- **Troubleshooting**:  
  If you encounter issues, check the [Docker Troubleshooting Guide](https://docs.docker.com/get-docker/) for help with installation or configuration.

---

**Why Docker?**  
Using Docker ensures that all the tools and dependencies are correctly set up in a container, so you don't need to worry about installing them individually on your system. It guarantees compatibility and simplifies the process of setting up Immunolyser 2.0.

---
