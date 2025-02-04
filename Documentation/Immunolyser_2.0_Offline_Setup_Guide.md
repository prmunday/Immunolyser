# Immunolyser Setup Guide

## Prerequisites

1. **Docker**: Ensure Docker Desktop is installed on your machine. You can download it from Docker Desktop. You will need machine administrator privileges to run it.

## Steps

1. **Download the Dockerfile**:
   - Download the Dockerfile from here.

2. **Create Required Folders**:
   - Create the following folders next to the Dockerfile:
     - `data_files`
     - `data_files_2`
     - `tools`

3. **Download and Place Tools**:

   - **Seq2Logo**:
     - Visit Seq2Logo 2.0.
     - Click on the "Download" tab and select version 2.1.
     - Download the `.gz` file and place it in the `tools` folder.

   - **GibbsCluster 2.0f**:
     - Download from GibbsCluster 2.0f.
     - Place the downloaded `.gz` file in the `tools` folder.

   - **NetMHCpan 4.1b**:
     - Download from NetMHCpan 4.1b.
     - Ensure the file name is `netMHCpan-4.1b.Linux.tar.gz` (remove any extra numbers if added).
     - Place the file in the `tools` folder.

   - **NetMHCIIpan 4.3e**:
     - Download from NetMHCIIpan 4.3e.
     - Place the file `netMHCIIpan-4.3e.Linux.tar.gz` in the `tools` folder.

4. **Build and Run Docker**:
   - Open a terminal and navigate to the directory containing the Dockerfile.
   - Run the following commands:
     ```sh
     docker-compose build
     docker-compose up
     ```
   - Note: The Docker build process can take 5 to 10 minutes depending on your system, and it will require sufficient space on your system.