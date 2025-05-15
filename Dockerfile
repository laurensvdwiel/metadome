FROM --platform=linux/amd64 python:3.13.3-bullseye

# Create directory to install required externals
RUN mkdir -p /usr/externals

# Install packages required for blast
Run apt-get update
RUN apt-get install libidn11

# Install BLAST+
RUN wget http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz
RUN tar -zxvf ncbi-blast-2.6.0+-x64-linux.tar.gz
RUN mv ncbi-blast-2.6.0+ /usr/externals/blast

# Install ClustalW
RUN wget http://www.clustal.org/download/2.1/clustalw-2.1-linux-x86_64-libcppstatic.tar.gz
RUN tar -zxvf clustalw-2.1-linux-x86_64-libcppstatic.tar.gz
RUN mv clustalw-2.1-linux-x86_64-libcppstatic /usr/externals/clustalw

# Install HMMER Tools
RUN wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz
RUN tar -zxvf hmmer-3.1b2-linux-intel-x86_64.tar.gz
RUN mv hmmer-3.1b2-linux-intel-x86_64 /usr/externals/hmmer

# Create directory for app
RUN mkdir -p /usr/src/app
WORKDIR /usr/src/app

# Configure Python environment
COPY requirements.txt /usr/src/app/
RUN pip install --no-cache-dir -r requirements.txt
COPY . /usr/src/app