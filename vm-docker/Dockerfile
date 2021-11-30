FROM 'google/cloud-sdk:slim'


#RUN apt-get update && apt-get install openjdk-8-jdk
# install packages 

RUN apt-get update && apt-get install -y --no-install-recommends \
    apt-utils \
    g++ \
    gawk \
    less \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libncurses5-dev \
    man-db \
    pkg-config \
    python3-venv \
    software-properties-common \
    unzip \
    wget \
    zlib1g-dev \
    && \
	# clean up apt cache
    rm -rf /var/lib/apt/lists/*

#install Java 8 for hail
RUN wget -qO - https://adoptopenjdk.jfrog.io/adoptopenjdk/api/gpg/key/public | apt-key add - \ 
    && add-apt-repository -y https://adoptopenjdk.jfrog.io/adoptopenjdk/deb/ \ 
    && apt-get update && apt-get install -y adoptopenjdk-8-hotspot

RUN mkdir tools 
WORKDIR /tools

# Install htslib
ENV SAMTOOLS_VERSION="1.11"
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd htslib-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean


# Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/${SAMTOOLS_VERSION}/bcftools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf bcftools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm bcftools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd bcftools-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd samtools-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean


# install vcftools
ENV VCFTOOLS_VERSION="0.1.16"
RUN wget https://github.com/vcftools/vcftools/releases/download/v${VCFTOOLS_VERSION}/vcftools-${VCFTOOLS_VERSION}.tar.gz \
    && tar xf vcftools-${VCFTOOLS_VERSION}.tar.gz \
    && rm vcftools-${VCFTOOLS_VERSION}.tar.gz \
	&& cd vcftools-${VCFTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean


# install picard
ENV PICARD_VERSION="2.25.0"
RUN wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar 

# install GATK
 ENV GATK_VERSION="4.2.0.0"
 RUN wget https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip \
    && unzip gatk-${GATK_VERSION}.zip \
    && rm gatk-${GATK_VERSION}.zip 


# Install python packages
RUN apt-get update \
	&& apt-get dist-upgrade -y \
	&& apt-get install -y --no-install-recommends\
    libc6-dev \
    libffi-dev \
    libgdbm-dev \
    liblapack-dev \
    liblapack3 \ 
    libncursesw5-dev \
    libopenblas-base \
    libopenblas-dev \
    libsqlite3-dev \
    libssl-dev \
    openssl \
    python-smbus \
    python3 python-pip 

# Upgrade pip to latest version
RUN python3 -m pip install --upgrade pip

# Install hail and other python libraries
ENV HAIL_VERSION="0.2.79"
RUN python3 --version
RUN python3 -m pip install \
    wheel \
    pypandoc \
    hail==${HAIL_VERSION} \
    scipy \
    numpy \
    pandas \
    matplotlib \
    seaborn \
    ipython \
    pybedtools \
    dill

# Install GCS Connector
RUN export SPARK_HOME=$(find_spark_home.py) && \
    curl https://storage.googleapis.com/hadoop-lib/gcs/gcs-connector-hadoop2-2.0.1.jar \
         >$SPARK_HOME/jars/gcs-connector-hadoop2-2.0.1.jar && \
    mkdir -p $SPARK_HOME/conf && \
    touch $SPARK_HOME/conf/spark-defaults.conf && \
    sed -i $SPARK_HOME/conf/spark-defaults.conf \
        -e 's:spark\.hadoop\.google\.cloud\.auth\.service\.account\.enable.*:spark.hadoop.google.cloud.auth.service.account.enable true:' \
        -e 's:spark\.hadoop\.google\.cloud\.auth\.service\.account\.json\.keyfile.*:spark\.hadoop\.google\.cloud\.auth\.service\.account\.json\.keyfile /gsa-key/key.json:'

# Install R
ENV R_VERSION="3.5.2-1"
RUN apt-get install -y \
    r-base=${R_VERSION} \
    r-base-dev=${R_VERSION} 

ENV BEDTOOLS_VERSION="2.29.0"
RUN wget https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz \
	&& tar xzf bedtools-${BEDTOOLS_VERSION}.tar.gz \
	&& rm bedtools-${BEDTOOLS_VERSION}.tar.gz \
	&& cd bedtools2 \
	&& make \
	&& make install \
	&& make clean

WORKDIR /home
