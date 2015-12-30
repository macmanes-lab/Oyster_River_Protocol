# How to set up AWS machine for assembly
---

If you are hoping to attempt a Trinity assembly, requirements for RAM = .5 * X million read pairs. For instance, to assemble 40 million paired-end reads using Trinity, you'll need a minimum of 20Gb of RAM.

These instructions work with a standard Ubuntu 14.04 machine available on AWS. Similar instructions should work for people on their own workstations, especially if you have `sudo` privileges. 


### Update Software and install things from apt-get

```
sudo apt-get update && sudo apt-get -y upgrade

sudo apt-get -y install cmake sparsehash valgrind libboost-atomic1.55-dev libibnetdisc-dev ruby-full gsl-bin \
      libgsl0-dev libgsl0ldbl libboost1.55-all-dev libboost1.55-dbg subversion tmux git curl bowtie \
      libncurses5-dev samtools gcc make g++ python-dev unzip dh-autoreconf default-jre python-pip zlib1g-dev \
      hmmer cmake libhdf5-dev r-base pkg-config libpng12-dev libfreetype6-dev python-sklearn build-essential \
      libsm6 libxrender1 libfontconfig1 liburi-escape-xs-perl emboss liburi-perl infernal python-pip python-dev python-numpy

```

### Format and Mount hard drive (if needed)

```
sudo mkfs -t ext4 /dev/xvdf
sudo mount /dev/xvdf /mnt
sudo chown -R ubuntu:ubuntu /mnt
```

### Install SolexaQA


```
curl -LO http://downloads.sourceforge.net/project/solexaqa/src/SolexaQA%2B%2B_v3.1.4.zip
unzip SolexaQA%2B%2B_v3.1.4.zip
cd Linux_x64
PATH=$PATH:$(pwd)
```

### Install Transdecoder

```
cd
curl -LO https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz
tar -xvzf 2.0.1.tar.gz
cd TransDecoder-2.0.1; make
export PATH=$PATH:$HOME/TransDecoder-2.0.1
```

### Install LAST

```
cd
curl -LO http://last.cbrc.jp/last-658.zip
unzip last-658.zip
cd last-658
make
export PATH=$PATH:$HOME/last-658/src
```

### Install dammit!

```
sudo gem install crb-blast
sudo pip install -U setuptools
sudo pip install numpy --upgrade
sudo pip install matplotlib --upgrade
sudo pip install dammit
```


### Install Perl Module
```
sudo cpan URI::Escape
```

### Install seqtk

```
cd $HOME
git clone https://github.com/lh3/seqtk.git
cd seqtk
make -j4
PATH=$PATH:$(pwd)
```

### Install bwa

```
cd $HOME
git clone https://github.com/lh3/bwa.git
cd bwa
make -j4
PATH=$PATH:$(pwd)
```

### Install khmer

```
sudo easy_install -U setuptools
sudo pip install khmer
```

### Install Kallisto

```
cd
git clone https://github.com/pachterlab/kallisto.git
cd kallisto
mkdir build
cd build
cmake ..
make
sudo make install
```


### Install Salmon

```
cd
curl -LO https://github.com/COMBINE-lab/salmon/archive/v0.5.1.tar.gz
tar -zxf v0.5.1.tar.gz
cd salmon-0.5.1/
mkdir build
cd build
cmake -DCMAKE_C_COMPILER=$(which gcc-4.9) -DCMAKE_CXX_COMPILER=$(which g++-4.9) ..
make -j6
sudo make all install
export LD_LIBRARY_PATH=/home/ubuntu/salmon-0.5.1/lib
```

### Install Transrate

```
cd
curl -LO https://bintray.com/artifact/download/blahah/generic/transrate-1.0.1-linux-x86_64.tar.gz
tar -zxf transrate-1.0.1-linux-x86_64.tar.gz
PATH=$PATH:/home/ubuntu/transrate-1.0.1-linux-x86_64
```

### Install BUSCO

```
cd
curl -LO http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz
tar -zxf BUSCO_v1.1b1.tar.gz
cd BUSCO_v1.1b1
chmod +x BUSCO_v1.1b1.py
PATH=$PATH:$(pwd)
curl -LO http://busco.ezlab.org/files/metazoa_buscos.tar.gz
curl -LO http://busco.ezlab.org/files/vertebrata_buscos.tar.gz
tar -zxf vertebrata_buscos.tar.gz
tar -zxf metazoa_buscos.tar.gz


cd
curl -LO http://bioinf.uni-greifswald.de/augustus/binaries/old/augustus-3.0.2.tar.gz
tar -zxf augustus-3.0.2.tar.gz
cd augustus-3.0.2/
make
PATH=$PATH:$(pwd)/bin:$(pwd)/scripts
export AUGUSTUS_CONFIG_PATH=/home/ubuntu/augustus-3.0.2/config/
```

### Install BLAST


```
cd
curl -LO ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.3.0+-x64-linux.tar.gz
tar -zxf ncbi-blast-2.3.0+-x64-linux.tar.gz
PATH=$PATH:/home/ubuntu/ncbi-blast-2.3.0+/bin
```

### Install Trinity

```
cd
git clone https://github.com/trinityrnaseq/trinityrnaseq.git
cd trinityrnaseq
make -j6
PATH=$PATH:$(pwd)
```

### Install bfc

```
cd $HOME
git clone https://github.com/lh3/bfc.git
cd bfc
make
PATH=$PATH:$(pwd)
```

###Install RCorrector

```
cd
git clone https://github.com/mourisl/Rcorrector.git
cd Rcorrector
make
PATH=$PATH:$(pwd)
```

### Add all these things to the permanent path

```
echo PATH=$PATH >> ~/.profile
source ~/.profile
```
