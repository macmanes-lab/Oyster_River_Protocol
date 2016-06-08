# How to set up AWS machine for assembly
---

If you are hoping to attempt a Trinity assembly, requirements for RAM = .5 * X million read pairs. For instance, to assemble 40 million paired-end reads using Trinity, you'll need a minimum of 20Gb of RAM. For BinPacker, you'll need substantially more, maybe as much as 2 * X million read pairs. 

These instructions work with a standard Ubuntu 14.04 machine available on AWS. Similar instructions should work for people on their own workstations, especially if you have `sudo` privileges. 


### Update Software and install things from apt-get

```
sudo apt-get update && sudo apt-get -y upgrade

sudo apt-get -y install cmake sparsehash valgrind libboost-atomic1.55-dev libibnetdisc-dev gsl-bin bowtie \
      libgsl0-dev libgsl0ldbl libboost1.55-all-dev libboost1.55-dbg subversion tmux git curl parallel \
      libncurses5-dev samtools gcc make g++ python-dev unzip dh-autoreconf default-jre python-pip zlib1g-dev \
      hmmer libhdf5-dev r-base pkg-config libpng12-dev libfreetype6-dev python-sklearn build-essential \
      libsm6 libxrender1 libfontconfig1 liburi-escape-xs-perl emboss python-biopython liburi-perl infernal python-numpy

```

### Format and Mount hard drive (if needed)

```
sudo mkfs -t ext4 /dev/xvdf
sudo mount /dev/xvdf /mnt
sudo chown -R ubuntu:ubuntu /mnt
```

### Install Adapter seqs and a few utility scripts

```
cd && mkdir share && cd share
curl -LO https://raw.githubusercontent.com/macmanes-lab/general/master/filter.py
chmox +x filter.py
curl -LO https://s3.amazonaws.com/gen711/TruSeq3-PE.fa
PATH=$PATH:$(pwd)
```

### Install Perl Module
```
sudo cpan URI::Escape
```

### Install Ruby 2.x

```
cd
wget https://keybase.io/mpapis/key.asc
gpg --import key.asc
\curl -sSL https://get.rvm.io | bash -s stable --ruby
```

### Install SolexaQA

```
curl -LO http://downloads.sourceforge.net/project/solexaqa/src/SolexaQA%2B%2B_v3.1.4.zip
unzip SolexaQA%2B%2B_v3.1.4.zip
cd Linux_x64
PATH=$PATH:$(pwd)
```


### Install Skewer

```
cd $HOME
git clone https://github.com/relipmoc/skewer.git
cd skewer
make
PATH=$PATH:$(pwd)
curl -LO https://s3.amazonaws.com/gen711/TruSeq3-PE.fa
```

### Install bfc
You don't need this is you're using Rcorrector

```
cd $HOME
git clone https://github.com/lh3/bfc.git
cd bfc
make
PATH=$PATH:$(pwd)
```

### Install RCorrector
You don't need this is you're using bfc

```
cd
git clone https://github.com/mourisl/Rcorrector.git
cd Rcorrector
make
PATH=$PATH:$(pwd)
```

### Install BLAST
```
cd
curl -LO ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.4.0+-x64-linux.tar.gz
tar -zxf ncbi-blast-2.4.0+-x64-linux.tar.gz
PATH=$PATH:/home/ubuntu/ncbi-blast-2.4.0+/bin
```

### Install Trinity

```
cd
git clone https://github.com/trinityrnaseq/trinityrnaseq.git
cd trinityrnaseq
make -j6
PATH=$PATH:$(pwd)
```

### Install BinPacker
```
cd
git clone https://github.com/macmanes-lab/BinPacker.git
cd BinPacker

#### Change install.sh (this will not be necessary if following instructions on AWS)
#### change line to ./configure --with-boost=/home/ubuntu/boost/
#### save file

sh install.sh
```


### Install Vsearch

```
cd
git clone https://github.com/torognes/vsearch.git
cd vsearch
sh autogen.sh
./configure
make -j4
PATH=$PATH:$(pwd)/bin
```


### Install TransFuse

```
gem install transfuse
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
curl -LO https://github.com/COMBINE-lab/salmon/archive/v0.6.0.tar.gz
tar -zxf v0.6.0.tar.gz && rm v0.6.0.tar.gz
cd salmon-0.6.0/
mkdir build
cd build
cmake ..
make -j6
sudo make all install
export LD_LIBRARY_PATH=/home/ubuntu/salmon-0.6.0/lib
```

### Install Transrate

```
cd
curl -LO https://bintray.com/artifact/download/blahah/generic/transrate-1.0.3-linux-x86_64.tar.gz
tar -zxf transrate-1.0.3-linux-x86_64.tar.gz
PATH=$PATH:/home/ubuntu/transrate-1.0.3-linux-x86_64:/home/ubuntu/transrate-1.0.3-linux-x86_64/bin
```

### Install BUSCO

```
cd
curl -LO http://busco.ezlab.org/files/BUSCO_v1.2.tar.gz
tar -zxf BUSCO_v1.2.tar.gz
cd BUSCO_v1.2
chmod +x BUSCO_v1.2.py
PATH=$PATH:$(pwd)
curl -LO http://busco.ezlab.org/files/metazoa_buscos.tar.gz
curl -LO http://busco.ezlab.org/files/vertebrata_buscos.tar.gz
tar -zxf vertebrata_buscos.tar.gz
tar -zxf metazoa_buscos.tar.gz
```


### Install Transdecoder

```
cd
curl -LO https://github.com/TransDecoder/TransDecoder/archive/v3.0.0.tar.gz
tar -xvzf v3.0.0.tar.gz
cd TransDecoder-3.0.0
make -j6
export PATH=$PATH:$HOME/TransDecoder-3.0.0
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
gem install crb-blast
sudo pip install -U setuptools
sudo pip install numpy --upgrade
sudo pip install matplotlib --upgrade
sudo pip install dammit
```




### Add all these things to the permanent path

```
echo PATH=$PATH >> ~/.profile
echo export LD_LIBRARY_PATH=/home/ubuntu/salmon-0.6.0/lib >> ~/.profile
echo source /home/ubuntu/.rvm/scripts/rvm >> ~/.profile
source ~/.profile
```
