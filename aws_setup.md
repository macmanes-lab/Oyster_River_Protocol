# How to set up AWS machine for assembly



### Update Software and install things from apt-get

```
sudo apt-get update && sudo apt-get -y upgrade

sudo apt-get -y install cmake sparsehash valgrind libboost-atomic1.55-dev libibnetdisc-dev ruby-full gsl-bin \
      libgsl0-dev libgsl0ldbl libboost1.55-all-dev libboost1.55-dbg subversion tmux git curl bowtie \
      libncurses5-dev samtools gcc make g++ python-dev unzip dh-autoreconf default-jre python-pip zlib1g-dev \
      hmmer cmake libhdf5-dev

```

### Format and Mount hard drive (if needed)

```
sudo mkfs -t ext4 /dev/xvdf
sudo mount /dev/xvdf /mnt
sudo chown -R ubuntu:ubuntu /mnt
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
cd $HOME
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
cd $HOME
curl -LO https://github.com/COMBINE-lab/salmon/archive/v0.5.1.tar.gz
tar -zxf v0.5.1.tar.gz
cd salmon-0.5.1/
mkdir build
cd build
cmake -DCMAKE_C_COMPILER=$(which gcc-4.9) -DCMAKE_CXX_COMPILER=$(which g++-4.9) ..
make -j6
sudo make all install
```

### Install Transrate

```
cd
curl -LO https://bintray.com/artifact/download/blahah/generic/transrate-1.0.1-linux-x86_64.tar.gz
tar -zxf transrate-1.0.1-linux-x86_64.tar.gz
PATH=$PATH:/home/ubuntu/transrate-1.0.1-linux-x86_64
```

### Install BUSCO

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
PATH:$PATH:$(pwd)
```
