# How to set up AWS machine for assembly
---

If you are hoping to attempt a Trinity assembly, requirements for RAM = .5 * X million read pairs. For instance, to assemble 40 million paired-end reads using Trinity, you'll need a minimum of 20Gb of RAM. For BinPacker, you'll need substantially more, maybe as much as 2 * X million read pairs. 

These instructions work with a standard Ubuntu 16.04 (for instance, ami-2ef48339) machine available on AWS. Similar instructions should work for people on their own workstations, especially if you have `sudo` privileges. 


### Update Software and install things from apt-get

```
sudo apt-get update && sudo apt-get -y upgrade && sudo apt-get -y dist-upgrade

sudo apt-get -y install build-essential git python-pip python-numpy python-matplotlib

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
chmod +x filter.py
curl -LO https://s3.amazonaws.com/gen711/TruSeq3-PE.fa
PATH=$PATH:$(pwd)
```

### Install Perl Module
```
sudo cpan URI::Escape
```

### Install Ruby and LinuxBrew

```
cd
wget https://keybase.io/mpapis/key.asc
gpg --import key.asc
\curl -sSL https://get.rvm.io | bash -s stable --ruby
source /home/ubuntu/.rvm/scripts/rvm


sudo mkdir /home/linuxbrew
sudo chown $USER:$USER /home/linuxbrew
git clone https://github.com/Linuxbrew/brew.git /home/linuxbrew/.linuxbrew
echo 'export PATH="/home/linuxbrew/.linuxbrew/bin:$PATH"' >> ~/.profile
echo 'export MANPATH="/home/linuxbrew/.linuxbrew/share/man:$MANPATH"' >> ~/.profile
echo 'export INFOPATH="/home/linuxbrew/.linuxbrew/share/info:$INFOPATH"' >> ~/.profile
source ~/.profile
brew tap homebrew/science
brew update
brew doctor
```

### Install SolexaQA

```
curl -LO http://downloads.sourceforge.net/project/solexaqa/src/SolexaQA%2B%2B_v3.1.4.zip
unzip SolexaQA%2B%2B_v3.1.4.zip
cd Linux_x64
PATH=$PATH:$(pwd)
```


### Install Software (gcc skewer seqtk python jellyfish bfc rcorrector trinity BLAST)

```
brew install skewer seqtk python jellyfish bfc rcorrector hmmer infernal \
trinity --without-express vsearch salmon kallisto transdecoder last
```



### Install BLAST
```
cd
curl -LO ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.5.0+-x64-linux.tar.gz
tar -zxf ncbi-blast-2.5.0+-x64-linux.tar.gz
PATH=$PATH:/home/ubuntu/ncbi-blast-2.5.0+/bin
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

### Install TransFuse

```
gem install transfuse
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
git clone https://gitlab.com/ezlab/busco.git
cd busco
PATH=$PATH:$(pwd)
wget http://cegg.unige.ch/pub/BUSCO2/mammalia_odb9.tar.gz && tar -zxf mammalia_odb9.tar.gz
wget http://cegg.unige.ch/pub/BUSCO2/eukaryota_odb9.tar.gz && eukaryota_odb9.tar.gz
wget http://cegg.unige.ch/pub/BUSCO2/metazoa_odb9.tar.gz && metazoa_odb9.tar.gz
```




### Install dammit!

```
gem install crb-blast
pip install -U setuptools
pip install pandas
pip install dammit
sed -i 's/BUSCO_v1.1b1/BUSCO/' /home/linuxbrew/.linuxbrew/lib/python2.7/site-packages/dammit/dependencies.py
sed -i 's/BUSCO_v1.1b1/BUSCO/' /home/linuxbrew/.linuxbrew/lib/python2.7/site-packages/dammit/tasks.py
```




### Add all these things to the permanent path

```
echo PATH=$PATH >> ~/.profile
echo export LD_LIBRARY_PATH=/home/ubuntu/salmon-0.6.0/lib >> ~/.profile
echo source /home/ubuntu/.rvm/scripts/rvm >> ~/.profile
source ~/.profile
```
