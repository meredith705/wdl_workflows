# WDL Workflows

WDL workflows for useful tools. 

# WDL Cromwell Dependencies 
Download cromwell tools to run wdl worflows locally:
```
wget https://github.com/broadinstitute/cromwell/releases/download/71/cromwell-71.jar
wget https://github.com/broadinstitute/cromwell/releases/download/71/womtool-71.jar
```
Install jdk
```
sudo apt-get update
sudo apt install openjdk-17-jre-headless
```
Install docker ( directions from here https://docs.docker.com/engine/install/ubuntu/ )
```
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io
```

# Installation
Clone this git repo to obtain the wdl workflows:
```
git clone https://github.com/meredith705/wdl_workflows.git
```
