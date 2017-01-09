## Mutagen Project

### Quick Start

Linux machine with Python 2.7

install pip 

#### make a virtual environment somewhere:
$ virtualenv web

It will create a directory web

#### Clone repository from Bitbucket to local machine
we will clone into web/mutagen subdirectory by running the following command
$ git clone git@bitbucket.org/agoncear/model_cancer_mutations.git

#### Activate virtual environment
$ source env/bin/activate


#### Install required packages
$ pip install -r mutagen/config/development/requirements.txt

### Running code
$ cd mutagen

$ python run.py

Open with a web browser http://127.0.0.1:5000

### Updating source code from the repository
$ cd mutagen

$ git pull

### Making changes and pushing them into the repository
cd mutagen

git add <changed filename>

git commit

git push

### Git manipulation with a GUI
Use SourceTree (free) application: http://sourcetreeapp.com
