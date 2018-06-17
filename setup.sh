virtualenv venv/
source venv/bin/activate
pip install -r requirements_venv.txt
nodeenv -p --requirements=requirements_venv_nodejs.txt --jobs=4 --force venv
deactivate
