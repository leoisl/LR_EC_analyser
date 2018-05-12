virtualenv venv/
source venv/bin/activate
pip install -r requirements_venv.txt
cp mpld3_fix/_display.py venv/lib/python2.7/site-packages/mpld3/
nodeenv -p --requirements=requirements_venv_nodejs.txt --jobs=4 --force venv
deactivate
