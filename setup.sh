virtualenv venv/
source venv/bin/activate
pip install -r requirements_venv.txt
cp mpld3_fix/_display.py venv/lib/python2.7/site-packages/mpld3/
deactivate