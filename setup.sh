set -eux
virtualenv venv/
source venv/bin/activate
pip install -r requirements_venv.txt
deactivate