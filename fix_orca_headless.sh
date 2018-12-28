mv venv/bin/orca venv/bin/orca_core
printf '#!/bin/bash\nxvfb-run -a venv/bin/orca_core "$@"' > venv/bin/orca
chmod +x venv/bin/orca