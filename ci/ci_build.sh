PREFIX="/opt/conda"
mkdir -p "${PREFIX}/bin"
cp PlastEDMA/workflow/plastedma.py plastedma/resources/* "${PREFIX}/bin"
chmod +x /opt/conda/bin/plastedma.py