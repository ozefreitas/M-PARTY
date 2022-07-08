PREFIX="/opt/conda"
mkdir -p "${PREFIX}/bin"
mkdir -p "${PREFIX}/bin/workflow"
cp PlastEDMA/workflow/plastedma.py "${PREFIX}/bin/workflow"
cp -r PlastEDMA/workflow/scripts "${PREFIX}/bin/workflow"
cp -r PlastEDMA/resources "${PREFIX}/bin"
cp -r PlastEDMA/results "${PREFIX}/bin"
cp -r PlastEDMA/config "${PREFIX}/bin"
chmod +x /opt/conda/bin/plastedma.py