PREFIX="/opt/conda"
mkdir -p "${PREFIX}/bin"
cp PlastEDMA/plastedma.py "${PREFIX}/bin/"
cp -r PlastEDMA/workflow "${PREFIX}/bin/"
# cp -r PlastEDMA/workflow/scripts "${PREFIX}/bin/"
cp -r PlastEDMA/resources "${PREFIX}/bin"
cp -r PlastEDMA/results "${PREFIX}/bin"
cp -r PlastEDMA/config "${PREFIX}/bin"
chmod +x /opt/conda/bin/plastedma.py
