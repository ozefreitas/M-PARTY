PREFIX="/opt/conda"
mkdir -p "${PREFIX}/bin"
cp M-PARTY/m_party.py "${PREFIX}/bin/"
cp -r M-PARTY/workflow "${PREFIX}/bin/"
# cp -r M-PARTY/workflow/scripts "${PREFIX}/bin/"
cp -r M-PARTY/resources "${PREFIX}/bin"
# cp -r M-PARTY/results "${PREFIX}/bin"
cp -r M-PARTY/config "${PREFIX}/bin"
chmod +x /opt/conda/bin/m_party.py