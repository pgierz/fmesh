#!/bin/bash

curl https://nextcloud.awi.de/s/yNyyKYiifbgbAKT/download/RTopo-2.0.1_30sec_bos_fix_lowres_D3.nc > ./topo/RTopo-2.0.1_30sec_bos_fix_lowres_D3.nc
curl https://nextcloud.awi.de/s/yDPwoX8iFcLyyf9/download/Archive.zip > ./coastlines/coastlines.zip
cd coastlines
unzip coastlines.zip

