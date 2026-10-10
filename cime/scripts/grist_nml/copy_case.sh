Grid_level=$1
# G8
if [ -z "$Grid_level" ]; then
    echo "Error: Grid_level is not set."
    echo "Example: ./copy_case.sh G8"
    exit 1
fi

Dir=`pwd`
nml_path=/home/export/online1/mdt00/shisuan/swgbcm/licom/data/GB25/Grist-Data/${Grid_level}/nml
table_path=/home/export/online1/mdt00/shisuan/swgbcm/licom/data/GB25/Grist-Data/table_data


echo 'Copy Grid_level      : '${Grid_level}
echo 'Copy NML file from   : '${nml_path}
echo 'Link Table File from : '${table_path}
echo 'Direct to            : '${Dir} 

sleep 5
echo 'Copying ......'
ls ${nml_path}
cp ${nml_path}/*  ${Dir}

echo 'Linking ......'
ls ${table_path}
ln -sf ${table_path}/*  ${Dir}

