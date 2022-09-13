#!/bin/bash

echo

realpath() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

ipoldir="$(realpath )"

ptype=$1
starg=$2
enarg=$3
raddir="$(realpath $4)"
tempdir="$(realpath $5)"
simdir=$6
doppdir=$7

mkdir -p $raddir $tempdir
if [ ! -z $simdir ]; then
    mkdir -p $simdir
fi
if [ ! -z $doppdir ]; then
    mkdir -p $doppdir
fi


st=$(echo $starg | sed 's/[^0-9]*//g')
stdate=${st:0:8}
sttime=${st:8:4}
en=$(echo $enarg | sed 's/[^0-9]*//g')
endate=${en:0:8}
entime=${en:8:4}
stt=${stdate}_${sttime}
edt=${endate}_${entime}

configdir=$ipoldir/configtxt/${stt}_${edt}
mkdir -p $configdir

echo Determining data agency and type...
echo
sleep 3

station=$(basename $raddir)
if [[ "$station" == "CASAG" ]]; then
    agency='cwr'
elif [[ "$station" == "NPOL" ]]; then
    agency='olympex'
else
    agency='nexrad'
fi
echo $agency

if [[ "$(ls $raddir/* | head -n 1 | xargs basename)" == "wrfout"* ]]; then
	data='wrf'
    mp=$(ls $raddir/* | head -n 1 | xargs basename | cut -d '_' -f2)
    inputfile=input_${data}_${mp}_${stt}_${edt}.txt
    configfile=config_${data}_${mp}_${stt}_${edt}.txt
else
	data='obs'
    inputfile=input_${data}_${stt}_${edt}.txt
    configfile=config_${data}_${stt}_${edt}.txt
fi
echo $data

echo
echo Selecting radar files for analysis in range $stt to $edt...
echo
sleep 3

tfile=input_${stt}_${edt}.txt

for filepath in $(ls $raddir/* | sort); do
    file=$(basename $filepath)
    if [[ $data == 'obs' ]]; then
        filedt=$(echo $file | cut -d '_' -f2)$(echo $file | cut -d '_' -f3)
    else
        filedt=$(echo $file | cut -d '_' -f4 | tr -d '-')$(echo $file | cut -d '_' -f5 | cut -d '.' -f1 | tr -d ':')
    fi
    if [ "$filedt" -ge "$(echo $stt | tr -d '_')00" ] && [ "$filedt" -lt "$(echo $edt | tr -d '_')00" ]; then 
        echo $filepath >> $configdir/$tfile
        echo $(basename $filepath)
    fi
    if [ "$filedt" -ge "$(echo $edt | tr -d '_')00" ]; then
        break
    fi
done

mv $configdir/$tfile $configdir/$inputfile

#latcen=$(ncdump $indir/$headfile | grep "latitude =" | cut -d '=' -f2 | cut -d ';' -f1 | xargs)
#loncen=$(ncdump $indir/$headfile | grep "longitude =" | cut -d '=' -f2 | cut -d ';' -f1 | xargs)
latcen=47.116943359375
loncen=-124.106666564941
#latcen=49.0164413452148
#loncen=-122.487358093262

echo
echo Selecting temperature files for analysis in range $stt to $edt...
echo
sleep 3

tfile=tmp_${stt}_${edt}.txt

for filepath in $(ls $tempdir/* | sort); do
    file=$(basename $filepath)
    filedt=$(echo $file | cut -d '_' -f4 | tr -d '-')$(echo $file | cut -d '_' -f5 | cut -d '.' -f1 | tr -d ':')
    if [ "$filedt" -ge "$(echo $stt | tr -d '_')00" ] && [ "$filedt" -le "$(echo $edt | tr -d '_')00" ]; then 
        echo $filepath >> $configdir/$tfile
        echo $(basename $filepath)
    fi
    if [ "$filedt" -gt "$(echo $edt | tr -d '_')00" ]; then
        break
    fi
done

if [[ "$(head -n 1 $configdir/$tfile | xargs basename)" == "wrfout"* ]]; then
	temp='wrf'
    mp=$(head -n 1 $configdir/$tfile | xargs basename | cut -d '_' -f2)
    tempfile=temp_${temp}_${mp}_${stt}_${edt}.txt
    snd_on='False'
    wrft_on='True'
else
	temp='snd'
    tempfile=temp_${temp}_${stt}_${edt}.txt
    snd_on='True'
    wrft_on='False'
fi
mv $configdir/$tfile $configdir/$tempfile

if [ -z $simdir ]; then
    fold=$data
else
    fold='obswrf'
fi

if [ -z $doppdir ]; then
    dd_on='False'
else
    dd_on='True'
fi

outfigdir=outputfig/${fold}_temp${temp}_${station}_${stt}_${edt}
outrrdir=$(cd $raddir/../../ && pwd)/radar_rainrates/$station
mkdir -p $outfigdir $outrrdir

echo
echo Creating config file for iPOLARRIS...
echo
sleep 3

template=$ipoldir/${agency}_${data}_config_template.txt
cp $template $configdir/$configfile

sed -i '' "s/^type ==.*/type == $data == # Type of input data: 'obs' OR 'wrf' (obs + simulated)/g" $configdir/$configfile
sed -i '' "s/.*mphys ==.*/mphys == $data == # Type of microphysics used in model: 'obs' OR '<scheme>' if type = 'wrf'/g" $configdir/$configfile
sed -i '' "s/.*ptype ==.*/ptype == '$ptype' == # Output figure file extenstion (i.e. png, jpg, mp4, ...)/g" $configdir/$configfile
sed -i '' "s/.*sdatetime ==.*/sdatetime == '$(echo $stt | tr '_' '-')' == # Start time of analysis of interest/g" $configdir/$configfile
sed -i '' "s/.*edatetime ==.*/edatetime == '$(echo $edt | tr '_' '-')' == # End time of analysis of interest/g" $configdir/$configfile
sed -i '' "s%.*rfiles ==.*%rfiles == '$configdir/$inputfile' == # Path to list of radar files to read in%g" $configdir/$configfile
sed -i '' "s%.*wfiles ==.*%wfiles == '$configdir/$tempfile' == # Path to list of WRF temperature files to read in%g" $configdir/$configfile
sed -i '' "s/.*exper ==.*/exper == $station == # Radar location/g" $configdir/$configfile
sed -i '' "s/lat ==  == #.*/lat == $latcen == # Latitude of the radar station/g" $configdir/$configfile
sed -i '' "s/lon ==  == #.*/lon == $loncen == # Longitude of the radar station/g" $configdir/$configfile
sed -i '' "s%.*image_dir ==.*%image_dir == '$outfigdir/' == # Output figure directory%g" $configdir/$configfile
sed -i '' "s%.*rr_dir ==.*%rr_dir == '$outrrdir/' == # Output rain rate netcdf directory%g" $configdir/$configfile
sed -i '' "s/.*dd_on ==.*/dd_on == $dd_on == # Doppler gridded velocity on/g" $configdir/$configfile
sed -i '' "s/.*snd_on ==.*/snd_on == $snd_on == # Sounding temperature on/g" $configdir/$configfile
sed -i '' "s/.*wrft_on ==.*/wrft_on == $wrft_on == # WRF temperature on/g" $configdir/$configfile

echo Running iPOLARRIS...
sleep 3

if [ ! -z $simdir]; then
    python run_ipolarris.py $configdir/$configfile $configdir/$configfile2
else
    python run_ipolarris.py $configdir/$configfile
fi
